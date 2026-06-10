// // The libMesh Finite Element Library.
// // Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
// 
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// 
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// // Lesser General Public License for more details.
// 
// // You should have received a copy of the GNU Lesser General Public
// // License along with this library; if not, write to the Free Software
// // Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// stdlib include ---------------------------------------------------->
#include <cmath>
#include <array>
#include <iomanip>


// Femus lib includes ----------------->
#include "MGFE_L2_FE_template.h"
#include "MGFE_L2_lagrange_shape_1D.h"



// template <> double  MGFE<3,LAGRANGE_XFEM>::Rec_Lin_DPhi(int nPhi, double point[], int dimension, int DirDer);
// template <> double  MGFE<3,LAGRANGE_XFEM>::Rec_Quad_DPhi(int nPhi, double point[], int dimension, int DirDer);
// // Anonymous namespace for functions shared by LAGRANGE_XFEM and
// // L2_LAGRANGE implementations. Implementations appear at the bottom
// // of this file.
// namespace
// {
// // FE class members
template <>
 MGFE<3,LAGRANGE_XFEM>::MGFE (const FEType & fet) :
   MGFEbase (3,fet),
   // _dim(3),             ///< space dimension
   //  _order(0),         ///< pol approx (2=quad  1=lin 0=k)
   //  _NoShape(0),      ///< number of shapes
   //  _GeomType(0),  ///< geom type 0=hex/quad  1=tetra/tri
   //  _FamType(0), 
   last_side(INVALID_ELEM),
   last_edge(INVALID_ELEM)
 {
   // Sanity check.  Make sure the
   // Family specified in the template instantiation
   // matches the one in the FEType object
   assert (LAGRANGE_XFEM == this->get_family());
 }
 
 template <>
 MGFE<3,LAGRANGE_XFEM>::~MGFE(){}
 
 
  template <>  const int MGFE<3,LAGRANGE_XFEM>::get_dim()       {return _dim;}  ///< Dimension (3 or 2) used in the 3D FEM (1< _dim < MDIM_FEM)
 template <>  const int MGFE<3,LAGRANGE_XFEM>::get_order()     {return _order;}///< Order of the shape functions 0=cpnst 1=linear 2=quad
 template <>  const int MGFE<3,LAGRANGE_XFEM>::get_NoShape()   {return _NoShape;}  ///< Number of shape functions of the fem
 template <>  const int MGFE<3,LAGRANGE_XFEM>::get_NoElNodes() {return _NoElNodes;}  ///< Number of nodes in the reference element (always ref quad)
 template <>  const int MGFE<3,LAGRANGE_XFEM>::get_GeomType()  {return _GeomType;}   ///< Geom type 1=HEXAHEDRAL 0=TETRAHEDRAL
 template <>  const int MGFE<3,LAGRANGE_XFEM>::get_FamType()   {return _FamType;}    ///< Family  (geom 1: Lagrangian =0  RT =1)
  
   // Gaussian points -----------------------------------------------------------------------------
 template <>  const int MGFE<3,LAGRANGE_XFEM>::get_NoGauss(){ return _NoGauss;}  ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
  template <> const double  MGFE<3,LAGRANGE_XFEM>::get_xg(int i){          return  _xg[i];}          ///< gaussian points
  template <> const double MGFE<3,LAGRANGE_XFEM>::get_weight(int i){      return _weight[i];}            ///< Weight
  template <> const double &MGFE<3,LAGRANGE_XFEM>::get_phi(){     return  *phi;}           ///< Shape functions
  template <> const double & MGFE<3,LAGRANGE_XFEM>::get_dphidxez_map(){return  *_dphidxez_map;}      ///< Shape derivative functions in gaussian points
  template <> const double & MGFE<3,LAGRANGE_XFEM>::get_dphidxx_map(){ return  *_dphidxx_map;}       ///< Second order shape derivatives in Gaussian points
  
   // Dofs points ----------------------------------------------------------------------------------
 template <>  const double &  MGFE<3,LAGRANGE_XFEM>::get_xnodes(){return   *_xnodes;}              ///< xnode (hex or tetra) in (xi,eta,seta)
 template <>  const int &  MGFE<3,LAGRANGE_XFEM>::get_fem_conn(){ return *_fem_conn;}              ///< dof node coords (respect to _xnodes)
  // template <> const int     MGFE<3,LAGRANGE_XFEM>::get_deg(){return _deg;}                  ///< Polinomial degree
  
   // Nodes (NDOF_FEM) -----------------------------------------------------------------------------
 template <>  const double & MGFE<3,LAGRANGE_XFEM>::get_phi_map_nodes(){return   *_phi_map_nodes;   }          ///< Shape functions
 template <>  const double & MGFE<3,LAGRANGE_XFEM>::get_dphidxez_map_nodes(){return   *_dphidxez_map_nodes;} 
 
template <> int MGFE<3,LAGRANGE_XFEM>::GetFamilyType(int elem_dof, int dim) {
  int FamilyType;

    if(elem_dof == 7 || elem_dof == 10) {      FamilyType = 0;    }
    else if(elem_dof == 8 || elem_dof == 27) {      FamilyType = 1;    }
    else   printf("MGFE1::GetFamilyType Unknown FamilyType for element with dimension %d and %d number of dof \n",
      dim, elem_dof);

  return FamilyType;
}
 
 
template <> void MGFE<3,LAGRANGE_XFEM>::print_rec_xdmf(
const std::string& namefile,  // file <-
const std::string& namefem,  // file <-
const std::string& namelem  // file <-
) {  //=================================================================================================
  // order xdmf only for 27 not regular -----------------------------------

//   int hex_shift=0;   
//   if(_dim==1) {hex_shift= 6;} //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//   if(_dim==2) {hex_shift= 12;}//  QUAD9  _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)]


  int *ord=new int[27]; for(int i=0; i<27; i++) ord[i]=i;
  if(_NoElNodes==27) {
    for(int i=12; i<=15; i++) ord[i]=i+4; 
    for(int i=16; i<=19; i++) ord[i]=i-4; 
    ord[20]=24;
    ord[21]=22; 
    ord[22]=21; 
    ord[23]=23; 
    ord[24]=20;
  }

  // file -------------------------------------------------------------------
  std::ofstream out(namefile.c_str());
  if(!out) {std::cout << " Gauss Outfile MGFE1::write_c not opened." << std::endl; exit(3); }
  out << "<?xml version=\"1.0\" ?> "<<"\n"
      <<"<Xdmf>      "<<"\n"
      << "<Domain> "<<"\n";
// ------------------------------------------------------------------
// Grid EL (NDOF_FEM)
// ---------------------------------------------------------------------------------------
 out     << "  <Grid Name= \"Geom element (surface view)\" > "<<"\n"
      << "    <Topology Type=\""<< namelem.c_str()<<"\"  Dimensions=\""<<1 <<"\" > </Topology> "<<"\n"
      <<"           <Geometry GeometryType=\"XYZ\"> "<<"\n"
      <<"               <DataItem Format=\"XML\" Dimensions=\""<< "1  " <<   _NoElNodes  << "  3 "<< "\"> "<<"\n";
  for(int is=0; is< _NoElNodes; is++) {
    int is1=ord[is];
    out << "  ";
    for(int idim=0; idim<_dim; idim++) out << _xnodes[is1+ _NoElNodes*idim] << " ";
    for(int idim=_dim; idim<3; idim++) out << "   0.   ";
    out << "  \n";
  }
  out <<"               </DataItem> "<<"\n"
      <<"           </Geometry> "<<"\n";
  //  Attribute  shape on nodes ----------------------------------------------------  
  for(int is=0; is<_NoShape; is++) { //
//   for(int dir=0;dir< _dim;dir++){
    out<< " <Attribute Name=\"" << "phi"<<is<< "  "  <<"\" Center=\"Node\"> "<<"\n"
       << "      <DataItem Format=\"XML\" Dimensions=\""<< " " <<  _NoElNodes << " " <<"\" >"<<"\n";
    for(int i=0; i<  _NoElNodes; i++) { // points
      int i1=ord[i];
      out <<        _phi_map_nodes[(is) *  _NoElNodes + i1] << " ";
    }
    out <<"     </DataItem> "<<"\n"
        <<"   </Attribute>   "<<"\n";
//   }// _dir
  } //is ----------------------------------------------------------------------------
  
 // Attribute derivatives on nodes---------------------------------------------------------
  for(int is=0; is<_NoShape; is++) { //shape
   for(int dir=0;dir< _dim;dir++){
    out<< " <Attribute Name=\"" << "dphi"<<is<< "/d"<<dir  <<"\" Center=\"Node\"> "<<"\n"
       << "      <DataItem Format=\"XML\" Dimensions=\""<< " " <<  _NoElNodes << " " <<"\" >"<<"\n";
    for(int i=0; i<  _NoElNodes; i++) { // points
      out << _dphidxez_map_nodes[(is+dir*_NoShape)* _NoElNodes+i] << " "; 
    }
    out <<"     </DataItem> "<<"\n"
        <<"   </Attribute>   "<<"\n";
    }// _dir
  } //is ----------------------------------------------------------------  
    
  
  
  out<<"       </Grid> "<<"\n";
// -----------------------------------------------------------------------------------------------
//    Grid dofs (_NoShape)
//   ---------------------------------------------------------------------
out  << "  <Grid Name= \"Dof points (surface view)\">       "<<"\n"
     <<"    <Topology Type= \""<< namefem.c_str()<<"\"  Dimensions= \"1\"> </Topology>  "<<"\n"
     <<"            <Geometry GeometryType=\"XYZ\">"<<"\n"
     <<"               <DataItem Format=\"XML\" Dimensions=\""
     << "1  " <<  _NoShape  << "  3 "<< "\"> "<<"\n";
  for(int is=0; is<_NoShape; is++) {
    int is1=ord[_fem_conn[is]];
    out << "  ";
    for(int idim=0; idim<_dim; idim++) out << _xnodes[is1+idim* _NoElNodes] << " ";
    for(int idim=_dim; idim<3; idim++) out << "   0.   ";
    out << "  \n";
  }
  out <<"               </DataItem> "<<"\n"
      <<"           </Geometry> "<<"\n"
      <<"       </Grid> "<<"\n";
// -----------------------------------------------------------------------
  //  Grid gauss
// -----------------------------------------------------------------------
  out   << "  <Grid Name= \"Gauss points (Point view)\"> "<<"\n"
        <<"    <Topology Type=\"" << "Polyvertex"  <<"\"  Dimensions=\""<< _NoGauss<<"\"> "<<"\n"
        <<"           </Topology> "<<"\n"
        <<"           <Geometry GeometryType=\"XYZ\"> "<<"\n"
        <<"               <DataItem Format=\"XML\" Dimensions=\""<< "1  " <<  _NoGauss  << "  3 "<< "\"> "<<"\n";
  for(int is=0; is<_NoGauss; is++) {
    out << "  ";
    for(int idim=0; idim<_dim; idim++) out << _xg[is+idim*_NoGauss] << " ";
    for(int idim=_dim; idim<3; idim++) out << "   0.   ";
    out << "  \n";
  }
  out <<"               </DataItem> "<<"\n"
      <<"           </Geometry>     "<<"\n";
  //  Attribute  shape on Gauss-----------------------------------------------------------
  for(int is=0; is<_NoShape; is++) {

//   for(int dir=0;dir< _dim;dir++){
    out<< " <Attribute Name=\"" << "phi"<<is<< "  "  <<"\" Center=\"Node\"> "<<"\n"
       << "      <DataItem Format=\"XML\" Dimensions=\""<< " " << _NoGauss << " " <<"\" >"<<"\n";
    for(int i=0; i< _NoGauss; i++) { // points
      out <<        phi[(is) *  _NoGauss + i] << " ";
      
      
    }
    out <<"     </DataItem> "<<"\n"
        <<"   </Attribute>   "<<"\n";
//   }// _dir
  } //is ----------------------------------------------------------------------
  
  // Attribute derivatives on Gauss---------------------------------------------------------
  for(int is=0; is<_NoShape; is++) { //shape
   for(int dir=0;dir< _dim;dir++){
    out<< " <Attribute Name=\"" << "dphi"<<is<< "/d"<<dir  <<"\" Center=\"Node\"> "<<"\n"
       << "      <DataItem Format=\"XML\" Dimensions=\""<< " " << _NoGauss << " " <<"\" >"<<"\n";
    for(int i=0; i< _NoGauss; i++) { // points
      out <<       _dphidxez_map[(is+dir*_NoShape)*_NoGauss+i] << " ";
    }
    out <<"     </DataItem> "<<"\n"
        <<"   </Attribute>   "<<"\n";
    }// _dir
  } //is ----------------------------------------------------------------
  
  
  out   <<"       </Grid>             "<<"\n"
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
        <<"   </Domain> "<<"\n"
        <<"</Xdmf>     "<<"\n";


  delete[] ord;
   
  
  return;
} 
 
 //=================================================================================================
/// This function writes shape and derivative values at the gaussian points
 template <> void MGFE<3,LAGRANGE_XFEM>::write_c(
std::ostream& out  // file <-
) {  //=================================================================================================

  if(!out) {
    std::cout << " Gauss Outfile MGFE1::write_c not opened." << std::endl;
    exit(3);
  }

  // heading
  out << "gpoints " << _NoGauss << std::endl << std::endl;
  double sum_g; double sum_gd[3];
  if(_FamType ==0) {
    // max_nodes (decreasing level order)
    for(int k = 0; k < _NoGauss; k++) {
      out <<"g= "<< k << " w= " << std::setprecision(20) << _weight[k] << std::endl;
      sum_g =0.; sum_gd[0] =0.; sum_gd[1] =0.; sum_gd[2] =0.;
      for(int s = 0; s < _NoShape; s++) {
        sum_g += phi[s * _NoGauss + k];
        out << std::setprecision(20) << phi[s * _NoGauss + k] << " ";
        for(int idim = 0; idim < _dim; idim++) {
          out << std::setprecision(20)
              << _dphidxez_map[(idim * _NoShape + s) * _NoGauss + k] << "  ";
          sum_gd[idim] +=_dphidxez_map[(idim * _NoShape + s) * _NoGauss + k];
        }
        out << std::endl;
      }
      out << std::setprecision(20) <<" test  sum phi ="<< sum_g<< " sum  dphi/dx= ";
      for(int idim = 0; idim < _dim; idim++) out << std::setprecision(20) <<"     "<< sum_gd[idim] <<"  ;  ";
      out << std::endl;
    }
  }

  if(_FamType ==1) {
    // max_nodes (decreasing level order)
    for(int k = 0; k < _NoGauss; k++) {
      out <<" Vectorial RT0 g= "<< k << " w= " << std::setprecision(20) << _weight[k] << std::endl;
      for(int s = 0; s < _NoShape; s++) {
// phi[N1x(g1),N1x(g2),N1x(g3),.... N1y(g1),N1y(g2),N1y(g3),....N1z(g1),N1z(g2),N1z(g3),.......
        //          N2x(g1),N2x(g2),N2x(g3),.... N2y(g1),N2y(g2),N2y(g3),....N2z(g1),N2z(g2),N2z(g3),.......
        //          .......................................................................................]
        out << std::setprecision(20)<< "(";
        for(int idim = 0; idim < _dim; idim++) {
          out << std::setprecision(20)<<  phi[s * _dim*_NoGauss + k+ idim*_NoGauss] << " , ";
        }
        out << std::setprecision(20)<< ")";
        out << std::endl;
      }
      out << std::endl;
    }
  }
#ifdef PRINT_INFO
  std::cout << " MGFE1::write_c:  " << " fem with " << _NoGauss
            << " gaussian points and " << _NoShape << " shape functions \n"
            << std::endl;
#endif
  return;
}  

   // ==================================================================
 template <> void MGFE<3,LAGRANGE_XFEM>::write(
const std::string& name  // file <-
) {                           // ================================
  std::ofstream in(name.c_str());
  this->write_c(in);
}  
  // ==================================================================
 template <>
void MGFE<3,LAGRANGE_XFEM>::attach_quadrature_rule (QBase *q)
 {
assert(q);
   this->_qrule = q;
   // make sure we don't cache results from a previous quadrature rule
   this->_elem_type = INVALID_ELEM;
   return;
 }
 template <>
 unsigned int MGFE<3,LAGRANGE_XFEM>::n_shape_functions () const {
    return MGFE<3,LAGRANGE_XFEM>::n_dofs (this->_elem_type, this->_fe_type.order + this->_p_level);
  }
// using namespace libMesh;
// 
// template <FEFamily T>
// double fe_lagrange_3D_shape(const ElemType,
//                           const Order order,
//                           const unsigned int i,
//                           const Point & p);
// 
 template <FEFamily T>
 double fe_lagrange_3D_shape_deriv(const ElemType type,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const double p[]);
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
template <FEFamily T>
 double fe_lagrange_3D_shape_second_deriv(const ElemType type,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const double  p[]);
// 
// #endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
// } // anonymous namespace
// 
// namespace libMesh
// {
// 
// 
// // TODO: If optimizations for LAGRANGE_XFEM work well we should do
// // L2_LAGRANGE too...
// LIBMESH_DEFAULT_VECTORIZED_FE(3,L2_LAGRANGE)
// 
// 
// template<>
// void MGFE<3,LAGRANGE_XFEM>::all_shapes
//   (const Elem * elem,
//    const Order o,
//    const std::vector<Point> & p,
//    std::vector<std::vector<OutputShape>> & v,
//    const bool add_p_level)
// {
//   const ElemType type = elem->type();
// 
//   // Just loop on the harder-to-optimize cases
//   if (type != HEX8 && type != HEX27)
//     {
//       MGFE<3,LAGRANGE_XFEM>::default_all_shapes
//         (elem,o,p,v,add_p_level);
//       return;
//     }
// 
// #if LIBMESH_DIM == 3
// 
//   const unsigned int n_sf = v.size();
// 
//   switch (o)
//     {
//       // linear Lagrange shape functions
//     case FIRST:
//       {
//         switch (type)
//           {
//             // trilinear hexahedral shape functions
//           case HEX8:
//           case HEX20:
//           case HEX27:
//             {
//               libmesh_assert_less_equal (n_sf, 8);
// 
//               //                                0  1  2  3  4  5  6  7
//               static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
//               static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
//               static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};
// 
//               for (auto qp : index_range(p))
//                 {
//                   const Point & q_point = p[qp];
//                   // Compute hex shape functions as a tensor-product
//                   const double xi   = q_point(0);
//                   const double eta  = q_point(1);
//                   const double zeta = q_point(2);
// 
//                   // one_d_shapes[dim][i] = phi_i(p(dim))
//                   double one_d_shapes[3][2] = {
//                     {fe_lagrange_1D_linear_shape(0, xi),
//                      fe_lagrange_1D_linear_shape(1, xi)},
//                     {fe_lagrange_1D_linear_shape(0, eta),
//                      fe_lagrange_1D_linear_shape(1, eta)},
//                     {fe_lagrange_1D_linear_shape(0, zeta),
//                      fe_lagrange_1D_linear_shape(1, zeta)}};
// 
//                   for (unsigned int i : make_range(n_sf))
//                     v[i][qp] = one_d_shapes[0][i0[i]] *
//                                one_d_shapes[1][i1[i]] *
//                                one_d_shapes[2][i2[i]];
//                 }
//               return;
//             }
// 
//           default:
//             libmesh_error(); // How did we get here?
//           }
//       }
// 
// 
//       // quadratic Lagrange shape functions
//     case SECOND:
//       {
//         switch (type)
//           {
//             // triquadratic hexahedral shape functions
//           case HEX8:
// // TODO: refactor to optimize this
// //            assert(T == L2_LAGRANGE);std::cerr<<
// //                               "High order on first order elements only supported for L2 families";
//             (void(0));
//           case HEX27:
//             {
//               libmesh_assert_less_equal (n_sf, 27);
// 
//               // The only way to make any sense of this
//               // is to look at the mgflo/mg2/mgf documentation
//               // and make the cut-out cube!
//               //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
//               static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
//               static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
//               static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};
// 
//               for (auto qp : index_range(p))
//                 {
//                   const Point & q_point = p[qp];
//                   // Compute hex shape functions as a tensor-product
//                   const double xi   = q_point(0);
//                   const double eta  = q_point(1);
//                   const double zeta = q_point(2);
// 
//                   // linear_shapes[dim][i] = phi_i(p(dim))
//                   double one_d_shapes[3][3] = {
//                     {fe_lagrange_1D_quadratic_shape(0, xi),
//                      fe_lagrange_1D_quadratic_shape(1, xi),
//                      fe_lagrange_1D_quadratic_shape(2, xi)},
//                     {fe_lagrange_1D_quadratic_shape(0, eta),
//                      fe_lagrange_1D_quadratic_shape(1, eta),
//                      fe_lagrange_1D_quadratic_shape(2, eta)},
//                     {fe_lagrange_1D_quadratic_shape(0, zeta),
//                      fe_lagrange_1D_quadratic_shape(1, zeta),
//                      fe_lagrange_1D_quadratic_shape(2, zeta)}};
// 
//                   for (unsigned int i : make_range(n_sf))
//                     v[i][qp] = one_d_shapes[0][i0[i]] *
//                                one_d_shapes[1][i1[i]] *
//                                one_d_shapes[2][i2[i]];
//                 }
//               return;
//             }
// 
//           default:
//             libmesh_error(); // How did we get here?
//           }
//       }
// 
//       // unsupported order
//     default:
//       std::cerr<<"ERROR: Unsupported 3D FE order on HEX!: " << o);
//     }
// #else // LIBMESH_DIM != 3
//   libmesh_ignore(elem, o, p, v, add_p_level);
//   libmesh_not_implemented();
// #endif // LIBMESH_DIM == 3
// }
// 
// template<>
// void MGFE<3,LAGRANGE_XFEM>::shapes
//   (const Elem * elem,
//    const Order o,
//    const unsigned int i,
//    const std::vector<Point> & p,
//    std::vector<OutputShape> & v,
//    const bool add_p_level)
// {
//   MGFE<3,LAGRANGE_XFEM>::default_shapes
//     (elem,o,i,p,v,add_p_level);
// }
// 
// template<>
// void MGFE<3,LAGRANGE_XFEM>::shape_derivs
//   (const Elem * elem,
//    const Order o,
//    const unsigned int i,
//    const unsigned int j,
//    const std::vector<Point> & p,
//    std::vector<OutputShape> & v,
//    const bool add_p_level)
// {
//   MGFE<3,LAGRANGE_XFEM>::default_shape_derivs
//     (elem,o,i,j,p,v,add_p_level);
// }
// 
// template<>
// void MGFE<3,LAGRANGE_XFEM>::all_shape_derivs
//   (const Elem * elem,
//    const Order o,
//    const std::vector<Point> & p,
//    std::vector<std::vector<OutputShape>> * comps[3],
//    const bool add_p_level)
// {
//   const ElemType type = elem->type();
// 
//   // Just loop on the harder-to-optimize cases
//   if (type != HEX8 && type != HEX27)
//     {
//       MGFE<3,LAGRANGE_XFEM>::default_all_shape_derivs
//         (elem,o,p,comps,add_p_level);
//       return;
//     }
// 
// #if LIBMESH_DIM == 3
// 
//   libmesh_assert(comps[0]);
//   libmesh_assert(comps[1]);
//   libmesh_assert(comps[2]);
//   const unsigned int n_sf = comps[0]->size();
// 
//   switch (o)
//     {
//       // linear Lagrange shape functions
//     case FIRST:
//       {
//         switch (type)
//           {
//             // trilinear hexahedral shape functions
//           case HEX8:
//           case HEX20:
//           case HEX27:
//             {
//               libmesh_assert_equal_to (n_sf, 8);
// 
//               //                                0  1  2  3  4  5  6  7
//               static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
//               static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
//               static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};
// 
//               for (auto qp : index_range(p))
//                 {
//                   const Point & q_point = p[qp];
//                   // Compute hex shape functions as a tensor-product
//                   const double xi   = q_point(0);
//                   const double eta  = q_point(1);
//                   const double zeta = q_point(2);
// 
//                   // one_d_shapes[dim][i] = phi_i(p(dim))
//                   double one_d_shapes[3][2] = {
//                     {fe_lagrange_1D_linear_shape(0, xi),
//                      fe_lagrange_1D_linear_shape(1, xi)},
//                     {fe_lagrange_1D_linear_shape(0, eta),
//                      fe_lagrange_1D_linear_shape(1, eta)},
//                     {fe_lagrange_1D_linear_shape(0, zeta),
//                      fe_lagrange_1D_linear_shape(1, zeta)}};
// 
//                   // one_d_derivs[dim][i] = dphi_i/dxi(p(dim))
//                   double one_d_derivs[3][2] = {
//                     {fe_lagrange_1D_linear_shape_deriv(0, 0, xi),
//                      fe_lagrange_1D_linear_shape_deriv(1, 0, xi)},
//                     {fe_lagrange_1D_linear_shape_deriv(0, 0, eta),
//                      fe_lagrange_1D_linear_shape_deriv(1, 0, eta)},
//                     {fe_lagrange_1D_linear_shape_deriv(0, 0, zeta),
//                      fe_lagrange_1D_linear_shape_deriv(1, 0, zeta)}};
// 
//                     for (unsigned int i : make_range(n_sf))
//                       {
//                         (*comps[0])[i][qp] = one_d_derivs[0][i0[i]] *
//                                              one_d_shapes[1][i1[i]] *
//                                              one_d_shapes[2][i2[i]];
//                         (*comps[1])[i][qp] = one_d_shapes[0][i0[i]] *
//                                              one_d_derivs[1][i1[i]] *
//                                              one_d_shapes[2][i2[i]];
//                         (*comps[2])[i][qp] = one_d_shapes[0][i0[i]] *
//                                              one_d_shapes[1][i1[i]] *
//                                              one_d_derivs[2][i2[i]];
//                       }
//                 }
//               return;
//             }
// 
//           default:
//             libmesh_error(); // How did we get here?
//           }
//       }
// 
// 
//       // quadratic Lagrange shape functions
//     case SECOND:
//       {
//         switch (type)
//           {
//             // triquadratic hexahedral shape functions
//           case HEX8:
// // TODO: refactor to optimize this
// //            assert(T == L2_LAGRANGE);std::cerr<<
// //                               "High order on first order elements only supported for L2 families";
//             (void(0));
//           case HEX27:
//             {
//               libmesh_assert_less_equal (n_sf, 27);
// 
//               // The only way to make any sense of this
//               // is to look at the mgflo/mg2/mgf documentation
//               // and make the cut-out cube!
//               //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
//               static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
//               static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
//               static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};
// 
//               for (auto qp : index_range(p))
//                 {
//                   const Point & q_point = p[qp];
//                   // Compute hex shape functions as a tensor-product
//                   const double xi   = q_point(0);
//                   const double eta  = q_point(1);
//                   const double zeta = q_point(2);
// 
//                   // one_d_shapes[dim][i] = phi_i(p(dim))
//                   double one_d_shapes[3][3] = {
//                     {fe_lagrange_1D_quadratic_shape(0, xi),
//                      fe_lagrange_1D_quadratic_shape(1, xi),
//                      fe_lagrange_1D_quadratic_shape(2, xi)},
//                     {fe_lagrange_1D_quadratic_shape(0, eta),
//                      fe_lagrange_1D_quadratic_shape(1, eta),
//                      fe_lagrange_1D_quadratic_shape(2, eta)},
//                     {fe_lagrange_1D_quadratic_shape(0, zeta),
//                      fe_lagrange_1D_quadratic_shape(1, zeta),
//                      fe_lagrange_1D_quadratic_shape(2, zeta)}};
// 
//                   // one_d_derivs[dim][i] = dphi_i/dxi(p(dim))
//                   double one_d_derivs[3][3] = {
//                     {fe_lagrange_1D_quadratic_shape_deriv(0, 0, xi),
//                      fe_lagrange_1D_quadratic_shape_deriv(1, 0, xi),
//                      fe_lagrange_1D_quadratic_shape_deriv(2, 0, xi)},
//                     {fe_lagrange_1D_quadratic_shape_deriv(0, 0, eta),
//                      fe_lagrange_1D_quadratic_shape_deriv(1, 0, eta),
//                      fe_lagrange_1D_quadratic_shape_deriv(2, 0, eta)},
//                     {fe_lagrange_1D_quadratic_shape_deriv(0, 0, zeta),
//                      fe_lagrange_1D_quadratic_shape_deriv(1, 0, zeta),
//                      fe_lagrange_1D_quadratic_shape_deriv(2, 0, zeta)}};
// 
//                     for (unsigned int i : make_range(n_sf))
//                       {
//                         (*comps[0])[i][qp] = one_d_derivs[0][i0[i]] *
//                                              one_d_shapes[1][i1[i]] *
//                                              one_d_shapes[2][i2[i]];
//                         (*comps[1])[i][qp] = one_d_shapes[0][i0[i]] *
//                                              one_d_derivs[1][i1[i]] *
//                                              one_d_shapes[2][i2[i]];
//                         (*comps[2])[i][qp] = one_d_shapes[0][i0[i]] *
//                                              one_d_shapes[1][i1[i]] *
//                                              one_d_derivs[2][i2[i]];
//                       }
//                 }
//               return;
//             }
// 
//           default:
//             libmesh_error(); // How did we get here?
//           }
//       }
// 
//       // unsupported order
//     default:
//       std::cerr<<"ERROR: Unsupported 3D FE order on HEX!: " << o);
//     }
// #else // LIBMESH_DIM != 3
//   libmesh_ignore(elem, o, p, v, add_p_level);
//   libmesh_not_implemented();
// #endif // LIBMESH_DIM == 3
// }
// 
// 
// 
// 

// 
// 
// 
// template <>
// double MGFE<3,L2_LAGRANGE>::shape(const ElemType type,
//                               const Order order,
//                               const unsigned int i,
//                               const Point & p)
// {
//   return fe_lagrange_3D_shape<L2_LAGRANGE>(type, order, i, p);
// }
// 
// 
// 
// template <>
// double MGFE<3,LAGRANGE_XFEM>::shape(const Elem * elem,
//                            const Order order,
//                            const unsigned int i,
//                            const Point & p,
//                            const bool add_p_level)
// {
//   libmesh_assert(elem);
// 
//   // call the orientation-independent shape functions
//   return fe_lagrange_3D_shape<LAGRANGE_XFEM>(elem->type(), order + add_p_level*elem->p_level(), i, p);
// }
// 
// 
// 
// template <>
// double MGFE<3,L2_LAGRANGE>::shape(const Elem * elem,
//                               const Order order,
//                               const unsigned int i,
//                               const Point & p,
//                               const bool add_p_level)
// {
//   libmesh_assert(elem);
// 
//   // call the orientation-independent shape functions
//   return fe_lagrange_3D_shape<L2_LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), i, p);
// }
// 
// 
// 
// template <>
// double MGFE<3,LAGRANGE_XFEM>::shape(const FEType fet,
//                            const Elem * elem,
//                            const unsigned int i,
//                            const Point & p,
//                            const bool add_p_level)
// {
//   libmesh_assert(elem);
//   return fe_lagrange_3D_shape<LAGRANGE_XFEM>(elem->type(), fet.order + add_p_level*elem->p_level(), i, p);
// }
// 
// 
// 
// template <>
// double MGFE<3,L2_LAGRANGE>::shape(const FEType fet,
//                               const Elem * elem,
//                               const unsigned int i,
//                               const Point & p,
//                               const bool add_p_level)
// {
//   libmesh_assert(elem);
//   return fe_lagrange_3D_shape<L2_LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), i, p);
// }
// 
 template <>
 double MGFE<3,LAGRANGE_XFEM>::shape_deriv(const ElemType type,
                                  const Order order,
                                  const unsigned int i,
                                  const unsigned int j,
                                    const double p[])
                                  // const Point & p)
 {
   return fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(type, order, i, j, p);
 }
// 
// 
// 
// template <>
// double MGFE<3,L2_LAGRANGE>::shape_deriv(const ElemType type,
//                                     const Order order,
//                                     const unsigned int i,
//                                     const unsigned int j,
//                                     const Point & p)
// {
//   return fe_lagrange_3D_shape_deriv<L2_LAGRANGE>(type, order, i, j, p);
// }
// 
// 
// 
// template <>
// double MGFE<3,LAGRANGE_XFEM>::shape_deriv(const Elem * elem,
//                                  const Order order,
//                                  const unsigned int i,
//                                  const unsigned int j,
//                                  const Point & p,
//                                  const bool add_p_level)
// {
//   libmesh_assert(elem);
// 
//   // call the orientation-independent shape function derivatives
//   return fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(elem->type(), order + add_p_level*elem->p_level(), i, j, p);
// }
// 
// 
// template <>
// double MGFE<3,L2_LAGRANGE>::shape_deriv(const Elem * elem,
//                                     const Order order,
//                                     const unsigned int i,
//                                     const unsigned int j,
//                                     const Point & p,
//                                     const bool add_p_level)
// {
//   libmesh_assert(elem);
// 
//   // call the orientation-independent shape function derivatives
//   return fe_lagrange_3D_shape_deriv<L2_LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), i, j, p);
// }
// 
// 
// template <>
// double MGFE<3,LAGRANGE_XFEM>::shape_deriv(const FEType fet,
//                                  const Elem * elem,
//                                  const unsigned int i,
//                                  const unsigned int j,
//                                  const Point & p,
//                                  const bool add_p_level)
// {
//   libmesh_assert(elem);
//   return fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(elem->type(), fet.order + add_p_level*elem->p_level(), i, j, p);
// }
// 
// 
// template <>
// double MGFE<3,L2_LAGRANGE>::shape_deriv(const FEType fet,
//                                     const Elem * elem,
//                                     const unsigned int i,
//                                     const unsigned int j,
//                                     const Point & p,
//                                     const bool add_p_level)
// {
//   libmesh_assert(elem);
//   return fe_lagrange_3D_shape_deriv<L2_LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), i, j, p);
// }
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
 template <>
 double MGFE<3,LAGRANGE_XFEM>::shape_second_deriv(const ElemType type,
                                         const Order order,
                                         const unsigned int i,
                                         const unsigned int j,
                                          const double  p[])
                                         // const Point & p)
 {
   return fe_lagrange_3D_shape_second_deriv<LAGRANGE_XFEM>(type, order, i, j, p);
 }
// 
// 
// 
// template <>
// double MGFE<3,L2_LAGRANGE>::shape_second_deriv(const ElemType type,
//                                            const Order order,
//                                            const unsigned int i,
//                                            const unsigned int j,
//                                            const Point & p)
// {
//   return fe_lagrange_3D_shape_second_deriv<L2_LAGRANGE>(type, order, i, j, p);
// }
// 
// 
// 
// template <>
// double MGFE<3,LAGRANGE_XFEM>::shape_second_deriv(const Elem * elem,
//                                         const Order order,
//                                         const unsigned int i,
//                                         const unsigned int j,
//                                         const Point & p,
//                                         const bool add_p_level)
// {
//   libmesh_assert(elem);
// 
//   // call the orientation-independent shape function derivatives
//   return fe_lagrange_3D_shape_second_deriv<LAGRANGE_XFEM>
//     (elem->type(), order + add_p_level*elem->p_level(), i, j, p);
// }
// 
// 
// 
// template <>
// double MGFE<3,L2_LAGRANGE>::shape_second_deriv(const Elem * elem,
//                                            const Order order,
//                                            const unsigned int i,
//                                            const unsigned int j,
//                                            const Point & p,
//                                            const bool add_p_level)
// {
//   libmesh_assert(elem);
// 
//   // call the orientation-independent shape function derivatives
//   return fe_lagrange_3D_shape_second_deriv<L2_LAGRANGE>
//     (elem->type(), order + add_p_level*elem->p_level(), i, j, p);
// }
// 
// 
// template <>
// double MGFE<3,LAGRANGE_XFEM>::shape_second_deriv(const FEType fet,
//                                         const Elem * elem,
//                                         const unsigned int i,
//                                         const unsigned int j,
//                                         const Point & p,
//                                         const bool add_p_level)
// {
//   libmesh_assert(elem);
//   return fe_lagrange_3D_shape_second_deriv<LAGRANGE_XFEM>
//     (elem->type(), fet.order + add_p_level*elem->p_level(), i, j, p);
// }
// 
// 
// 
// template <>
// double MGFE<3,L2_LAGRANGE>::shape_second_deriv(const FEType fet,
//                                            const Elem * elem,
//                                            const unsigned int i,
//                                            const unsigned int j,
//                                            const Point & p,
//                                            const bool add_p_level)
// {
//   libmesh_assert(elem);
//   return fe_lagrange_3D_shape_second_deriv<L2_LAGRANGE>
//     (elem->type(), fet.order + add_p_level*elem->p_level(), i, j, p);
// }
// 
// 
// #endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
// } // namespace libMesh
// 
// 
// 
// namespace
// {
// using namespace libMesh;
// =====================================================================
template <FEFamily T>
double fe_lagrange_3D_shape(
  const ElemType type,
  const Order order,
  const unsigned int i,
  const double  p[]
){ // ==================================================================
// #if LIBMESH_DIM == 3
  switch (order)  {
    case CONSTANT:  return 1.; // ccccccccccccccccccccccccccccccccccccccccc 
      // linear Lagrange shape functions
    case FIRST: { //ffffffffffffffffffffffffffffffffffffffffffffffffffffffff>
        switch (type) {
            // trilinear hexahedral shape functions
          case HEX8: case HEX20:  case HEX27: { assert(i< 8);
              // Compute hex shape functions as a tensor-product
              const double xi   = p[0];  const double eta  = p[1];  const double zeta = p[2];
              //                                0  1  2  3  4  5  6  7
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};
              return (fe_lagrange_1D_linear_shape(i0[i], xi)*
                      fe_lagrange_1D_linear_shape(i1[i], eta)*
                      fe_lagrange_1D_linear_shape(i2[i], zeta));
            }
            // linear tetrahedral shape functions
          case TET4:  case TET10:  case TET14: {
              assert (i< 4);
              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const double zeta1 = p[0]; const double zeta2 = p[1]; const double zeta3 = p[2];
              const double zeta0 = 1. - zeta1 - zeta2 - zeta3;
              switch(i) {
                case 0: return zeta0;
                case 1: return zeta1;
                case 2: return zeta2;
                case 3: return zeta3;
                default:  std::cerr<<"Invalid i = " << i;
                }
            }
            // linear prism shape functions
          case PRISM6: case PRISM15:  case PRISM18: case PRISM20:  case PRISM21:  {  assert (i< 6);
              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;// Point p2d(p[0],p[1]);
              double p1d = p[2];
              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};
              return (MGFE<2,LAGRANGE_XFEM>::shape(TRI3,  FIRST, i1[i], p2d)*
                      fe_lagrange_1D_linear_shape(i0[i], p1d));
            }
            // linear pyramid shape functions -----------------------------------------
          case PYRAMID5:   case PYRAMID13:   case PYRAMID14:  case PYRAMID18:  {   assert (i< 5);
              const double xi   = p[0];  const double eta  = p[1];  const double zeta = p[2];
              const double eps  = 1.e-35;
              switch(i) {
                case 0:  return .25*(zeta + xi - 1.)*(zeta + eta - 1.)/((1. - zeta) + eps);
                case 1:  return .25*(zeta - xi - 1.)*(zeta + eta - 1.)/((1. - zeta) + eps);
                case 2:  return .25*(zeta - xi - 1.)*(zeta - eta - 1.)/((1. - zeta) + eps);
                case 3:  return .25*(zeta + xi - 1.)*(zeta - eta - 1.)/((1. - zeta) + eps);
                case 4:  return zeta;
                default:  std::cerr<<"Invalid i = " << i;
                }
            }
          default:  std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
          }
      } // <--ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


      // quadratic Lagrange shape functions sssssssssssssssssssssssssssssssssssssssssssssssssss->
    case SECOND:   {
        switch (type) {
            // serendipity hexahedral quadratic shape functions
          case HEX20:  {  assert (i< 20);
              const double xi   = p[0];  const double eta  = p[1];  const double zeta = p[2];
              // these functions are defined for (x,y,z) in [0,1]^3
              // so transform the locations
              const double x = .5*(xi   + 1.); const double y = .5*(eta  + 1.);   const double z = .5*(zeta + 1.);
              switch (i)  {
                case 0: return (1. - x)*(1. - y)*(1. - z)*(1. - 2.*x - 2.*y - 2.*z);
                case 1: return x*(1. - y)*(1. - z)*(2.*x - 2.*y - 2.*z - 1.);
                case 2: return x*y*(1. - z)*(2.*x + 2.*y - 2.*z - 3.);
                case 3: return (1. - x)*y*(1. - z)*(2.*y - 2.*x - 2.*z - 1.);
                case 4: return (1. - x)*(1. - y)*z*(2.*z - 2.*x - 2.*y - 1.);
                case 5: return x*(1. - y)*z*(2.*x - 2.*y + 2.*z - 3.);
                case 6: return x*y*z*(2.*x + 2.*y + 2.*z - 5.);
                case 7: return (1. - x)*y*z*(2.*y - 2.*x + 2.*z - 3.);
                case 8: return 4.*x*(1. - x)*(1. - y)*(1. - z);
                case 9: return 4.*x*y*(1. - y)*(1. - z);
                case 10: return 4.*x*(1. - x)*y*(1. - z);
                case 11: return 4.*(1. - x)*y*(1. - y)*(1. - z);
                case 12: return 4.*(1. - x)*(1. - y)*z*(1. - z);
                case 13: return 4.*x*(1. - y)*z*(1. - z);
                case 14: return 4.*x*y*z*(1. - z);
                case 15: return 4.*(1. - x)*y*z*(1. - z);
                case 16: return 4.*x*(1. - x)*(1. - y)*z;
                case 17: return 4.*x*y*(1. - y)*z;
                case 18: return 4.*x*(1. - x)*y*z;
                case 19: return 4.*(1. - x)*y*(1. - y)*z;
                default:  std::cerr<<"Invalid i = " << i;
                }
            }
            // triquadratic hexahedral shape functions
          case HEX8:  assert(T == L2_LAGRANGE); std::cerr<<  "High order on first order elements only supported for L2 families";
            (void(0));
          case HEX27: {   assert (i< 27);
              // Compute hex shape functions as a tensor-product
              const double xi   = p[0];  const double eta  = p[1];  const double zeta = p[2];
              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};

              return (fe_lagrange_1D_quadratic_shape(i0[i], xi)*
                      fe_lagrange_1D_quadratic_shape(i1[i], eta)*
                      fe_lagrange_1D_quadratic_shape(i2[i], zeta));
            }

            // quadratic tetrahedral shape functions
          case TET4:  assert(T == L2_LAGRANGE); std::cerr<<  "High order on first order elements only supported for L2 families";
            (void(0));
          case TET10:  assert (i< 10);
            (void(0));
          case TET14:  {  assert (i< 14);
              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const double zeta1 = p[0];  const double zeta2 = p[1];   const double zeta3 = p[2];
              const double zeta0 = 1. - zeta1 - zeta2 - zeta3;
              switch(i){
                case 0: return zeta0*(2.*zeta0 - 1.);
                case 1: return zeta1*(2.*zeta1 - 1.);
                case 2: return zeta2*(2.*zeta2 - 1.);
                case 3: return zeta3*(2.*zeta3 - 1.);
                case 4: return 4.*zeta0*zeta1;
                case 5: return 4.*zeta1*zeta2;
                case 6: return 4.*zeta2*zeta0;
                case 7: return 4.*zeta0*zeta3;
                case 8: return 4.*zeta1*zeta3;
                case 9: return 4.*zeta2*zeta3;
                default:  std::cerr<<"Invalid i = " << i;
                }
            }
            // "serendipity" prism
          case PRISM15: {  assert (i< 15);
              const double xi   = p[0]; const double eta  = p[1];  const double zeta = p[2];
              switch(i)  {
                case 0:  return (1. - zeta)*(xi + eta - 1.)*(xi + eta + 0.5*zeta);
                case 1:  return (1. - zeta)*xi*(xi - 1. - 0.5*zeta);
                case 2: // phi1 with xi <- eta
                  return (1. - zeta)*eta*(eta - 1. - 0.5*zeta);
                case 3: // phi0 with zeta <- (-zeta)
                  return (1. + zeta)*(xi + eta - 1.)*(xi + eta - 0.5*zeta);
                case 4: // phi1 with zeta <- (-zeta)
                  return (1. + zeta)*xi*(xi - 1. + 0.5*zeta);
                case 5: // phi4 with xi <- eta
                  return (1. + zeta)*eta*(eta - 1. + 0.5*zeta);
                case 6:  return 2.*(1. - zeta)*xi*(1. - xi - eta);
                case 7:  return 2.*(1. - zeta)*xi*eta;
                case 8:  return 2.*(1. - zeta)*eta*(1. - xi - eta);
                case 9:  return (1. - zeta)*(1. + zeta)*(1. - xi - eta);
                case 10: return (1. - zeta)*(1. + zeta)*xi;
                case 11: // phi10 with xi <-> eta
                  return (1. - zeta)*(1. + zeta)*eta;
                case 12: // phi6 with zeta <- (-zeta)
                  return 2.*(1. + zeta)*xi*(1. - xi - eta);
                case 13: // phi7 with zeta <- (-zeta)
                  return 2.*(1. + zeta)*xi*eta;
                case 14: // phi8 with zeta <- (-zeta)
                  return 2.*(1. + zeta)*eta*(1. - xi - eta);
                default:  std::cerr<<"Invalid i = " << i;
                }
            }
            // quadratic prism shape functions
          case PRISM6:   assert(T == L2_LAGRANGE);std::cerr<< "High order on first order elements only supported for L2 families";
            (void(0));
          case PRISM18: case PRISM20:  case PRISM21: {  assert (i< 18);
              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;// Point p2d(p[0],p[1]);
              double p1d = p[2];
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};
              return (MGFE<2,LAGRANGE_XFEM>::shape(TRI6,  SECOND, i1[i], p2d)*
                      fe_lagrange_1D_quadratic_shape(i0[i], p1d));
            }
            // G. Bedrosian, "Shape functions and integration formulas for
            // three-dimensional finite element analysis", Int. J. Numerical
            // Methods Engineering, vol 35, p. 95-108, 1992.
          case PYRAMID13:  {  assert (i< 13);
              const double xi   = p[0];  const double eta  = p[1];   const double zeta = p[2];
              const double eps  = 1.e-35;
              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              double den = (1. - zeta + eps);
              switch(i) {
                case 0:  return 0.25*(-xi - eta - 1.)*((1. - xi)*(1. - eta) - zeta + xi*eta*zeta/den);
                case 1:  return 0.25*(-eta + xi - 1.)*((1. + xi)*(1. - eta) - zeta - xi*eta*zeta/den);
                case 2:  return 0.25*(xi + eta - 1.)*((1. + xi)*(1. + eta) - zeta + xi*eta*zeta/den);
                case 3:  return 0.25*(eta - xi - 1.)*((1. - xi)*(1. + eta) - zeta - xi*eta*zeta/den);
                case 4:  return zeta*(2.*zeta - 1.);
                case 5:  return 0.5*(1. + xi - zeta)*(1. - xi - zeta)*(1. - eta - zeta)/den;
                case 6:  return 0.5*(1. + eta - zeta)*(1. - eta - zeta)*(1. + xi - zeta)/den;
                case 7:  return 0.5*(1. + xi - zeta)*(1. - xi - zeta)*(1. + eta - zeta)/den;
                case 8:  return 0.5*(1. + eta - zeta)*(1. - eta - zeta)*(1. - xi - zeta)/den;
                case 9:  return zeta*(1. - xi - zeta)*(1. - eta - zeta)/den;
                case 10:  return zeta*(1. + xi - zeta)*(1. - eta - zeta)/den;
                case 11:  return zeta*(1. + eta - zeta)*(1. + xi - zeta)/den;
                case 12:  return zeta*(1. - xi - zeta)*(1. + eta - zeta)/den;
                default:  std::cerr<<"Invalid i = " << i;
                }
            }
            // Quadratic shape functions, as defined in R. Graglia, "Higher order
            // bases on pyramidal elements", IEEE Trans Antennas and Propagation,
            // vol 47, no 5, May 1999.
          case PYRAMID5:  assert(T == L2_LAGRANGE);std::cerr<< "High order on first order elements only supported for L2 families";
            (void(0));
          case PYRAMID14: case PYRAMID18:   {     assert (i< 14);
              const double xi   = p[0];  const double eta  = p[1];  const double zeta = p[2];
              const double eps  = 1.e-35;
              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              double
                p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right
              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              double  den = (-1. + zeta + eps),  den2 = den*den;
              switch(i)  {
                case 0:  return p4*p1*(xi*eta - zeta + zeta*zeta)/den2;
                case 1:  return -p1*p2*(xi*eta + zeta - zeta*zeta)/den2;
                case 2:  return p2*p3*(xi*eta - zeta + zeta*zeta)/den2;
                case 3:  return -p3*p4*(xi*eta + zeta - zeta*zeta)/den2;
                case 4:  return zeta*(2.*zeta - 1.);
                case 5:  return -4.*p2*p1*p4*eta/den2;
                case 6:  return 4.*p1*p2*p3*xi/den2;
                case 7:  return 4.*p2*p3*p4*eta/den2;
                case 8:  return -4.*p3*p4*p1*xi/den2;
                case 9:  return -4.*p1*p4*zeta/den;
                case 10:  return -4.*p2*p1*zeta/den;
                case 11:  return -4.*p3*p2*zeta/den;
                case 12:  return -4.*p4*p3*zeta/den;
                case 13:  return 16.*p1*p2*p3*p4/den2;
                default:   std::cerr<<"Invalid i = " << i;
                }
            }
          default:   std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
          }
      }
    case THIRD: {  //ttttttttttttttttttttttttttttttttttttttttttttttttttttt-->
        switch (type) {
            // quadratic Lagrange shape functions with cubic bubbles
          case PRISM20:  { assert (i< 20);
              // Compute Prism21 shape functions as a tensor-product
              // of a triangle and an edge, then redistribute the
              // central bubble function over the other Tri6 nodes
              // around it (in a way consistent with the Tri7 shape
              // function definitions).
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1};
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;// Point p2d(p[0],p[1]);
              double p1d = p[2];
              const double mainval = MGFE<3,LAGRANGE_XFEM>::shape(PRISM21, THIRD, i, p);
              if (i0[i] != 2)  return mainval;
              const double bubbleval =  MGFE<2,LAGRANGE_XFEM>::shape(TRI7, THIRD, 6, p2d) *  fe_lagrange_1D_quadratic_shape(2, p1d);
              if (i < 12)  return mainval - bubbleval / 9; // vertices
              return mainval + bubbleval * (double(4) / 9);
            }
            // quadratic Lagrange shape functions with cubic bubbles
          case PRISM21: {  assert (i< 21);
              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;// Point p2d(p[0],p[1]);
              double p1d = p[2];
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5, 6, 6, 6};
              return (MGFE<2,LAGRANGE_XFEM>::shape(TRI7, THIRD, i1[i], p2d)* fe_lagrange_1D_quadratic_shape(i0[i], p1d)); }
            // Weird rational shape functions with weirder bubbles...
          case PYRAMID18: {   assert (i< 18);
              const double xi   = p[0]; const double eta  = p[1];  const double zeta = p[2];  const double eps  = 1.e-35;
              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              const double  p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right
              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              const double  den = (-1. + zeta + eps),  den2 = den*den;
              // Bubble functions on triangular sides.  We actually
              // have a degree of freedom to play with here, and I'm
              // not certain how best to use it, so let's leave it
              // as a variable in case we figure that out later.
              constexpr double alpha = 0.5;
              const double  bub_f1 = ((1-alpha)*(1-zeta) + alpha*(-eta)),
                bub_f2 = ((1-alpha)*(1-zeta) + alpha*(xi)),
                bub_f3 = ((1-alpha)*(1-zeta) + alpha*(eta)),
                bub_f4 = ((1-alpha)*(1-zeta) + alpha*(-xi));

              const double  bub1 = bub_f1*p1*p2*p4*zeta/den2,
                bub2 = bub_f2*p1*p2*p3*zeta/den2,
                bub3 = bub_f3*p2*p3*p4*zeta/den2,
                bub4 = bub_f4*p1*p3*p4*zeta/den2;

              switch(i) {
                case 0: return p4*p1*(xi*eta - zeta + zeta*zeta)/den2 + 3*(bub1+bub4);
                case 1: return -p1*p2*(xi*eta + zeta - zeta*zeta)/den2 + 3*(bub1+bub2);
                case 2: return p2*p3*(xi*eta - zeta + zeta*zeta)/den2 + 3*(bub2+bub3);
                case 3: return -p3*p4*(xi*eta + zeta - zeta*zeta)/den2 + 3*(bub3+bub4);
                case 4: return zeta*(2.*zeta - 1.) + 3*(bub1+bub2+bub3+bub4);
                case 5: return -4.*p2*p1*p4*eta/den2 - 12*bub1;
                case 6: return 4.*p1*p2*p3*xi/den2 - 12*bub2;
                case 7: return 4.*p2*p3*p4*eta/den2 - 12*bub3;
                case 8: return -4.*p3*p4*p1*xi/den2 - 12*bub4;
                case 9: return -4.*p1*p4*zeta/den - 12*(bub1+bub4);
                case 10: return -4.*p2*p1*zeta/den - 12*(bub1+bub2);
                case 11: return -4.*p3*p2*zeta/den - 12*(bub2+bub3);
                case 12: return -4.*p4*p3*zeta/den - 12*(bub3+bub4);
                case 13: return 16.*p1*p2*p3*p4/den2;
                case 14: return 27*bub1;
                case 15: return 27*bub2;
                case 16: return 27*bub3;
                case 17: return 27*bub4;
                default:  std::cerr<<"Invalid i = " << i;
                }
            }

            // quadratic Lagrange shape functions with cubic bubbles
          case TET14: {  assert (i< 14);
              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const double zeta1 = p[0]; const double zeta2 = p[1];  const double zeta3 = p[2];
              const double zeta0 = 1. - zeta1 - zeta2 - zeta3;
              // Bubble functions (not yet scaled) on side nodes
              const double bubble_012 = zeta0*zeta1*zeta2;   const double bubble_013 = zeta0*zeta1*zeta3;
              const double bubble_123 = zeta1*zeta2*zeta3;   const double bubble_023 = zeta0*zeta2*zeta3;
              switch(i)                {
                case 0: return zeta0*(2.*zeta0 - 1.) + 3.*(bubble_012+bubble_013+bubble_023);
                case 1: return zeta1*(2.*zeta1 - 1.) + 3.*(bubble_012+bubble_013+bubble_123);
                case 2: return zeta2*(2.*zeta2 - 1.) + 3.*(bubble_012+bubble_023+bubble_123);
                case 3: return zeta3*(2.*zeta3 - 1.) + 3.*(bubble_013+bubble_023+bubble_123);
                case 4:   return 4.*zeta0*zeta1 - 12.*(bubble_012+bubble_013);
                case 5:   return 4.*zeta1*zeta2 - 12.*(bubble_012+bubble_123);
                case 6:   return 4.*zeta2*zeta0 - 12.*(bubble_012+bubble_023);
                case 7:   return 4.*zeta0*zeta3 - 12.*(bubble_013+bubble_023);
                case 8:   return 4.*zeta1*zeta3 - 12.*(bubble_013+bubble_123);
                case 9: return 4.*zeta2*zeta3 - 12.*(bubble_023+bubble_123);
                case 10: return 27.*bubble_012;
                case 11:  return 27.*bubble_013;
                case 12:  return 27.*bubble_123;
                case 13:  return 27.*bubble_023;
                default:  std::cerr<<"Invalid i = " << i;
                }
            }
          default: std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
          }
      }//  <---------ttttttttttttttttttttttttttttttttttttttttttttttttttttt thirdorder
      // unsupported order
    default:  std::cerr<<"ERROR: Unsupported 3D FE order!: " << order; return 0.;
    }
}
 

// ======================================================================================================
template <FEFamily T>
double fe_lagrange_3D_shape_deriv(
  const ElemType type,
  const Order order,
  const unsigned int i,
  const unsigned int j,
  const double  p[]
){// =============================================================================
// #if LIBMESH_DIM == 3
  assert(j< 3);
  switch (order)  {
     case CONSTANT: return 0.;
    case FIRST: {// linear Lagrange shape functions -------------------------
        switch (type) {
            // trilinear hexahedral shape functions
          case HEX8:  case HEX20:  case HEX27:   {   assert(i< 8);
              // Compute hex shape functions as a tensor-product
              const double xi   = p[0];  const double eta  = p[1];  const double zeta = p[2];
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};
              switch(j) {
                case 0: return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)*
                                fe_lagrange_1D_linear_shape      (i1[i], eta)*
                                fe_lagrange_1D_linear_shape      (i2[i], zeta));
                case 1:
                  return (fe_lagrange_1D_linear_shape      (i0[i], xi)*
                          fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_linear_shape      (i2[i], zeta));
                case 2:
                  return (fe_lagrange_1D_linear_shape      (i0[i], xi)*
                          fe_lagrange_1D_linear_shape      (i1[i], eta)*
                          fe_lagrange_1D_linear_shape_deriv(i2[i], 0, zeta));
                default:   std::cerr<<"Invalid j = " << j;
                }
            }
            // linear tetrahedral shape functions --------------------------
          case TET4:   case TET10:  case TET14:  {
              assert(i< 4);
              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const double dzeta0dxi = -1.; const double dzeta1dxi =  1.; const double dzeta2dxi =  0.;  const double dzeta3dxi =  0.;
              const double dzeta0deta = -1.; const double dzeta1deta =  0.; const double dzeta2deta =  1.; const double dzeta3deta =  0.;
              const double dzeta0dzeta = -1.;const double dzeta1dzeta =  0.; const double dzeta2dzeta =  0.; const double dzeta3dzeta =  1.;
              switch (j) {
                case 0: {// d()/dxi
                    switch(i) {
                      case 0: return dzeta0dxi;
                      case 1: return dzeta1dxi;
                      case 2: return dzeta2dxi;
                      case 3: return dzeta3dxi;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 1: {// d()/deta
                    switch(i)  {
                      case 0: return dzeta0deta;
                      case 1: return dzeta1deta;
                      case 2: return dzeta2deta;
                      case 3: return dzeta3deta;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 2:  {
                    switch(i){  // d()/dzeta
                      case 0:   return dzeta0dzeta;
                      case 1:   return dzeta1dzeta;
                      case 2:   return dzeta2dzeta;
                      case 3:  return dzeta3dzeta;
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }
                default:  std::cerr<<"Invalid shape function derivative j = " << j;
                }
            }
            // linear prism shape functions  ------------------------
          case PRISM6:  case PRISM15: case PRISM18:   case PRISM20:  case PRISM21:   {   assert(i< 6);
              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;
              // Point p2d(p[0],p[1]);
              double p1d = p[2];
              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};
              switch (j)  {
                case 0: // d()/dxi
                  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI3,FIRST,i1[i],0,p2d)*fe_lagrange_1D_linear_shape(i0[i], p1d));
                case 1: // d()/deta
                  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI3,FIRST,i1[i],1,p2d)*fe_lagrange_1D_linear_shape(i0[i], p1d));
                case 2:// d()/dzeta
                  return (MGFE<2,LAGRANGE_XFEM>::shape(TRI3,  FIRST, i1[i], p2d)*fe_lagrange_1D_linear_shape_deriv(i0[i], 0, p1d));

                default:
                  std::cerr<<"Invalid shape function derivative j = " << j;
                }
            }
            // linear pyramid shape functions
          case PYRAMID5:  case PYRAMID13: case PYRAMID14:  case PYRAMID18:  {  assert(i< 5);
              const double xi   = p[0];  const double eta  = p[1]; const double zeta = p[2];   const double eps  = 1.e-35;
              switch (j) {
                case 0: // d/dxi
                  switch(i){
                    case 0:  return  .25*(zeta + eta - 1.)/((1. - zeta) + eps);
                    case 1:  return -.25*(zeta + eta - 1.)/((1. - zeta) + eps);
                    case 2:  return -.25*(zeta - eta - 1.)/((1. - zeta) + eps);
                    case 3:  return  .25*(zeta - eta - 1.)/((1. - zeta) + eps);
                    case 4:  return 0;
                    default:  std::cerr<<"Invalid i = " << i;
                    }
                case 1:// d/deta
                  switch(i)  {
                    case 0: return  .25*(zeta + xi - 1.)/((1. - zeta) + eps);
                    case 1: return  .25*(zeta - xi - 1.)/((1. - zeta) + eps);
                    case 2: return -.25*(zeta - xi - 1.)/((1. - zeta) + eps);
                    case 3: return -.25*(zeta + xi - 1.)/((1. - zeta) + eps);
                    case 4: return 0;
                    default: std::cerr<<"Invalid i = " << i;
                    }
                  
                case 2: {// d/dzeta
                    // We computed the derivatives with general eps and
                    // then let eps tend to zero in the numerators...
                    double num = zeta*(2. - zeta) - 1.,  den = (1. - zeta + eps)*(1. - zeta + eps);
                    switch(i) {
                      case 0:  case 2:  return .25*(num + xi*eta)/den;
                      case 1:  case 3:  return .25*(num - xi*eta)/den;
                      case 4:  return 1.;
                      default: std::cerr<<"Invalid i = " << i;
                      }
                  }
                default: std::cerr<<"Invalid j = " << j;
                }
            }
          default:  std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
          }
      }
    case SECOND:  { // quadratic Lagrange shape functions -------------------
        switch (type)  {
            // serendipity hexahedral quadratic shape functions
          case HEX20:  { assert(i< 20);
              const double xi   = p[0]; const double eta  = p[1];  const double zeta = p[2];
              // these functions are defined for (x,y,z) in [0,1]^3
              // so transform the locations
              const double x = .5*(xi   + 1.); const double y = .5*(eta  + 1.); const double z = .5*(zeta + 1.);
              // and don't forget the chain rule!
              switch (j) {
                case 0: // d/dx*dx/dxi
                  switch (i) {
                    case 0: return .5*(1. - y)*(1. - z)*((1. - x)*(-2.) + (-1.)*(1. - 2.*x - 2.*y - 2.*z));
                    case 1: return .5*(1. - y)*(1. - z)*(x*(2.) + (1.)*(2.*x - 2.*y - 2.*z - 1.));
                    case 2: return .5*y*(1. - z)*(x*(2.) + (1.)*(2.*x + 2.*y - 2.*z - 3.));
                    case 3: return .5*y*(1. - z)*((1. - x)*(-2.) + (-1.)*(2.*y - 2.*x - 2.*z - 1.));
                    case 4: return .5*(1. - y)*z*((1. - x)*(-2.) + (-1.)*(2.*z - 2.*x - 2.*y - 1.));
                    case 5: return .5*(1. - y)*z*(x*(2.) +  (1.)*(2.*x - 2.*y + 2.*z - 3.));
                    case 6: return .5*y*z*(x*(2.) + (1.)*(2.*x + 2.*y + 2.*z - 5.));
                    case 7: return .5*y*z*((1. - x)*(-2.) +  (-1.)*(2.*y - 2.*x + 2.*z - 3.));
                    case 8: return 2.*(1. - y)*(1. - z)*(1. - 2.*x);
                    case 9: return 2.*y*(1. - y)*(1. - z);
                    case 10: return 2.*y*(1. - z)*(1. - 2.*x);
                    case 11: return 2.*y*(1. - y)*(1. - z)*(-1.);
                    case 12: return 2.*(1. - y)*z*(1. - z)*(-1.);
                    case 13: return 2.*(1. - y)*z*(1. - z);
                    case 14: return 2.*y*z*(1. - z);
                    case 15: return 2.*y*z*(1. - z)*(-1.);
                    case 16: return 2.*(1. - y)*z*(1. - 2.*x);
                    case 17: return 2.*y*(1. - y)*z;
                    case 18: return 2.*y*z*(1. - 2.*x);
                    case 19: return 2.*y*(1. - y)*z*(-1.);
                    default:   std::cerr<<"Invalid i = " << i;
                    }
                  // d/dy*dy/deta
                case 1:
                  switch (i)
                    {
                    case 0:
                      return .5*(1. - x)*(1. - z)*((1. - y)*(-2.) +
                                                   (-1.)*(1. - 2.*x - 2.*y - 2.*z));

                    case 1:
                      return .5*x*(1. - z)*((1. - y)*(-2.) +
                                            (-1.)*(2.*x - 2.*y - 2.*z - 1.));

                    case 2:
                      return .5*x*(1. - z)*(y*(2.) +
                                            (1.)*(2.*x + 2.*y - 2.*z - 3.));

                    case 3:
                      return .5*(1. - x)*(1. - z)*(y*(2.) +
                                                   (1.)*(2.*y - 2.*x - 2.*z - 1.));

                    case 4:
                      return .5*(1. - x)*z*((1. - y)*(-2.) +
                                            (-1.)*(2.*z - 2.*x - 2.*y - 1.));

                    case 5:
                      return .5*x*z*((1. - y)*(-2.) +
                                     (-1.)*(2.*x - 2.*y + 2.*z - 3.));

                    case 6:
                      return .5*x*z*(y*(2.) +
                                     (1.)*(2.*x + 2.*y + 2.*z - 5.));

                    case 7:
                      return .5*(1. - x)*z*(y*(2.) +
                                            (1.)*(2.*y - 2.*x + 2.*z - 3.));

                    case 8:
                      return 2.*x*(1. - x)*(1. - z)*(-1.);

                    case 9:
                      return 2.*x*(1. - z)*(1. - 2.*y);

                    case 10:
                      return 2.*x*(1. - x)*(1. - z);

                    case 11:
                      return 2.*(1. - x)*(1. - z)*(1. - 2.*y);

                    case 12:
                      return 2.*(1. - x)*z*(1. - z)*(-1.);

                    case 13:
                      return 2.*x*z*(1. - z)*(-1.);

                    case 14:
                      return 2.*x*z*(1. - z);

                    case 15:
                      return 2.*(1. - x)*z*(1. - z);

                    case 16:
                      return 2.*x*(1. - x)*z*(-1.);

                    case 17:
                      return 2.*x*z*(1. - 2.*y);

                    case 18:
                      return 2.*x*(1. - x)*z;

                    case 19:
                      return 2.*(1. - x)*z*(1. - 2.*y);

                    default:
                      std::cerr<<"Invalid i = " << i;
                    }


                  // d/dz*dz/dzeta
                case 2:
                  switch (i)
                    {
                    case 0:
                      return .5*(1. - x)*(1. - y)*((1. - z)*(-2.) +
                                                   (-1.)*(1. - 2.*x - 2.*y - 2.*z));

                    case 1:
                      return .5*x*(1. - y)*((1. - z)*(-2.) +
                                            (-1.)*(2.*x - 2.*y - 2.*z - 1.));

                    case 2:
                      return .5*x*y*((1. - z)*(-2.) +
                                     (-1.)*(2.*x + 2.*y - 2.*z - 3.));

                    case 3:
                      return .5*(1. - x)*y*((1. - z)*(-2.) +
                                            (-1.)*(2.*y - 2.*x - 2.*z - 1.));

                    case 4:
                      return .5*(1. - x)*(1. - y)*(z*(2.) +
                                                   (1.)*(2.*z - 2.*x - 2.*y - 1.));

                    case 5:
                      return .5*x*(1. - y)*(z*(2.) +
                                            (1.)*(2.*x - 2.*y + 2.*z - 3.));

                    case 6:
                      return .5*x*y*(z*(2.) +
                                     (1.)*(2.*x + 2.*y + 2.*z - 5.));

                    case 7:
                      return .5*(1. - x)*y*(z*(2.) +
                                            (1.)*(2.*y - 2.*x + 2.*z - 3.));

                    case 8:
                      return 2.*x*(1. - x)*(1. - y)*(-1.);

                    case 9:
                      return 2.*x*y*(1. - y)*(-1.);

                    case 10:
                      return 2.*x*(1. - x)*y*(-1.);

                    case 11:
                      return 2.*(1. - x)*y*(1. - y)*(-1.);

                    case 12:
                      return 2.*(1. - x)*(1. - y)*(1. - 2.*z);

                    case 13:
                      return 2.*x*(1. - y)*(1. - 2.*z);

                    case 14:
                      return 2.*x*y*(1. - 2.*z);

                    case 15:
                      return 2.*(1. - x)*y*(1. - 2.*z);

                    case 16:
                      return 2.*x*(1. - x)*(1. - y);

                    case 17:
                      return 2.*x*y*(1. - y);

                    case 18:
                      return 2.*x*(1. - x)*y;

                    case 19:
                      return 2.*(1. - x)*y*(1. - y);

                    default:
                      std::cerr<<"Invalid i = " << i;
                    }

                default:
                  std::cerr<<"Invalid shape function derivative j = " << j;
                }
            }

            // triquadratic hexahedral shape functions
          case HEX8:  assert(T == L2_LAGRANGE);std::cerr<< "High order on first order elements only supported for L2 families";
            (void(0));
          case HEX27: {  assert(i< 27);
              // Compute hex shape functions as a tensor-product
              const double xi   = p[0];  const double eta  = p[1];  const double zeta = p[2];
              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};
              switch(j) {
                case 0: return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*
                                fe_lagrange_1D_quadratic_shape      (i1[i], eta)*
                                fe_lagrange_1D_quadratic_shape      (i2[i], zeta));
                case 1: return (fe_lagrange_1D_quadratic_shape      (i0[i], xi)*
                                fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta)*
                                fe_lagrange_1D_quadratic_shape      (i2[i], zeta));
                case 2:
                  return (fe_lagrange_1D_quadratic_shape      (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape      (i1[i], eta)*
                          fe_lagrange_1D_quadratic_shape_deriv(i2[i], 0, zeta));
                default:
                  std::cerr<<"Invalid j = " << j;
                }
            }

            // quadratic tetrahedral shape functions
          case TET4:  assert(T == L2_LAGRANGE);std::cerr<< "High order on first order elements only supported for L2 families";
            (void(0));
          case TET10:   case TET14:  {   assert(i< 10);
              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const double zeta1 = p[0];  const double zeta2 = p[1]; const double zeta3 = p[2];
              const double zeta0 = 1. - zeta1 - zeta2 - zeta3;
              const double dzeta0dxi = -1.; const double dzeta1dxi =  1.;
              const double dzeta2dxi =  0.; const double dzeta3dxi =  0.;
              const double dzeta0deta = -1.; const double dzeta1deta =  0.;
              const double dzeta2deta =  1.; const double dzeta3deta =  0.;
              const double dzeta0dzeta = -1.; const double dzeta1dzeta =  0.;
              const double dzeta2dzeta =  0.; const double dzeta3dzeta =  1.;
              switch (j)   {
                  // d()/dxi
                case 0: {
                    switch(i)   {
                      case 0:  return (4.*zeta0 - 1.)*dzeta0dxi;
                      case 1:  return (4.*zeta1 - 1.)*dzeta1dxi;
                      case 2:  return (4.*zeta2 - 1.)*dzeta2dxi;
                      case 3:  return (4.*zeta3 - 1.)*dzeta3dxi;
                      case 4:  return 4.*(zeta0*dzeta1dxi + dzeta0dxi*zeta1);
                      case 5:  return 4.*(zeta1*dzeta2dxi + dzeta1dxi*zeta2);
                      case 6:  return 4.*(zeta0*dzeta2dxi + dzeta0dxi*zeta2);
                      case 7:  return 4.*(zeta0*dzeta3dxi + dzeta0dxi*zeta3);
                      case 8:  return 4.*(zeta1*dzeta3dxi + dzeta1dxi*zeta3);
                      case 9:  return 4.*(zeta2*dzeta3dxi + dzeta2dxi*zeta3);
                      default:    std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 1:  {   // d()/deta
                    switch(i)  {
                      case 0: return (4.*zeta0 - 1.)*dzeta0deta;
                      case 1: return (4.*zeta1 - 1.)*dzeta1deta;
                      case 2: return (4.*zeta2 - 1.)*dzeta2deta;
                      case 3: return (4.*zeta3 - 1.)*dzeta3deta;
                      case 4: return 4.*(zeta0*dzeta1deta + dzeta0deta*zeta1);
                      case 5: return 4.*(zeta1*dzeta2deta + dzeta1deta*zeta2);
                      case 6: return 4.*(zeta0*dzeta2deta + dzeta0deta*zeta2);
                      case 7: return 4.*(zeta0*dzeta3deta + dzeta0deta*zeta3);
                      case 8: return 4.*(zeta1*dzeta3deta + dzeta1deta*zeta3);
                      case 9: return 4.*(zeta2*dzeta3deta + dzeta2deta*zeta3);
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }
                 
                case 2: { // d()/dzeta
                    switch(i) {
                      case 0: return (4.*zeta0 - 1.)*dzeta0dzeta;
                      case 1: return (4.*zeta1 - 1.)*dzeta1dzeta;
                      case 2: return (4.*zeta2 - 1.)*dzeta2dzeta;
                      case 3: return (4.*zeta3 - 1.)*dzeta3dzeta;
                      case 4: return 4.*(zeta0*dzeta1dzeta + dzeta0dzeta*zeta1);
                      case 5: return 4.*(zeta1*dzeta2dzeta + dzeta1dzeta*zeta2);
                      case 6: return 4.*(zeta0*dzeta2dzeta + dzeta0dzeta*zeta2);
                      case 7: return 4.*(zeta0*dzeta3dzeta + dzeta0dzeta*zeta3);
                      case 8: return 4.*(zeta1*dzeta3dzeta + dzeta1dzeta*zeta3);
                      case 9: return 4.*(zeta2*dzeta3dzeta + dzeta2dzeta*zeta3);
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                default:  std::cerr<<"Invalid j = " << j;
                }
            }
            // "serendipity" prism
          case PRISM15: {  assert(i< 15);
              const double xi   = p[0];  const double eta  = p[1]; const double zeta = p[2];

              switch (j) {
                  // d()/dxi
                case 0: {
                    switch(i)   {
                      case 0:  return (2.*xi + 2.*eta + 0.5*zeta - 1.)*(1. - zeta);
                      case 1:  return (2.*xi - 1. - 0.5*zeta)*(1. - zeta);
                      case 2:  return 0.;
                      case 3:  return (2.*xi + 2.*eta - 0.5*zeta - 1.)*(1. + zeta);
                      case 4:  return (2.*xi - 1. + 0.5*zeta)*(1. + zeta);
                      case 5:  return 0.;
                      case 6:  return (4.*xi + 2.*eta - 2.)*(zeta - 1.);
                      case 7:  return -2.*(zeta - 1.)*eta;
                      case 8:  return 2.*(zeta - 1.)*eta;
                      case 9:  return (zeta - 1.)*(1. + zeta);
                      case 10: return (1. - zeta)*(1. + zeta);
                      case 11: return 0.;
                      case 12: return (-4.*xi - 2.*eta + 2.)*(1. + zeta);
                      case 13: return 2.*(1. + zeta)*eta;
                      case 14: return -2.*(1. + zeta)*eta;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 1:   {   // d()/deta
                    switch(i) {
                      case 0:  return (2.*xi + 2.*eta + 0.5*zeta - 1.)*(1. - zeta);
                      case 1:  return 0.;
                      case 2:  return (2.*eta - 1. - 0.5*zeta)*(1. - zeta);
                      case 3:  return (2.*xi + 2.*eta - 0.5*zeta - 1.)*(1. + zeta);
                      case 4:  return 0.;
                      case 5:  return (2.*eta - 1. + 0.5*zeta)*(1. + zeta);
                      case 6:  return 2.*(zeta - 1.)*xi;
                      case 7:  return 2.*(1. - zeta)*xi;
                      case 8:  return (2.*xi + 4.*eta - 2.)*(zeta - 1.);
                      case 9:  return (zeta - 1.)*(1. + zeta);
                      case 10: return 0.;
                      case 11: return (1. - zeta)*(1. + zeta);
                      case 12: return -2.*(1. + zeta)*xi;
                      case 13: return 2.*(1. + zeta)*xi;
                      case 14: return (-2.*xi - 4.*eta + 2.)*(1. + zeta);
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 2: {     // d()/dzeta
                    switch(i) {
                      case 0:  return (-xi - eta - zeta + 0.5)*(xi + eta - 1.);
                      case 1:  return -0.5*xi*(2.*xi - 1. - 2.*zeta);
                      case 2:  return -0.5*eta*(2.*eta - 1. - 2.*zeta);
                      case 3:  return (xi + eta - zeta - 0.5)*(xi + eta - 1.);
                      case 4:  return 0.5*xi*(2.*xi - 1. + 2.*zeta);
                      case 5:  return 0.5*eta*(2.*eta - 1. + 2.*zeta);
                      case 6:  return 2.*xi*(xi + eta - 1.);
                      case 7:  return -2.*xi*eta;
                      case 8:  return 2.*eta*(xi + eta - 1.);
                      case 9:  return 2.*zeta*(xi + eta - 1.);
                      case 10: return -2.*xi*zeta;
                      case 11: return -2.*eta*zeta;
                      case 12: return 2.*xi*(1. - xi - eta);
                      case 13: return 2.*xi*eta;
                      case 14: return 2.*eta*(1. - xi - eta);
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                default:  std::cerr<<"Invalid j = " << j;
                }
            }



            // quadratic prism shape functions
          case PRISM6:   assert(T == L2_LAGRANGE);std::cerr<<  "High order on first order elements only supported for L2 families";
            (void(0));
          case PRISM18:
          case PRISM20:
          case PRISM21:
            {
              assert(i< 18);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;
              double p1d = p[2];

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};

              switch (j)
                {
                  // d()/dxi
                case 0:
                  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI6,  SECOND, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/deta
                case 1:
                  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI6,  SECOND, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/dzeta
                case 2:
                  return (MGFE<2,LAGRANGE_XFEM>::shape(TRI6,  SECOND, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                default:
                  std::cerr<<"Invalid shape function derivative j = " << j;
                }
            }

            // G. Bedrosian, "Shape functions and integration formulas for
            // three-dimensional finite element analysis", Int. J. Numerical
            // Methods Engineering, vol 35, p. 95-108, 1992.
          case PYRAMID13:
            {
              assert(i< 13);

              const double xi   = p[0];
              const double eta  = p[1];
              const double zeta = p[2];
              const double eps  = 1.e-35;

              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              double
                den = (-1. + zeta + eps),
                den2 = den*den,
                xi2 = xi*xi,
                eta2 = eta*eta,
                zeta2 = zeta*zeta,
                zeta3 = zeta2*zeta;

              switch (j)
                {
                  // d/dxi
                case 0:
                  switch(i)
                    {
                    case 0:
                      return 0.25*(-zeta - eta + 2.*eta*zeta - 2.*xi + 2.*zeta*xi + 2.*eta*xi + zeta2 + eta2)/den;

                    case 1:
                      return -0.25*(-zeta - eta + 2.*eta*zeta + 2.*xi - 2.*zeta*xi - 2.*eta*xi + zeta2 + eta2)/den;

                    case 2:
                      return -0.25*(-zeta + eta - 2.*eta*zeta + 2.*xi - 2.*zeta*xi + 2.*eta*xi + zeta2 + eta2)/den;

                    case 3:
                      return 0.25*(-zeta + eta - 2.*eta*zeta - 2.*xi + 2.*zeta*xi - 2.*eta*xi + zeta2 + eta2)/den;

                    case 4:
                      return 0.;

                    case 5:
                      return -(-1. + eta + zeta)*xi/den;

                    case 6:
                      return 0.5*(-1. + eta + zeta)*(1. + eta - zeta)/den;

                    case 7:
                      return (1. + eta - zeta)*xi/den;

                    case 8:
                      return -0.5*(-1. + eta + zeta)*(1. + eta - zeta)/den;

                    case 9:
                      return -(-1. + eta + zeta)*zeta/den;

                    case 10:
                      return (-1. + eta + zeta)*zeta/den;

                    case 11:
                      return -(1. + eta - zeta)*zeta/den;

                    case 12:
                      return (1. + eta - zeta)*zeta/den;

                    default:
                      std::cerr<<"Invalid i = " << i;
                    }

                  // d/deta
                case 1:
                  switch(i)
                    {
                    case 0:
                      return 0.25*(-zeta - 2.*eta + 2.*eta*zeta - xi + 2.*zeta*xi + 2.*eta*xi + zeta2 + xi2)/den;

                    case 1:
                      return -0.25*(zeta + 2.*eta - 2.*eta*zeta - xi + 2.*zeta*xi + 2.*eta*xi - zeta2 - xi2)/den;

                    case 2:
                      return -0.25*(-zeta + 2.*eta - 2.*eta*zeta + xi - 2.*zeta*xi + 2.*eta*xi + zeta2 + xi2)/den;

                    case 3:
                      return 0.25*(zeta - 2.*eta + 2.*eta*zeta + xi - 2.*zeta*xi + 2.*eta*xi - zeta2 - xi2)/den;

                    case 4:
                      return 0.;

                    case 5:
                      return -0.5*(-1. + xi + zeta)*(1. + xi - zeta)/den;

                    case 6:
                      return (1. + xi - zeta)*eta/den;

                    case 7:
                      return 0.5*(-1. + xi + zeta)*(1. + xi - zeta)/den;

                    case 8:
                      return -(-1. + xi + zeta)*eta/den;

                    case 9:
                      return -(-1. + xi + zeta)*zeta/den;

                    case 10:
                      return (1. + xi - zeta)*zeta/den;

                    case 11:
                      return -(1. + xi - zeta)*zeta/den;

                    case 12:
                      return (-1. + xi + zeta)*zeta/den;

                    default:
                      std::cerr<<"Invalid i = " << i;
                    }

                 
                case 2: { // d/dzeta
                    switch(i)  {
                      case 0:  return -0.25*(xi + eta + 1.)*(-1. + 2.*zeta - zeta2 + eta*xi)/den2;
                      case 1:  return 0.25*(eta - xi + 1.)*(1. - 2.*zeta + zeta2 + eta*xi)/den2;
                      case 2:  return 0.25*(xi + eta - 1.)*(-1. + 2.*zeta - zeta2 + eta*xi)/den2;
                      case 3:  return -0.25*(eta - xi - 1.)*(1. - 2.*zeta + zeta2 + eta*xi)/den2;
                      case 4:  return 4.*zeta - 1.;
                      case 5:  return 0.5*(-2 + eta + 6.*zeta + eta*xi2 + eta*zeta2 - 6.*zeta2 + 2.*zeta3 - 2.*eta*zeta)/den2;
                      case 6:  return -0.5*(2 - 6.*zeta + xi + xi*zeta2 + eta2*xi + 6.*zeta2 - 2.*zeta3 - 2.*zeta*xi)/den2;
                      case 7:  return -0.5*(2 + eta - 6.*zeta + eta*xi2 + eta*zeta2 + 6.*zeta2 - 2.*zeta3 - 2.*eta*zeta)/den2;
                      case 8:  return 0.5*(-2 + 6.*zeta + xi + xi*zeta2 + eta2*xi - 6.*zeta2 + 2.*zeta3 - 2.*zeta*xi)/den2;
                      case 9: return (1. - eta - 4.*zeta - xi - xi*zeta2 - eta*zeta2 + eta*xi + 5.*zeta2 - 2.*zeta3 + 2.*eta*zeta + 2.*zeta*xi)/den2;
                      case 10: return -(-1. + eta + 4.*zeta - xi - xi*zeta2 + eta*zeta2 + eta*xi - 5.*zeta2 + 2.*zeta3 - 2.*eta*zeta + 2.*zeta*xi)/den2;
                      case 11:  return (1. + eta - 4.*zeta + xi + xi*zeta2 + eta*zeta2 + eta*xi + 5.*zeta2 - 2.*zeta3 - 2.*eta*zeta - 2.*zeta*xi)/den2;
                      case 12: return -(-1. - eta + 4.*zeta + xi + xi*zeta2 - eta*zeta2 + eta*xi - 5.*zeta2 + 2.*zeta3 + 2.*eta*zeta - 2.*zeta*xi)/den2;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                default:  std::cerr<<"Invalid j = " << j;
                }
            }

            // Quadratic shape functions, as defined in R. Graglia, "Higher order
            // bases on pyramidal elements", IEEE Trans Antennas and Propagation,
            // vol 47, no 5, May 1999.
          case PYRAMID5:
            assert(T == L2_LAGRANGE);std::cerr<< "High order on first order elements only supported for L2 families";
            (void(0));
          case PYRAMID14:   case PYRAMID18:  {    assert(i< 14);
              const double xi   = p[0]; const double eta  = p[1];  const double zeta = p[2]; const double eps  = 1.e-35;
              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              double
                p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right
              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              double den = (-1. + zeta + eps),  den2 = den*den,  den3 = den2*den;

              switch (j){
                  // d/dxi
                case 0:
                  switch(i)
                    {
                    case 0:
                      return 0.5*p1*(-xi*eta + zeta - zeta*zeta + 2.*p4*eta)/den2;

                    case 1:
                      return -0.5*p1*(xi*eta + zeta - zeta*zeta + 2.*p2*eta)/den2;

                    case 2:
                      return 0.5*p3*(xi*eta - zeta + zeta*zeta + 2.*p2*eta)/den2;

                    case 3:
                      return -0.5*p3*(-xi*eta - zeta + zeta*zeta + 2.*p4*eta)/den2;

                    case 4:
                      return 0.;

                    case 5:
                      return 2.*p1*eta*xi/den2;

                    case 6:
                      return 2.*p1*p3*(xi + 2.*p2)/den2;

                    case 7:
                      return -2.*p3*eta*xi/den2;

                    case 8:
                      return -2.*p1*p3*(-xi + 2.*p4)/den2;

                    case 9:
                      return 2.*p1*zeta/den;

                    case 10:
                      return -2.*p1*zeta/den;

                    case 11:
                      return -2.*p3*zeta/den;

                    case 12:
                      return 2.*p3*zeta/den;

                    case 13:
                      return -8.*p1*p3*xi/den2;

                    default:
                      std::cerr<<"Invalid i = " << i;
                    }

                  // d/deta
                case 1:
                  switch(i)
                    {
                    case 0:
                      return -0.5*p4*(xi*eta - zeta + zeta*zeta - 2.*p1*xi)/den2;

                    case 1:
                      return 0.5*p2*(xi*eta + zeta - zeta*zeta - 2.*p1*xi)/den2;

                    case 2:
                      return 0.5*p2*(xi*eta - zeta + zeta*zeta + 2.*p3*xi)/den2;

                    case 3:
                      return -0.5*p4*(xi*eta + zeta - zeta*zeta + 2.*p3*xi)/den2;

                    case 4:
                      return 0.;

                    case 5:
                      return 2.*p2*p4*(eta - 2.*p1)/den2;

                    case 6:
                      return -2.*p2*xi*eta/den2;

                    case 7:
                      return 2.*p2*p4*(eta + 2.*p3)/den2;

                    case 8:
                      return 2.*p4*xi*eta/den2;

                    case 9:
                      return 2.*p4*zeta/den;

                    case 10:
                      return 2.*p2*zeta/den;

                    case 11:
                      return -2.*p2*zeta/den;

                    case 12:
                      return -2.*p4*zeta/den;

                    case 13:
                      return -8.*p2*p4*eta/den2;

                    default:
                      std::cerr<<"Invalid i = " << i;
                    }


                  // d/dzeta
                case 2:
                  {
                    switch(i)
                      {
                      case 0:
                        return -0.5*p1*(xi*eta - zeta + zeta*zeta)/den2
                          - 0.5*p4*(xi*eta - zeta + zeta*zeta)/den2
                          + p4*p1*(2.*zeta - 1)/den2
                          - 2.*p4*p1*(xi*eta - zeta + zeta*zeta)/den3;

                      case 1:
                        return 0.5*p2*(xi*eta + zeta - zeta*zeta)/den2
                          + 0.5*p1*(xi*eta + zeta - zeta*zeta)/den2
                          - p1*p2*(1 - 2.*zeta)/den2
                          + 2.*p1*p2*(xi*eta + zeta - zeta*zeta)/den3;

                      case 2:
                        return -0.5*p3*(xi*eta - zeta + zeta*zeta)/den2
                          - 0.5*p2*(xi*eta - zeta + zeta*zeta)/den2
                          + p2*p3*(2.*zeta - 1)/den2
                          - 2.*p2*p3*(xi*eta - zeta + zeta*zeta)/den3;

                      case 3:
                        return 0.5*p4*(xi*eta + zeta - zeta*zeta)/den2
                          + 0.5*p3*(xi*eta + zeta - zeta*zeta)/den2
                          - p3*p4*(1 - 2.*zeta)/den2
                          + 2.*p3*p4*(xi*eta + zeta - zeta*zeta)/den3;

                      case 4:
                        return 4.*zeta - 1.;

                      case 5:
                        return 2.*p4*p1*eta/den2
                          + 2.*p2*p4*eta/den2
                          + 2.*p1*p2*eta/den2
                          + 8.*p2*p1*p4*eta/den3;

                      case 6:
                        return -2.*p2*p3*xi/den2
                          - 2.*p1*p3*xi/den2
                          - 2.*p1*p2*xi/den2
                          - 8.*p1*p2*p3*xi/den3;

                      case 7:
                        return -2.*p3*p4*eta/den2
                          - 2.*p2*p4*eta/den2
                          - 2.*p2*p3*eta/den2
                          - 8.*p2*p3*p4*eta/den3;

                      case 8:
                        return 2.*p4*p1*xi/den2
                          + 2.*p1*p3*xi/den2
                          + 2.*p3*p4*xi/den2
                          + 8.*p3*p4*p1*xi/den3;

                      case 9:
                        return 2.*p4*zeta/den
                          + 2.*p1*zeta/den
                          - 4.*p1*p4/den
                          + 4.*p1*p4*zeta/den2;

                      case 10:
                        return 2.*p1*zeta/den
                          + 2.*p2*zeta/den
                          - 4.*p2*p1/den
                          + 4.*p2*p1*zeta/den2;

                      case 11:
                        return 2.*p2*zeta/den
                          + 2.*p3*zeta/den
                          - 4.*p3*p2/den
                          + 4.*p3*p2*zeta/den2;

                      case 12:
                        return 2.*p3*zeta/den
                          + 2.*p4*zeta/den
                          - 4.*p4*p3/den
                          + 4.*p4*p3*zeta/den2;

                      case 13:
                        return -8.*p2*p3*p4/den2
                          - 8.*p3*p4*p1/den2
                          - 8.*p2*p1*p4/den2
                          - 8.*p1*p2*p3/den2
                          - 32.*p1*p2*p3*p4/den3;

                      default:
                        std::cerr<<"Invalid i = " << i;
                      }
                  }

                default:
                  std::cerr<<"Invalid j = " << j;
                }
            }


          default:
            std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
          }
      }

    case THIRD: {
        switch (type) {
            // quadratic Lagrange shape functions with a cubic bubble
          case TET14:  {  assert(i< 14);
              // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
              const double zeta1 = p[0]; const double zeta2 = p[1]; const double zeta3 = p[2];
              const double zeta0 = 1. - zeta1 - zeta2 - zeta3;   const double dzeta0dxi = -1.;
              const double dzeta1dxi =  1.;  const double dzeta2dxi =  0.;   const double dzeta3dxi =  0.;
              const double dbubble012dxi = (zeta0-zeta1)*zeta2;   const double dbubble013dxi = (zeta0-zeta1)*zeta3;
              const double dbubble123dxi = zeta2*zeta3;           const double dbubble023dxi = -zeta2*zeta3;

              const double dzeta0deta = -1.;  const double dzeta1deta =  0.;  const double dzeta2deta =  1.;  const double dzeta3deta =  0.;
              const double dbubble012deta = (zeta0-zeta2)*zeta1;    const double dbubble013deta = -zeta1*zeta3;
              const double dbubble123deta = zeta1*zeta3;            const double dbubble023deta = (zeta0-zeta2)*zeta3;

              const double dzeta0dzeta = -1.; const double dzeta1dzeta =  0.; const double dzeta2dzeta =  0.; const double dzeta3dzeta =  1.;
              const double dbubble012dzeta = -zeta1*zeta2;              const double dbubble013dzeta = (zeta0-zeta3)*zeta1;
              const double dbubble123dzeta = zeta1*zeta2;               const double dbubble023dzeta = (zeta0-zeta3)*zeta2;

              switch (j)
                {
                  // d()/dxi
                case 0:
                  {
                    switch(i)
                      {
                      case 0:
                        return (4.*zeta0 - 1.)*dzeta0dxi + 3.*(dbubble012dxi+dbubble013dxi+dbubble023dxi);

                      case 1:
                        return (4.*zeta1 - 1.)*dzeta1dxi + 3.*(dbubble012dxi+dbubble013dxi+dbubble123dxi);

                      case 2:
                        return (4.*zeta2 - 1.)*dzeta2dxi + 3.*(dbubble012dxi+dbubble023dxi+dbubble123dxi);

                      case 3:
                        return (4.*zeta3 - 1.)*dzeta3dxi + 3.*(dbubble013dxi+dbubble023dxi+dbubble123dxi);

                      case 4:
                        return 4.*(zeta0*dzeta1dxi + dzeta0dxi*zeta1) - 12.*(dbubble012dxi+dbubble013dxi);

                      case 5:
                        return 4.*(zeta1*dzeta2dxi + dzeta1dxi*zeta2) - 12.*(dbubble012dxi+dbubble123dxi);

                      case 6:
                        return 4.*(zeta0*dzeta2dxi + dzeta0dxi*zeta2) - 12.*(dbubble012dxi+dbubble023dxi);

                      case 7:
                        return 4.*(zeta0*dzeta3dxi + dzeta0dxi*zeta3) - 12.*(dbubble013dxi+dbubble023dxi);

                      case 8:
                        return 4.*(zeta1*dzeta3dxi + dzeta1dxi*zeta3) - 12.*(dbubble013dxi+dbubble123dxi);

                      case 9:
                        return 4.*(zeta2*dzeta3dxi + dzeta2dxi*zeta3) - 12.*(dbubble023dxi+dbubble123dxi);

                      case 10:
                        return 27.*dbubble012dxi;

                      case 11:
                        return 27.*dbubble013dxi;

                      case 12:
                        return 27.*dbubble123dxi;

                      case 13:
                        return 27.*dbubble023dxi;

                      default:
                        std::cerr<<"Invalid i = " << i;
                      }
                  }

                  // d()/deta
                case 1:
                  {
                    switch(i)
                      {
                      case 0:
                        return (4.*zeta0 - 1.)*dzeta0deta + 3.*(dbubble012deta+dbubble013deta+dbubble023deta);;

                      case 1:
                        return (4.*zeta1 - 1.)*dzeta1deta + 3.*(dbubble012deta+dbubble013deta+dbubble123deta);

                      case 2:
                        return (4.*zeta2 - 1.)*dzeta2deta + 3.*(dbubble012deta+dbubble023deta+dbubble123deta);

                      case 3:
                        return (4.*zeta3 - 1.)*dzeta3deta + 3.*(dbubble013deta+dbubble023deta+dbubble123deta);

                      case 4:
                        return 4.*(zeta0*dzeta1deta + dzeta0deta*zeta1) - 12.*(dbubble012deta+dbubble013deta);

                      case 5:
                        return 4.*(zeta1*dzeta2deta + dzeta1deta*zeta2) - 12.*(dbubble012deta+dbubble123deta);

                      case 6:
                        return 4.*(zeta0*dzeta2deta + dzeta0deta*zeta2) - 12.*(dbubble012deta+dbubble023deta);

                      case 7:
                        return 4.*(zeta0*dzeta3deta + dzeta0deta*zeta3) - 12.*(dbubble013deta+dbubble023deta);

                      case 8:
                        return 4.*(zeta1*dzeta3deta + dzeta1deta*zeta3) - 12.*(dbubble013deta+dbubble123deta);

                      case 9:
                        return 4.*(zeta2*dzeta3deta + dzeta2deta*zeta3) - 12.*(dbubble023deta+dbubble123deta);

                      case 10:
                        return 27.*dbubble012deta;

                      case 11:
                        return 27.*dbubble013deta;

                      case 12:
                        return 27.*dbubble123deta;

                      case 13:
                        return 27.*dbubble023deta;

                      default:
                        std::cerr<<"Invalid i = " << i;
                      }
                  }

                  
                case 2:  {// d()/dzeta
                    switch(i) {
                      case 0: return (4.*zeta0 - 1.)*dzeta0dzeta + 3.*(dbubble012dzeta+dbubble013dzeta+dbubble023dzeta);
                      case 1: return (4.*zeta1 - 1.)*dzeta1dzeta + 3.*(dbubble012dzeta+dbubble013dzeta+dbubble123dzeta);
                      case 2: return (4.*zeta2 - 1.)*dzeta2dzeta + 3.*(dbubble012dzeta+dbubble023dzeta+dbubble123dzeta);
                      case 3: return (4.*zeta3 - 1.)*dzeta3dzeta + 3.*(dbubble013dzeta+dbubble023dzeta+dbubble123dzeta);
                      case 4: return 4.*(zeta0*dzeta1dzeta + dzeta0dzeta*zeta1) - 12.*(dbubble012dzeta+dbubble013dzeta);
                      case 5: return 4.*(zeta1*dzeta2dzeta + dzeta1dzeta*zeta2) - 12.*(dbubble012dzeta+dbubble123dzeta);
                      case 6: return 4.*(zeta0*dzeta2dzeta + dzeta0dzeta*zeta2) - 12.*(dbubble012dzeta+dbubble023dzeta);
                      case 7: return 4.*(zeta0*dzeta3dzeta + dzeta0dzeta*zeta3) - 12.*(dbubble013dzeta+dbubble023dzeta);
                      case 8: return 4.*(zeta1*dzeta3dzeta + dzeta1dzeta*zeta3) - 12.*(dbubble013dzeta+dbubble123dzeta);
                      case 9: return 4.*(zeta2*dzeta3dzeta + dzeta2dzeta*zeta3) - 12.*(dbubble023dzeta+dbubble123dzeta);
                      case 10: return 27.*dbubble012dzeta;
                      case 11: return 27.*dbubble013dzeta;
                      case 12: return 27.*dbubble123dzeta;
                      case 13: return 27.*dbubble023dzeta;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                default:  std::cerr<<"Invalid j = " << j;
                }
            }

          case PRISM20: {  assert(i< 20);
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1, 2};
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.; double p1d = p[2];
              const double mainval = MGFE<3,LAGRANGE_XFEM>::shape_deriv(PRISM21, THIRD, i, j, p);
              if (i0[i] != 2)  return mainval;
              double bubbleval = 0;
              switch (j) {
                  // d()/dxi
                case 0: bubbleval =  MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI7, THIRD, 6, 0, p2d)* fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;
                  // d()/deta
                case 1:  bubbleval =   MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI7, THIRD, 6, 1, p2d)*  fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;
                  // d()/dzeta
                case 2:  bubbleval = MGFE<2,LAGRANGE_XFEM>::shape(TRI7, THIRD, 6, p2d)* fe_lagrange_1D_quadratic_shape_deriv(2, 0, p1d);
                  break;
                default: std::cerr<<"Invalid shape function derivative j = " << j;
                }

              if (i < 12) // vertices
                return mainval - bubbleval / 9;

              return mainval + bubbleval * (double(4) / 9);
            }

          case PRISM21:{ assert(i< 21);
              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;  double p1d = p[2];
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5, 6, 6, 6};
              switch (j) {
                  // d()/dxi
                case 0: return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI7, THIRD, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/deta
                case 1: return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI7, THIRD, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));

                  // d()/dzeta
                case 2: return (MGFE<2,LAGRANGE_XFEM>::shape(TRI7, THIRD, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));

                default: std::cerr<<"Invalid shape function derivative j = " << j;
                }
            }

          // case PYRAMID18: { return fe_fdm_deriv(type, order, i, j, p, fe_lagrange_3D_shape<T>);  }
             // case PYRAMID18: { return fe_fdm_deriv(type, order, i, j, p, fe_lagrange_3D_shape<LAGRANGE_XFEM>);  }
          default: std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
          }
      }
      // unsupported order
    default:  std::cerr<<"ERROR: Unsupported 3D FE order!: " << order; return 0.;
    }


}
// 
// 
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// ==========================================================================
template <FEFamily T>
double fe_lagrange_3D_shape_second_deriv(
  const ElemType type,
  const Order order,
  const unsigned int i,
  const unsigned int j,
  const double  p[]
){//==================================================================================
// #if LIBMESH_DIM == 3
  assert (j< 6);
  switch (order)  {
    case CONSTANT: return 0.;  //ccccccccccccccccccccccccccccccccccccccccccccccccc  
      // linear Lagrange shape functions
    case FIRST: { // ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff->
        switch (type) {
            // Linear tets have all second derivatives = 0
          case TET4:   case TET10:  case TET14:  {    return 0.; }
            // The following elements use either tensor product or
            // rational basis functions, and therefore probably have
            // second derivatives, but we have not implemented them
            // yet...
          case PRISM6:  case PRISM15:   case PRISM18:   case PRISM20:     case PRISM21:     {            assert(i< 6);
              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3]; p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;
              double p1d = p[2];
              //                                0  1  2  3  4  5
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2};
              switch (j) {
                  // All repeated second derivatives and the xi-eta derivative are zero on PRISMs
                case 0: // d^2()/dxi^2
                case 1: // d^2()/dxideta
                case 2: // d^2()/deta^2
                case 5: // d^2()/dzeta^2
                  { return 0.; }
                case 3: // d^2()/dxidzeta
                  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI3,  FIRST, i1[i], 0, p2d)*fe_lagrange_1D_linear_shape_deriv(i0[i], 0, p1d));
                case 4: // d^2()/detadzeta
                  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI3,  FIRST, i1[i], 1, p2d)* fe_lagrange_1D_linear_shape_deriv(i0[i], 0, p1d));
                default:
                  std::cerr<<"Invalid j = " << j;
                }
            }
          case PYRAMID5:  case PYRAMID13:  case PYRAMID14:    case PYRAMID18:    {assert(i< 5);
              const double xi   = p[0];   const double eta  = p[1];    const double zeta = p[2]; const double eps  = 1.e-35;
              switch (j)  {
                  // xi-xi and eta-eta derivatives are all zero for PYRAMID5.
                case 0: // d^2()/dxi^2
                case 2: // d^2()/deta^2
                  return 0.;

                case 1: {// d^2()/dxideta
                     switch (i)  {
                      case 0:   case 2:   return 0.25/(1. - zeta + eps);
                      case 1:    case 3:  return -0.25/(1. - zeta + eps);
                      case 4:  return 0.;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 3:   { double den = (1. - zeta + eps)*(1. - zeta + eps);// d^2()/dxidzeta
                    switch (i) {
                      case 0:  case 2:   return 0.25*eta/den;
                      case 1:  case 3:   return -0.25*eta/den;
                      case 4:  return 0.;
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }

                case 4:  {// d^2()/detadzeta
                    double den = (1. - zeta + eps)*(1. - zeta + eps);
                    switch (i) {
                      case 0:   case 2:   return 0.25*xi/den;
                      case 1:   case 3:   return -0.25*xi/den;
                      case 4:   return 0.;
                      default: std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 5: { // d^2()/dzeta^2
                    double den = (1. - zeta + eps)*(1. - zeta + eps)*(1. - zeta + eps);
                    switch (i)    {
                      case 0: case 2:  return 0.5*xi*eta/den;
                      case 1: case 3:  return -0.5*xi*eta/den;
                      case 4: return 0.;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                default: std::cerr<<"Invalid j = " << j;
                }
            }
            // Trilinear shape functions on HEX8s have nonzero mixed second derivatives
          case HEX8:  case HEX20:   case HEX27:{   assert (i< 8);
              // Compute hex shape functions as a tensor-product
              const double xi   = p[0];  const double eta  = p[1];const double zeta = p[2];
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1};
              switch (j) {
                  // All repeated second derivatives are zero on HEX8
                case 0: // d^2()/dxi^2
                case 2: // d^2()/deta^2
                case 5: // d^2()/dzeta^2
                  {  return 0.; }
                case 1: // d^2()/dxideta
                  return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_linear_shape      (i2[i], zeta));
                case 3: // d^2()/dxidzeta
                  return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_linear_shape      (i1[i], eta)*
                          fe_lagrange_1D_linear_shape_deriv(i2[i], 0, zeta));
                case 4: // d^2()/detadzeta
                  return (fe_lagrange_1D_linear_shape      (i0[i], xi)*
                          fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_linear_shape_deriv(i2[i], 0, zeta));
                default:   std::cerr<<"Invalid j = " << j;
                }
            }

          default:
            std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
          }
      } // <----fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

      // quadratic Lagrange shape functions
    case SECOND:    { // ssssssssssssssssssssssssssssssssssssssssssssssssssssssss---> 
       switch (type)   {   // serendipity hexahedral quadratic shape functions
          case HEX20: { assert(i< 20);
              const double xi   = p[0];  const double eta  = p[1];  const double zeta = p[2];
              // these functions are defined for (x,y,z) in [0,1]^3
              // so transform the locations
              const double x = .5*(xi   + 1.); const double y = .5*(eta  + 1.);   const double z = .5*(zeta + 1.);
              switch(j)              {
                case 0: { // d^2()/dxi^2
                    switch(i)   {
                      case 0:  case 1: return (1. - y) * (1. - z);
                      case 2:  case 3: return y * (1. - z);
                      case 4:  case 5: return (1. - y) * z;
                      case 6:  case 7:  return y * z;
                      case 8: return -2. * (1. - y) * (1. - z);
                      case 10: return -2. * y * (1. - z);
                      case 16:  return -2. * (1. - y) * z;
                      case 18:  return -2. * y * z;
                      case 9:  case 11:
                      case 12:   case 13:   case 14:   case 15:      case 17:   case 19:    return 0.;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 1: {// d^2()/dxidet
                    switch(i) {
                      case 0: return (1.25 - x - y - .5*z) * (1. - z);
                      case 1: return (-x + y + .5*z - .25) * (1. - z);
                      case 2: return (x + y - .5*z - .75) * (1. - z);
                      case 3: return (-y + x + .5*z - .25) * (1. - z);
                      case 4:  return -.25*z * (4.*x + 4.*y - 2.*z - 3);
                      case 5:  return -.25*z * (-4.*y + 4.*x + 2.*z - 1.);
                      case 6:  return .25*z * (-5 + 4.*x + 4.*y + 2.*z);
                      case 7:   return .25*z * (4.*x - 4.*y - 2.*z + 1.);
                      case 8: return (-1. + 2.*x) * (1. - z);
                      case 9:  return (1. - 2.*y) * (1. - z);
                      case 10: return (1. - 2.*x) * (1. - z);
                      case 11:  return (-1. + 2.*y) * (1. - z);
                      case 12: return z * (1. - z);
                      case 13: return -z * (1. - z);
                      case 14: return z * (1. - z);
                      case 15: return -z * (1. - z);
                      case 16: return (-1. + 2.*x) * z;
                      case 17:  return (1. - 2.*y) * z;
                      case 18: return (1. - 2.*x) * z;
                      case 19: return (-1. + 2.*y) * z;
                      default: std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 2: // d^2()/deta^2
                  switch(i) {
                    case 0:   case 3:   return (1. - x) * (1. - z);
                    case 1:  case 2:  return x * (1. - z);
                    case 4:  case 7: return (1. - x) * z;
                    case 5: case 6:   return x * z;
                    case 9: return -2. * x * (1. - z);
                    case 11:   return -2. * (1. - x) * (1. - z);
                    case 17: return -2. * x * z;
                    case 19:  return -2. * (1. - x) * z;
                    case 8:    case 10:    case 12:    case 13:      case 14:    case 15:   case 16:   case 18:  return 0.;
                    default:
                      std::cerr<<"Invalid i = " << i;
                    }
                case 3: // d^2()/dxidzeta
                  switch(i)   {
                    case 0:   return (1.25 - x - .5*y - z) * (1. - y);
                    case 1:   return (-x + .5*y + z - .25) * (1. - y);
                    case 2:   return -.25*y * (2.*y + 4.*x - 4.*z - 1.);
                    case 3:   return -.25*y * (-2.*y + 4.*x + 4.*z - 3);
                    case 4:   return (-z + x + .5*y - .25) * (1. - y);
                    case 5:   return (x - .5*y + z - .75) * (1. - y);
                    case 6:   return .25*y * (2.*y + 4.*x + 4.*z - 5);
                    case 7:   return .25*y * (-2.*y + 4.*x - 4.*z + 1.);
                    case 8:   return (-1. + 2.*x) * (1. - y);
                    case 9:   return -y * (1. - y);
                    case 10:   return (-1. + 2.*x) * y;
                    case 11:   return y * (1. - y);
                    case 12:   return (-1. + 2.*z) * (1. - y);
                    case 13:   return (1. - 2.*z) * (1. - y);
                    case 14:   return (1. - 2.*z) * y;
                    case 15:   return (-1. + 2.*z) * y;
                    case 16:   return (1. - 2.*x) * (1. - y);
                    case 17:   return y * (1. - y);
                    case 18:   return (1. - 2.*x) * y;
                    case 19:   return -y * (1. - y);
                    default:   std::cerr<<"Invalid i = " << i;
                    }
                case 4: // d^2()/detadzeta
                  switch(i)  {
                    case 0:  return (1.25 - .5*x - y - z) * (1. - x);
                    case 1:  return .25*x * (2.*x - 4.*y - 4.*z + 3.);
                    case 2:  return -.25*x * (2.*x + 4.*y - 4.*z - 1.);
                    case 3:  return (-y + .5*x + z - .25) * (1. - x);
                    case 4:  return (-z + .5*x + y - .25) * (1. - x);
                    case 5:  return -.25*x * (2.*x - 4.*y + 4.*z - 1.);
                    case 6:  return .25*x * (2.*x + 4.*y + 4.*z - 5);
                    case 7:  return (y - .5*x + z - .75) * (1. - x);
                    case 8: return x * (1. - x);
                    case 9:  return (-1. + 2.*y) * x;
                    case 10:  return -x * (1. - x);
                    case 11:  return (-1. + 2.*y) * (1. - x);
                    case 12:  return (-1. + 2.*z) * (1. - x);
                    case 13:  return (-1. + 2.*z) * x;
                    case 14:  return (1. - 2.*z) * x;
                    case 15:  return (1. - 2.*z) * (1. - x);
                    case 16:  return -x * (1. - x);
                    case 17:  return (1. - 2.*y) * x;
                    case 18:  return x * (1. - x);
                    case 19:  return (1. - 2.*y) * (1. - x);
                    default:  std::cerr<<"Invalid i = " << i;
                    }
                case 5: // d^2()/dzeta^2
                  switch(i)  {
                    case 0:  case 4:  return (1. - x) * (1. - y);
                    case 1:  case 5:  return x * (1. - y);
                    case 2:  case 6:   return x * y;
                    case 3:   case 7:  return (1. - x) * y;
                    case 12:  return -2. * (1. - x) * (1. - y);
                    case 13:   return -2. * x * (1. - y);
                    case 14:   return -2. * x * y;
                    case 15:  return -2. * (1. - x) * y;
                    case 8:  case 9:   case 10:  case 11:   case 16:   case 17: case 18:  case 19:  return 0.;
                    default:   std::cerr<<"Invalid i = " << i;
                    }
                default: std::cerr<<"Invalid j = " << j;
                }
            }

            // triquadratic hexahedral shape functions
          case HEX8:   assert(T == L2_LAGRANGE);std::cerr<< "High order on first order elements only supported for L2 families"; (void(0));
          case HEX27:{ assert (i <27);
              // Compute hex shape functions as a tensor-product
              const double xi   = p[0]; const double eta  = p[1]; const double zeta = p[2];
              // The only way to make any sense of this
              // is to look at the mgflo/mg2/mgf documentation
              // and make the cut-out cube!
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
              static const unsigned int i0[] = {0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 0, 2, 2, 1, 2, 0, 2, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 1, 2, 0, 0, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2};
              static const unsigned int i2[] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1, 2};

              switch(j){
                  // d^2()/dxi^2
                case 0:   return (fe_lagrange_1D_quadratic_shape_second_deriv(i0[i], 0, xi)*
                                  fe_lagrange_1D_quadratic_shape             (i1[i], eta)*
                                  fe_lagrange_1D_quadratic_shape             (i2[i], zeta));
                  // d^2()/dxideta
                case 1:  return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_quadratic_shape      (i2[i], zeta));
                  // d^2()/deta^2
                case 2:  return (fe_lagrange_1D_quadratic_shape             (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_quadratic_shape             (i2[i], zeta));
                  // d^2()/dxidzeta
                case 3: return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape      (i1[i], eta)*
                          fe_lagrange_1D_quadratic_shape_deriv(i2[i], 0, zeta));
                  // d^2()/detadzeta
                case 4:
                  return (fe_lagrange_1D_quadratic_shape      (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta)*
                          fe_lagrange_1D_quadratic_shape_deriv(i2[i], 0, zeta));
                  // d^2()/dzeta^2
                case 5:
                  return (fe_lagrange_1D_quadratic_shape             (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape             (i1[i], eta)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i2[i], 0, zeta));
                default:  std::cerr<<"Invalid j = " << j;
                }
            }
            // quadratic tetrahedral shape functions
          case TET4:  assert(T == L2_LAGRANGE);std::cerr<<  "High order on first order elements only supported for L2 families";
            (void(0));
          case TET10: case TET14: {
              // The area coordinates are the same as used for the
              // shape() and shape_deriv() functions.
              // const double zeta0 = 1. - zeta1 - zeta2 - zeta3;
              // const double zeta1 = p[0];
              // const double zeta2 = p[1];
              // const double zeta3 = p[2];
              static const double dzetadxi[4][3] =   {
                  {-1., -1., -1.},
                  {1.,   0.,  0.},
                  {0.,   1.,  0.},
                  {0.,   0.,  1.}    };

              // Convert from j -> (j,k) indices for independent variable
              // (0=xi, 1=eta, 2=zeta)
              static const unsigned short int independent_var_indices[6][2] =  {
                  {0, 0}, // d^2 phi / dxi^2
                  {0, 1}, // d^2 phi / dxi deta
                  {1, 1}, // d^2 phi / deta^2
                  {0, 2}, // d^2 phi / dxi dzeta
                  {1, 2}, // d^2 phi / deta dzeta
                  {2, 2}  // d^2 phi / dzeta^2
                };

              // Convert from i -> zeta indices.  Each quadratic shape
              // function for the Tet10 depends on up to two of the zeta
              // area coordinate functions (see the shape() function above).
              // This table just tells which two area coords it uses.
              static const unsigned short int zeta_indices[10][2] =  {
                  {0, 0},
                  {1, 1},
                  {2, 2},
                  {3, 3},
                  {0, 1},
                  {1, 2},
                  {2, 0},
                  {0, 3},
                  {1, 3},
                  {2, 3},
                };

              // Look up the independent variable indices for this value of j.
              const unsigned int my_j = independent_var_indices[j][0];
              const unsigned int my_k = independent_var_indices[j][1];
              if (i<4)  {  return 4.*dzetadxi[i][my_j]*dzetadxi[i][my_k];   }
              else if (i<10)  {
                  const unsigned short int my_m = zeta_indices[i][0];
                  const unsigned short int my_n = zeta_indices[i][1];
                  return 4.*(dzetadxi[my_n][my_j]*dzetadxi[my_m][my_k] +
                             dzetadxi[my_m][my_j]*dzetadxi[my_n][my_k] );
              }
              else std::cerr<<"Invalid shape function index " << i;
            }
            // "serendipity" prism
          case PRISM15:  { assert(i< 15);
              const double xi   = p[0];   const double eta  = p[1];    const double zeta = p[2];
              switch (j)   {
                  // d^2()/dxi^2
                case 0:     {   switch(i)     {
                      case 0:   case 1:   return 2.*(1. - zeta);
                      case 2:   case 5:    case 7:  case 8:  case 9: case 10: case 11: case 13: case 14:  return 0.;
                      case 3:  case 4:    return 2.*(1. + zeta);
                      case 6:   return 4.*(zeta - 1);
                      case 12:  return -4.*(1. + zeta);
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }
                  // d^2()/dxideta
                case 1:                  {         switch(i)       {
                      case 0:  case 7:  return 2.*(1. - zeta);
                      case 1:  case 2:  case 4:   case 5:    case 9:   case 10:  case 11:  return 0.;
                      case 3:  case 13:  return 2.*(1. + zeta);
                      case 6:  case 8: return 2.*(zeta - 1.);
                      case 12:     case 14:  return -2.*(1. + zeta);
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                  // d^2()/deta^2
                case 2:   {      switch(i)  {
                      case 0: case 2:    return 2.*(1. - zeta);
                      case 1:  case 4:  case 6: case 7: case 9:  case 10:   case 11:  case 12:  case 13:  return 0.;
                      case 3:  case 5:   return 2.*(1. + zeta);
                      case 8: return 4.*(zeta - 1.);
                      case 14:  return -4.*(1. + zeta);
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                  // d^2()/dxidzeta
                case 3:  {
                    switch(i)  {
                      case 0:  return 1.5 - zeta - 2.*xi - 2.*eta;
                      case 1:  return 0.5 + zeta - 2.*xi;
                      case 2: case 5: case 11:   return 0.;
                      case 3:  return -1.5 - zeta + 2.*xi + 2.*eta;
                      case 4:  return -0.5 + zeta + 2.*xi;
                      case 6:  return 4.*xi + 2.*eta - 2.;
                      case 7:  return -2.*eta;
                      case 8:  return 2.*eta;
                      case 9:  return 2.*zeta;
                      case 10: return -2.*zeta;
                      case 12: return -4.*xi - 2.*eta + 2.;
                      case 13: return 2.*eta;
                      case 14: return -2.*eta;
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }
                  // d^2()/detadzeta
                case 4:   {
                    switch(i)   {
                      case 0:  return 1.5 - zeta - 2.*xi - 2.*eta;
                      case 1:  case 4:  case 10:  return 0.;
                      case 2:  return .5 + zeta - 2.*eta;
                      case 3: return -1.5 - zeta + 2.*xi + 2.*eta;
                      case 5: return -.5 + zeta + 2.*eta;
                      case 6: return 2.*xi;
                      case 7:  return -2.*xi;
                      case 8:  return 2.*xi + 4.*eta - 2.;
                      case 9:  return 2.*zeta;
                      case 11:   return -2.*zeta;
                      case 12:   return -2.*xi;
                      case 13:   return 2.*xi;
                      case 14:    return -2.*xi - 4.*eta + 2.;
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }
                  // d^2()/dzeta^2
                case 5:    {
                    switch(i)  {
                      case 0:    case 3:  return 1. - xi - eta;
                      case 1:   case 4:   return xi;
                      case 2:    case 5:   return eta;
                      case 6:   case 7:   case 8:   case 12:    case 13:   case 14:  return 0.;
                      case 9:  return 2.*xi + 2.*eta - 2.;
                      case 10:  return -2.*xi;
                      case 11: return -2.*eta;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                default:  std::cerr<<"Invalid j = " << j;
                }
            }
            // quadratic prism shape functions
          case PRISM6: assert(T == L2_LAGRANGE);std::cerr<< "High order on first order elements only supported for L2 families";
            (void(0));
          case PRISM18: case PRISM20:  case PRISM21:  { assert (i< 18);
              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;
              double p1d = p[2];
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5};
              switch (j){
                  // d^2()/dxi^2
                case 0:return (MGFE<2,LAGRANGE_XFEM>::shape_second_deriv(TRI6, SECOND, i1[i], 0, p2d)*fe_lagrange_1D_quadratic_shape(i0[i], p1d));
                  // d^2()/dxideta
                case 1:  return (MGFE<2,LAGRANGE_XFEM>::shape_second_deriv(TRI6, SECOND, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));
                  // d^2()/deta^2
                case 2:   return (MGFE<2,LAGRANGE_XFEM>::shape_second_deriv(TRI6, SECOND, i1[i], 2, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));
                  // d^2()/dxidzeta
                case 3:  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI6,  SECOND, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));
                  // d^2()/detadzeta
                case 4:  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI6,  SECOND, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));
                  // d^2()/dzeta^2
                case 5:  return (MGFE<2,LAGRANGE_XFEM>::shape(TRI6,  SECOND, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i0[i], 0, p1d));
                default:  std::cerr<<"Invalid shape function derivative j = " << j;
                }
            }
            // Quadratic shape functions, as defined in R. Graglia, "Higher order
            // bases on pyramidal elements", IEEE Trans Antennas and Propagation,
            // vol 47, no 5, May 1999.
          case PYRAMID5:  assert(T == L2_LAGRANGE);std::cerr<<  "High order on first order elements only supported for L2 families"; (void(0));
          case PYRAMID14:  case PYRAMID18:{  assert (i< 14);
              const double xi   = p[0];  const double eta  = p[1];  const double zeta = p[2]; const double eps  = 1.e-35;
              // The "normalized coordinates" defined by Graglia.  These are
              // the planes which define the faces of the pyramid.
              double   p1 = 0.5*(1. - eta - zeta), // back
                p2 = 0.5*(1. + xi  - zeta), // left
                p3 = 0.5*(1. + eta - zeta), // front
                p4 = 0.5*(1. - xi  - zeta); // right
              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              double  den = (-1. + zeta + eps),  den2 = den*den,  den3 = den2*den,    den4 = den2*den2;
              // These terms are used in several of the derivatives
              double  numer_mp = xi*eta - zeta + zeta*zeta,   numer_pm = xi*eta + zeta - zeta*zeta;

              switch (j) {
                case 0:   {// d^2()/dxi^2
                    switch(i) {
                      case 0:  case 1: return -p1*eta/den2;
                      case 2:  case 3:   return p3*eta/den2;
                      case 4:   case 9:     case 10:   case 11:   case 12:  return 0.;
                      case 5:  return 2.*p1*eta/den2;
                      case 6:  case 8: return 4.*p1*p3/den2;
                      case 7:  return -2.*p3*eta/den2;
                      case 13:   return -8.*p1*p3/den2;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 1: {// d^2()/dxideta                  
                    switch(i) {
                      case 0:  return 0.25*numer_mp/den2 - 0.5*p1*xi/den2- 0.5*p4*eta/den2 + p4*p1/den2;
                      case 1: return 0.25*numer_pm/den2 - 0.5*p1*xi/den2 + 0.5*p2*eta/den2 - p1*p2/den2;
                      case 2:  return 0.25*numer_mp/den2+ 0.5*p3*xi/den2 + 0.5*p2*eta/den2 + p2*p3/den2;
                      case 3:  return 0.25*numer_pm/den2 + 0.5*p3*xi/den2- 0.5*p4*eta/den2 - p3*p4/den2;
                      case 4:     return 0.;
                      case 5:   return p4*eta/den2 - 2.*p4*p1/den2 - p2*eta/den2  + 2.*p1*p2/den2;
                      case 6:   return -p3*xi/den2 + p1*xi/den2    - 2.*p2*p3/den2+ 2.*p1*p2/den2;
                      case 7:   return p4*eta/den2 + 2.*p3*p4/den2  - p2*eta/den2 - 2.*p2*p3/den2;
                      case 8:   return -p3*xi/den2  + p1*xi/den2 - 2.*p4*p1/den2 + 2.*p3*p4/den2;
                      case 9:                      case 11:       return -zeta/den;
                      case 10:     case 12:        return zeta/den;
                      case 13:   return 4.*p4*p1/den2- 4.*p3*p4/den2+ 4.*p2*p3/den2  - 4.*p1*p2/den2;
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 2:  {// d^2()/deta^2
                    switch(i)   {
                      case 0:  case 3:  return -p4*xi/den2;
                      case 1:  case 2:  return p2*xi/den2;
                      case 4:  case 9:  case 10:  case 11:   case 12:   return 0.;
                      case 5: case 7:   return 4.*p2*p4/den2;
                      case 6: return -2.*p2*xi/den2;
                      case 8:   return 2.*p4*xi/den2;
                      case 13:  return -8.*p2*p4/den2;
                      default: std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 3: // d^2()/dxidzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.25*numer_mp/den2
                          - 0.5*p1*(2.*zeta - 1.)/den2
                          + p1*numer_mp/den3
                          - 0.5*p1*eta/den2
                          - 0.5*p4*eta/den2
                          - 2.*p4*p1*eta/den3;

                      case 1:
                        return 0.25*numer_pm/den2
                          - 0.5*p1*(1 - 2.*zeta)/den2
                          + p1*numer_pm/den3
                          + 0.5*p2*eta/den2
                          + 0.5*p1*eta/den2
                          + 2.*p1*p2*eta/den3;

                      case 2:
                        return -0.25*numer_mp/den2
                          + 0.5*p3*(2.*zeta - 1.)/den2
                          - p3*numer_mp/den3
                          - 0.5*p3*eta/den2
                          - 0.5*p2*eta/den2
                          - 2.*p2*p3*eta/den3;

                      case 3:
                        return -0.25*numer_pm/den2
                          + 0.5*p3*(1 - 2.*zeta)/den2
                          - p3*numer_pm/den3
                          + 0.5*p4*eta/den2
                          + 0.5*p3*eta/den2
                          + 2.*p3*p4*eta/den3;

                      case 4:
                        return 0.;

                      case 5:
                        return p4*eta/den2
                          + 4.*p4*p1*eta/den3
                          - p2*eta/den2
                          - 4.*p1*p2*eta/den3;

                      case 6:
                        return -p3*xi/den2
                          - p1*xi/den2
                          - 4.*p1*p3*xi/den3
                          - 2.*p2*p3/den2
                          - 2.*p1*p3/den2
                          - 2.*p1*p2/den2
                          - 8.*p1*p2*p3/den3;

                      case 7:
                        return -p4*eta/den2
                          - 4.*p3*p4*eta/den3
                          + p2*eta/den2
                          + 4.*p2*p3*eta/den3;

                      case 8:
                        return -p3*xi/den2
                          - p1*xi/den2
                          - 4.*p1*p3*xi/den3
                          + 2.*p4*p1/den2
                          + 2.*p1*p3/den2
                          + 2.*p3*p4/den2
                          + 8.*p3*p4*p1/den3;

                      case 9:
                        return -zeta/den
                          + 2.*p1/den
                          - 2.*p1*zeta/den2;

                      case 10:
                        return zeta/den
                          - 2.*p1/den
                          + 2.*p1*zeta/den2;

                      case 11:
                        return zeta/den
                          - 2.*p3/den
                          + 2.*p3*zeta/den2;

                      case 12:
                        return -zeta/den
                          + 2.*p3/den
                          - 2.*p3*zeta/den2;

                      case 13:
                        return -4.*p4*p1/den2
                          - 4.*p3*p4/den2
                          - 16.*p3*p4*p1/den3
                          + 4.*p2*p3/den2
                          + 4.*p1*p2/den2
                          + 16.*p1*p2*p3/den3;

                      default:
                        std::cerr<<"Invalid i = " << i;
                      }
                  }

                case 4: // d^2()/detadzeta
                  {
                    switch(i)
                      {
                      case 0:
                        return 0.25*numer_mp/den2
                          - 0.5*p4*(2.*zeta - 1.)/den2
                          + p4*numer_mp/den3
                          - 0.5*p1*xi/den2
                          - 0.5*p4*xi/den2
                          - 2.*p4*p1*xi/den3;

                      case 1:
                        return -0.25*numer_pm/den2
                          + 0.5*p2*(1. - 2.*zeta)/den2
                          - p2*numer_pm/den3
                          + 0.5*p2*xi/den2
                          + 0.5*p1*xi/den2
                          + 2.*p1*p2*xi/den3;

                      case 2:
                        return -0.25*numer_mp/den2
                          + 0.5*p2*(2.*zeta - 1.)/den2
                          - p2*numer_mp/den3
                          - 0.5*p3*xi/den2
                          - 0.5*p2*xi/den2
                          - 2.*p2*p3*xi/den3;

                      case 3:
                        return 0.25*numer_pm/den2
                          - 0.5*p4*(1. - 2.*zeta)/den2
                          + p4*numer_pm/den3
                          + 0.5*p4*xi/den2
                          + 0.5*p3*xi/den2
                          + 2.*p3*p4*xi/den3;

                      case 4:
                        return 0.;

                      case 5:
                        return -p4*eta/den2
                          - p2*eta/den2
                          - 4.*p2*p4*eta/den3
                          + 2.*p4*p1/den2
                          + 2.*p2*p4/den2
                          + 2.*p1*p2/den2
                          + 8.*p2*p1*p4/den3;

                      case 6:
                        return p3*xi/den2
                          + 4.*p2*p3*xi/den3
                          - p1*xi/den2
                          - 4.*p1*p2*xi/den3;

                      case 7:
                        return -p4*eta/den2
                          - p2*eta/den2
                          - 4.*p2*p4*eta/den3
                          - 2.*p3*p4/den2
                          - 2.*p2*p4/den2
                          - 2.*p2*p3/den2
                          - 8.*p2*p3*p4/den3;

                      case 8:
                        return p1*xi/den2
                          + 4.*p4*p1*xi/den3
                          - p3*xi/den2
                          - 4.*p3*p4*xi/den3;

                      case 9:
                        return -zeta/den
                          + 2.*p4/den
                          - 2.*p4*zeta/den2;

                      case 10:
                        return -zeta/den
                          + 2.*p2/den
                          - 2.*p2*zeta/den2;

                      case 11:
                        return zeta/den
                          - 2.*p2/den
                          + 2.*p2*zeta/den2;

                      case 12:
                        return zeta/den
                          - 2.*p4/den
                          + 2.*p4*zeta/den2;

                      case 13:
                        return 4.*p3*p4/den2
                          + 4.*p2*p3/den2
                          + 16.*p2*p3*p4/den3
                          - 4.*p4*p1/den2
                          - 4.*p1*p2/den2
                          - 16.*p2*p1*p4/den3;

                      default:      std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 5: { // d^2()/dzeta^2
                    switch(i)   {
                      case 0:   return 0.5*numer_mp/den2 - p1*(2.*zeta - 1.)/den2  + 2.*p1*numer_mp/den3   - p4*(2.*zeta - 1.)/den2
                          + 2.*p4*numer_mp/den3           + 2.*p4*p1/den2     - 4.*p4*p1*(2.*zeta - 1.)/den3  + 6.*p4*p1*numer_mp/den4;
                      case 1:   return -0.5*numer_pm/den2 + p2*(1 - 2.*zeta)/den2 - 2.*p2*numer_pm/den3   + p1*(1 - 2.*zeta)/den2
                          - 2.*p1*numer_pm/den3           + 2.*p1*p2/den2         + 4.*p1*p2*(1 - 2.*zeta)/den3  - 6.*p1*p2*numer_pm/den4;
                      case 2:   return 0.5*numer_mp/den2 - p3*(2.*zeta - 1.)/den2 + 2.*p3*numer_mp/den3- p2*(2.*zeta - 1.)/den2
                          + 2.*p2*numer_mp/den3 + 2.*p2*p3/den2  - 4.*p2*p3*(2.*zeta - 1.)/den3    + 6.*p2*p3*numer_mp/den4;
                      case 3:   return -0.5*numer_pm/den2   + p4*(1 - 2.*zeta)/den2- 2.*p4*numer_pm/den3   + p3*(1 - 2.*zeta)/den2
                          - 2.*p3*numer_pm/den3     + 2.*p3*p4/den2 + 4.*p3*p4*(1 - 2.*zeta)/den3   - 6.*p3*p4*numer_pm/den4;
                      case 4:   return 4.;
                      case 5:  return -2.*p1*eta/den2 - 2.*p4*eta/den2 - 8.*p4*p1*eta/den3 - 2.*p2*eta/den2  - 8.*p2*p4*eta/den3
                          - 8.*p1*p2*eta/den3    - 24.*p2*p1*p4*eta/den4;
                      case 6:  return 2.*p3*xi/den2 + 2.*p2*xi/den2   + 8.*p2*p3*xi/den3  + 2.*p1*xi/den2 + 8.*p1*p3*xi/den3
                          + 8.*p1*p2*xi/den3   + 24.*p1*p2*p3*xi/den4;
                      case 7:  return 2.*p4*eta/den2 + 2.*p3*eta/den2   + 8.*p3*p4*eta/den3      + 2.*p2*eta/den2
                          + 8.*p2*p4*eta/den3  + 8.*p2*p3*eta/den3  + 24.*p2*p3*p4*eta/den4;
                      case 8: return -2.*p1*xi/den2  - 2.*p4*xi/den2 - 8.*p4*p1*xi/den3    - 2.*p3*xi/den2
                          - 8.*p1*p3*xi/den3   - 8.*p3*p4*xi/den3                  - 24.*p3*p4*p1*xi/den4;
                      case 9:    return -2.*zeta/den   + 4.*p4/den    - 4.*p4*zeta/den2   + 4.*p1/den
                          - 4.*p1*zeta/den2            + 8.*p4*p1/den2- 8.*p1*p4*zeta/den3;
                      case 10:  return -2.*zeta/den + 4.*p1/den  - 4.*p1*zeta/den2    + 4.*p2/den
                          - 4.*p2*zeta/den2  + 8.*p1*p2/den2   - 8.*p2*p1*zeta/den3;
                      case 11:   return -2.*zeta/den  + 4.*p2/den  - 4.*p2*zeta/den2  + 4.*p3/den
                          - 4.*p3*zeta/den2           + 8.*p2*p3/den2   - 8.*p3*p2*zeta/den3;
                      case 12: return -2.*zeta/den  + 4.*p3/den    - 4.*p3*zeta/den2  + 4.*p4/den
                          - 4.*p4*zeta/den2   + 8.*p3*p4/den2    - 8.*p4*p3*zeta/den3;
                      case 13:  return 8.*p3*p4/den2 + 8.*p2*p4/den2
                          + 8.*p2*p3/den2  + 32.*p2*p3*p4/den3   + 8.*p4*p1/den2  + 8.*p1*p3/den2
                          + 32.*p3*p4*p1/den3 + 8.*p1*p2/den2  + 32.*p2*p1*p4/den3   + 32.*p1*p2*p3/den3
                          + 96.*p1*p2*p3*p4/den4;
                      default:    std::cerr<<"Invalid i = " << i;
                      }
                  }
                default:  std::cerr<<"Invalid j = " << j;
                }
            }
            // G. Bedrosian, "Shape functions and integration formulas for
            // three-dimensional finite element analysis", Int. J. Numerical
            // Methods Engineering, vol 35, p. 95-108, 1992.
          case PYRAMID13: { assert (i< 13);
              const double xi   = p[0];  const double eta  = p[1];   const double zeta = p[2];  const double eps  = 1.e-35;
              // Denominators are perturbed by epsilon to avoid
              // divide-by-zero issues.
              double  den = (-1. + zeta + eps),  den2 = den*den, den3 = den2*den,   xi2 = xi*xi,
                eta2 = eta*eta,  zeta2 = zeta*zeta,   zeta3 = zeta2*zeta;
              switch (j)  {
                case 0:{ // d^2()/dxi^2
                    switch(i)  {
                      case 0:  case 1:  return 0.5*(-1. + zeta + eta)/den;
                      case 2:  case 3:  return 0.5*(-1. + zeta - eta)/den;
                      case 4:  case 6:   case 8:    case 9:  case 10:  case 11: case 12:   return 0.;
                      case 5:   return (1. - eta - zeta)/den;
                      case 7:   return (1. + eta - zeta)/den;
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 1:   {// d^2()/dxideta                
                    switch(i)  {
                      case 0:  return  0.25*(-1. + 2.*zeta + 2.*xi + 2.*eta)/den;
                      case 1:  return -0.25*(-1. + 2.*zeta - 2.*xi + 2.*eta)/den;
                      case 2:  return -0.25*(1. - 2.*zeta + 2.*xi + 2.*eta)/den;
                      case 3:  return  0.25*(1. - 2.*zeta - 2.*xi + 2.*eta)/den;
                      case 4:  return 0.;
                      case 5:  return -xi/den;
                      case 6:  return eta/den;
                      case 7:  return xi/den;
                      case 8:  return -eta/den;
                      case 9:  return -zeta/den;
                      case 10:  return zeta/den;
                      case 11:  return -zeta/den;
                      case 12: return zeta/den;
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 2:  {// d^2()/deta^2
                    switch(i) {
                      case 0:    case 3:   return 0.5*(-1. + zeta + xi)/den;
                      case 1:    case 2:  return 0.5*(-1. + zeta - xi)/den;
                      case 4:   case 5:   case 7:   case 9:  case 10:   case 11:   case 12: return 0.;
                      case 6:  return (1. + xi - zeta)/den;
                      case 8:  return (1. - xi - zeta)/den;
                      default:   std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 3: {// d^2()/dxidzeta
                    switch(i) {
                      case 0:  return -0.25*(-1. + 2.*zeta - zeta2 + eta + 2.*eta*xi + eta2)/den2;
                      case 1: return 0.25*(-1. + 2.*zeta - zeta2 + eta - 2.*eta*xi + eta2)/den2;
                      case 2:   return 0.25*(-1. + 2.*zeta - zeta2 - eta + 2.*eta*xi + eta2)/den2;
                      case 3:  return -0.25*(-1. + 2.*zeta - zeta2 - eta - 2.*eta*xi + eta2)/den2;
                      case 4:  return 0.;
                      case 5:  return eta*xi/den2;
                      case 6:  return -0.5*(1. + zeta2 + eta2 - 2.*zeta)/den2;
                      case 7:  return -eta*xi/den2;
                      case 8:  return 0.5*(1. + zeta2 + eta2 - 2.*zeta)/den2;
                      case 9: return (-1. - zeta2 + eta + 2.*zeta)/den2;
                      case 10: return -(-1. - zeta2 + eta + 2.*zeta)/den2;
                      case 11: return (1. + zeta2 + eta - 2.*zeta)/den2;
                      case 12: return -(1. + zeta2 + eta - 2.*zeta)/den2;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 4: {// d^2()/detadzeta
                    switch(i) {
                      case 0:  return -0.25*(-1. + 2.*zeta - zeta2 + xi + 2.*eta*xi + xi2)/den2;
                      case 1:  return 0.25*(1. - 2.*zeta + zeta2 + xi + 2.*eta*xi - xi2)/den2;
                      case 2:  return 0.25*(-1. + 2.*zeta - zeta2 - xi + 2.*eta*xi + xi2)/den2;
                      case 3:  return -0.25*(1. - 2.*zeta + zeta2 - xi + 2.*eta*xi - xi2)/den2;
                      case 4:  return 0.;
                      case 5:  return 0.5*(1. + xi2 + zeta2 - 2.*zeta)/den2;
                      case 6:  return -eta*xi/den2;
                      case 7:  return -0.5*(1. + xi2 + zeta2 - 2.*zeta)/den2;
                      case 8:  return eta*xi/den2;
                      case 9:  return (-1. - zeta2 + xi + 2.*zeta)/den2;
                      case 10: return -(1. + zeta2 + xi - 2.*zeta)/den2;
                      case 11: return (1. + zeta2 + xi - 2.*zeta)/den2;
                      case 12: return -(-1. - zeta2 + xi + 2.*zeta)/den2;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                case 5:   {// d^2()/dzeta^2
                    switch(i) {
                      case 0: return 0.5*(xi + eta + 1.)*eta*xi/den3;
                      case 1: return -0.5*(eta - xi + 1.)*eta*xi/den3;
                      case 2: return -0.5*(xi + eta - 1.)*eta*xi/den3;
                      case 3: return 0.5*(eta - xi - 1.)*eta*xi/den3;
                      case 4: return 4.;
                      case 5: return -(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta*xi2)/den3;
                      case 6: return (-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta2*xi)/den3;
                      case 7: return (-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta*xi2)/den3;
                      case 8: return -(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta2*xi)/den3;
                      case 9: return -2.*(-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta*xi)/den3;
                      case 10: return 2.*(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta*xi)/den3;
                      case 11: return -2.*(-1. + 3.*zeta - 3.*zeta2 + zeta3 + eta*xi)/den3;
                      case 12: return 2.*(1. - 3.*zeta + 3.*zeta2 - zeta3 + eta*xi)/den3;
                      default:  std::cerr<<"Invalid i = " << i;
                      }
                  }
                default:  std::cerr<<"Invalid j = " << j;
                }
            }
          default:  std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
          }
      }

    case THIRD:  {
        switch (type) {
          case TET14:  { assert (i< 14); // quadratic Lagrange shape functions with a cubic bubble
              // The area coordinates are the same as used for the
              // shape() and shape_deriv() functions.
              // const double zeta0 = 1. - zeta1 - zeta2 - zeta3;
              // const double zeta1 = p[0];
              // const double zeta2 = p[1];
              // const double zeta3 = p[2];
              static const double dzetadxi[4][3] = {
                  {-1., -1., -1.},
                  {1.,   0.,  0.},
                  {0.,   1.,  0.},
                  {0.,   0.,  1.}
                };
              // Convert from j -> (j,k) indices for independent variable
              // (0=xi, 1=eta, 2=zeta)
              static const unsigned short int independent_var_indices[6][2] =  {
                  {0, 0}, // d^2 phi / dxi^2
                  {0, 1}, // d^2 phi / dxi deta
                  {1, 1}, // d^2 phi / deta^2
                  {0, 2}, // d^2 phi / dxi dzeta
                  {1, 2}, // d^2 phi / deta dzeta
                  {2, 2}  // d^2 phi / dzeta^2
                };
              // Convert from i -> zeta indices.  Each quadratic shape
              // function for the Tet10 depends on up to two of the zeta
              // area coordinate functions (see the shape() function above).
              // This table just tells which two area coords it uses.
              static const unsigned short int zeta_indices[10][2] = {
                  {0, 0},
                  {1, 1},
                  {2, 2},
                  {3, 3},
                  {0, 1},
                  {1, 2},
                  {2, 0},
                  {0, 3},
                  {1, 3},
                  {2, 3},
                };

              // Look up the independent variable indices for this value of j.
              const unsigned int my_j = independent_var_indices[j][0];
              const unsigned int my_k = independent_var_indices[j][1];
              double returnval = 0;
              if (i<4)  returnval = 4.*dzetadxi[i][my_j]*dzetadxi[i][my_k];
              else if (i<10)  {
                  const unsigned short int my_m = zeta_indices[i][0];
                  const unsigned short int my_n = zeta_indices[i][1];
                  returnval =  4.*(dzetadxi[my_n][my_j]*dzetadxi[my_m][my_k] +
                                  dzetadxi[my_m][my_j]*dzetadxi[my_n][my_k] );
                }
              const double zeta1 = p[0];  const double zeta2 = p[1];  const double zeta3 = p[2];
              const double zeta0 = 1. - zeta1 - zeta2 - zeta3;
              // Fill these with whichever derivative we're concerned
              // with
              double d2bubble012, d2bubble013, d2bubble023, d2bubble123;
              switch (j)
                {
                  // d^2()/dxi^2
                case 0:
                  {
                    d2bubble012 = -2.*zeta2;
                    d2bubble013 = -2.*zeta3;
                    d2bubble023 = 0.;
                    d2bubble123 = 0.;
                    break;
                  }

                  // d^2()/dxideta
                case 1:
                  {
                    d2bubble012 = (zeta0-zeta1)-zeta2;
                    d2bubble013 = -zeta3;
                    d2bubble123 = zeta3;
                    d2bubble023 = -zeta3;
                    break;
                  }

                  // d^2()/deta^2
                case 2:
                  {
                    d2bubble012 = -2.*zeta1;
                    d2bubble013 = 0.;
                    d2bubble123 = 0.;
                    d2bubble023 = -2.*zeta3;
                    break;
                  }

                  // d^2()/dxi dzeta
                case 3:
                  {
                    d2bubble012 = -zeta2;
                    d2bubble013 = (zeta0-zeta3)-zeta1;
                    d2bubble123 = zeta2;
                    d2bubble023 = -zeta2;
                    break;
                  }

                  // d^2()/deta dzeta
                case 4:
                  {
                    d2bubble012 = -zeta1;
                    d2bubble013 = -zeta1;
                    d2bubble123 = zeta1;
                    d2bubble023 = (zeta0-zeta3)-zeta2;
                    break;
                  }

                  // d^2()/dzeta^2
                case 5:
                  {
                    d2bubble012 = 0.;
                    d2bubble013 = -2.*zeta1;
                    d2bubble123 = 0.;
                    d2bubble023 = -2.*zeta2;
                    break;
                  }

                default:
                  std::cerr<<"Invalid j = " << j;
                }

              switch (i)
                {
                case 0:
                  return returnval + 3.*(d2bubble012+d2bubble013+d2bubble023);

                case 1:
                  return returnval + 3.*(d2bubble012+d2bubble013+d2bubble123);

                case 2:
                  return returnval + 3.*(d2bubble012+d2bubble023+d2bubble123);

                case 3:
                  return returnval + 3.*(d2bubble013+d2bubble023+d2bubble123);

                case 4:
                  return returnval - 12.*(d2bubble012+d2bubble013);

                case 5:
                  return returnval - 12.*(d2bubble012+d2bubble123);

                case 6:
                  return returnval - 12.*(d2bubble012+d2bubble023);

                case 7:
                  return returnval - 12.*(d2bubble013+d2bubble023);

                case 8:
                  return returnval - 12.*(d2bubble013+d2bubble123);

                case 9:
                  return returnval - 12.*(d2bubble023+d2bubble123);

                case 10:
                  return 27.*d2bubble012;

                case 11:
                  return 27.*d2bubble013;

                case 12:
                  return 27.*d2bubble123;

                case 13:
                  return 27.*d2bubble023;

                default:
                  std::cerr<<"Invalid i = " << i;
                }
            }

          case PRISM20:  {  assert(i< 20);

              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;// Point p2d(p[0],p[1]);
              double p1d = p[2];

              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1};
              const double mainval = MGFE<3,LAGRANGE_XFEM>::shape_second_deriv(PRISM21, THIRD, i, j, p);

              if (i0[i] != 2)     return mainval;
              double bubbleval = 0;

              switch (j)   {
                  // d^2()/dxi^2
                case 0:
                  bubbleval =
                    MGFE<2,LAGRANGE_XFEM>::shape_second_deriv(TRI7, THIRD, 6, 0, p2d)*
                    fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;

                  // d^2()/dxideta
                case 1:
                  bubbleval =
                    MGFE<2,LAGRANGE_XFEM>::shape_second_deriv(TRI7, THIRD, 6, 1, p2d)*
                    fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;

                  // d^2()/deta^2
                case 2:
                  bubbleval =
                    MGFE<2,LAGRANGE_XFEM>::shape_second_deriv(TRI7, THIRD, 6, 2, p2d)*
                    fe_lagrange_1D_quadratic_shape(2, p1d);
                  break;

                  // d^2()/dxidzeta
                case 3:
                  bubbleval =
                    MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI7, THIRD, 6, 0, p2d)*
                    fe_lagrange_1D_quadratic_shape_deriv(2, 0, p1d);
                  break;

                  // d^2()/detadzeta
                case 4:
                  bubbleval =
                    MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI7, THIRD, 6, 1, p2d)*
                    fe_lagrange_1D_quadratic_shape_deriv(2, 0, p1d);
                  break;

                  // d^2()/dzeta^2
                case 5:
                  bubbleval =
                    MGFE<2,LAGRANGE_XFEM>::shape(TRI7, THIRD, 6, p2d)*
                    fe_lagrange_1D_quadratic_shape_second_deriv(2, 0, p1d);
                  break;

                default:
                  std::cerr<<"Invalid shape function derivative j = " << j;
                }

              if (i < 12) // vertices
                return mainval - bubbleval / 9;

              return mainval + bubbleval * (double(4) / 9);
            }

          case PRISM21: {  assert(i< 21);
              // Compute prism shape functions as a tensor-product
              // of a triangle and an edge
              double p2d[3];p2d[0]=p[0];p2d[1]=p[1];p2d[2]=0.;// Point p2d(p[0],p[1]);
              double p1d = p[2];
              //                                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
              static const unsigned int i0[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1, 2, 2, 2, 0, 1, 2};
              static const unsigned int i1[] = {0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5, 6, 6, 6};
              switch (j) {
                case 0:// d^2()/dxi^2
                  return (MGFE<2,LAGRANGE_XFEM>::shape_second_deriv(TRI7, THIRD, i1[i], 0, p2d)*
                                            fe_lagrange_1D_quadratic_shape(i0[i], p1d));
                case 1: // d^2()/dxideta
                  return (MGFE<2,LAGRANGE_XFEM>::shape_second_deriv(TRI7, THIRD, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));
                case 2:// d^2()/deta^2
                  return (MGFE<2,LAGRANGE_XFEM>::shape_second_deriv(TRI7, THIRD, i1[i], 2, p2d)*
                          fe_lagrange_1D_quadratic_shape(i0[i], p1d));
                case 3: // d^2()/dxidzeta
                  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI7, THIRD, i1[i], 0, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));
                case 4:// d^2()/detadzeta
                  return (MGFE<2,LAGRANGE_XFEM>::shape_deriv(TRI7, THIRD, i1[i], 1, p2d)*
                          fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, p1d));
                case 5:// d^2()/dzeta^2
                  return (MGFE<2,LAGRANGE_XFEM>::shape(TRI7, THIRD, i1[i], p2d)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i0[i], 0, p1d));
                default:  std::cerr<<"Invalid shape function derivative j = " << j;
                }
            }

          case PYRAMID18:
            {std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
              // return fe_fdm_second_deriv(type, order, i, j, p, fe_lagrange_3D_shape_deriv<T>);
            }
          default:   std::cerr<<"ERROR: Unsupported 3D element type!: " << type;
          }
      }
      // unsupported order
    default: std::cerr<<"ERROR: Unsupported 3D MGFE order!: " << order; return 0.;
    }

// #else // LIBMESH_DIM != 3
//   libmesh_ignore(type, order, i, j, p);
//   libmesh_not_implemented();
// #endif
}
// 
// #endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
// 
// } // anonymous namespace



template <> double   MGFE<3,LAGRANGE_XFEM>:: compute_fe_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz[], // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
   double InvJac[]
 ){
  
  return 0.;
}

template <> double   MGFE<3,LAGRANGE_XFEM>:: compute_fe_face_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz[] // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
 ){
  
  return 0.;
}
// ========================================
/// This function computes the 3D derivatives at
/// the gaussian point ng:
template <> double MGFE<3,LAGRANGE_XFEM>::Jac(  // 3D  jacobean ->
const int ng,    // gaussian point <-
const double xyz[],    // coordinates <-
double InvJac[]  // inv Jacobean ->
) {                  // =====================================

  double x_xi = 0., x_eta = 0., x_zeta = 0.;
  double y_xi = 0., y_eta = 0., y_zeta = 0.;
  double z_xi = 0., z_eta = 0., z_zeta = 0.;
  int nshape = _NoShape;
  int offset = nshape * _NoGauss;

  for(int s = 0; s < nshape; s++) {
    int sng = s * _NoGauss + ng;
    double dphidx = _dphidxez_map[sng];
    double dphideta = _dphidxez_map[sng + offset];
    double dphidzeta = _dphidxez_map[sng + 2 * offset];

    x_xi += xyz[s] * dphidx;
    x_eta += xyz[s] * dphideta;
    x_zeta += xyz[s] * dphidzeta;
    y_xi += xyz[s + nshape] * dphidx;
    y_eta += xyz[s + nshape] * dphideta;
    y_zeta += xyz[s + nshape] * dphidzeta;
    z_xi += xyz[s + 2 * nshape] * dphidx;
    z_eta += xyz[s + 2 * nshape] * dphideta;
    z_zeta += xyz[s + 2 * nshape] * dphidzeta;
  }

  double det = x_xi * (y_eta * z_zeta - y_zeta * z_eta) - x_eta * (y_xi * z_zeta - y_zeta * z_xi) +
               x_zeta * (y_xi * z_eta - y_eta * z_xi);
  double invdet = 1. / det;
  InvJac[0] = (y_eta * z_zeta - y_zeta * z_eta) * invdet;
  InvJac[3] = -(x_eta * z_zeta - x_zeta * z_eta) * invdet;  //
  InvJac[6] = (y_zeta * x_eta - x_zeta * y_eta) * invdet;
  InvJac[1] = -(y_xi * z_zeta - y_zeta * z_xi) * invdet;  //
  InvJac[4] = (x_xi * z_zeta - x_zeta * z_xi) * invdet;
  InvJac[7] = -(y_zeta * x_xi - y_xi * x_zeta) * invdet;  //
  InvJac[2] = (z_eta * y_xi - y_eta * z_xi) * invdet;
  InvJac[5] = -(x_xi * z_eta - x_eta * z_xi) * invdet;  //
  InvJac[8] = (x_xi * y_eta - y_xi * x_eta) * invdet;

  return (det);
}

// ===============================================
/// This function computes the surface Jacobian at the gaussian point ng:
template <> double MGFE<3,LAGRANGE_XFEM>::JacSur(
const int ng,    // gaussian point <-
const double x[],      // coordinates <-
double InvJac[]  // Jacobean ->
) {        // ======================================
int NDOF_FEMB=9;
  // Values to compute at gaussian points
  double dxdxi = 0.;
  double dxdeta = 0.;  //
  double dydxi = 0.;
  double dydeta = 0.;  // d(x,y,z)d(xi,eta)
  double dzdxi = 0.;
  double dzdeta = 0.;  //

  //   int nshape=_NoShape[1];
  int offset = _NoShape * _NoGauss;

  for(int s = 0; s < (int)_NoShape; s++) {
    int sng = s * _NoGauss + ng;
    double dphidxi = _dphidxez_map[sng];
    double dphideta = _dphidxez_map[sng + offset];

    dxdxi += x[s] * dphidxi;
    dxdeta += x[s] * dphideta;
    dydxi += x[s + NDOF_FEMB] * dphidxi;
    dydeta += x[s + NDOF_FEMB] * dphideta;
    dzdxi += x[s + 2 * NDOF_FEMB] * dphidxi;
    dzdeta += x[s + 2 * NDOF_FEMB] * dphideta;
  }

  // surface weighted jacobean
  double det = sqrt(
               (dxdxi * dydeta - dxdeta * dydxi) * (dxdxi * dydeta - dxdeta * dydxi) +
               (dydxi * dzdeta - dydeta * dzdxi) * (dydxi * dzdeta - dydeta * dzdxi) +
               (dzdxi * dxdeta - dzdeta * dxdxi) * (dzdxi * dxdeta - dzdeta * dxdxi));
  double idet = 1. / det;
  InvJac[0] = dydeta * idet;   // dxi dx
  InvJac[1] = -dydxi * idet;   // deta dx
  InvJac[2] = -dxdeta * idet;  // dxi dy
  InvJac[3] = dxdxi * idet;    // deta dy
  return (det);
}

// ========================================
/// This function computes the 3D derivatives at
/// the gaussian point ng. It computes the Jacobian (f=Jv+b  J=(nabla_v f)^T)
template <> double MGFE<3,LAGRANGE_XFEM>::JacobianOnGauss(  // 3D  jacobean ->
const int dim,   ///< dimension
const int ng,    ///< gaussian point <-
const double xyz[],    ///< coord vect <-
double Jac[]     ///< Jacobean ->
) {                  // =====================================
// Note for c++ array[2][3] = {a11,a12,a13,a21,a22,a23}
//   const int order=dim-1;
  const int nshape = _NoShape;
  const int offset = nshape * _NoGauss;

  // (Jac^T)_ij= d x_j/deta_i
  // (Jac^T)_ij=[d x_j/deta_0,d x_j/deta_1,d x_j/deta_2]
  for(int i=0; i<dim; i++) {
    for(int j=0; j<dim; j++) {
      Jac[i*dim+j]=0;

      for(int s = 0; s < nshape; s++) {
        const int  sng = s * _NoGauss + ng;
        // Transp JacT[i*dim+j]=NablaCoord += xyz[s+j*nshape] *_dphidxez_map1[order][sng+ i * offset];
        Jac[i*dim+j] += xyz[s+i*nshape] *_dphidxez_map[sng+j*offset];
      }
    }
  }

  double det = Jac[0];



 det = Jac[0]*(Jac[4]*Jac[8]-Jac[5]*Jac[7])
                     -Jac[3]*(Jac[1]*Jac[8]-Jac[2]*Jac[7])
                     +Jac[6]*(Jac[1]*Jac[5]-Jac[2]*Jac[4]);

  return (det);
}
// template <> double  MGFE<3,LAGRANGE_XFEM>::FirstDerivateOfLocalPhi(int nPhi, double point[], int Dimension, int DirDer, int FamilyType) {
//   double PhiDer;
// 
//   if(FamilyType == 0) {   // TRIANGULAR ELEMENTS
//     if(_order == 1) {
//         PhiDer = Tri_3d_LinearDerPhi(nPhi, point, DirDer);
//     }
// 
//     if(_order == 2) {
//    
//         PhiDer = Tri_3d_QuadraticDerPhi(nPhi, point, DirDer);
//    
//     }
//   }
// 
//   if(FamilyType == 1) {   // QUADRANGULAR ELEMENTS
//     if(_order == 1) {
//       PhiDer = Rec_Lin_DPhi(nPhi, point, Dimension, DirDer);
//     }
// 
//     if(_order == 2) {
//       PhiDer = Rec_Quad_DPhi(nPhi, point, Dimension, DirDer);
//     }
//   }
// 
//   return PhiDer;
// }
// ========================================
template <> double MGFE<3,LAGRANGE_XFEM>::JacobianOnPoint(
const int dim,const double xyz[], double ElemCoords[], double Jac[], int FamilyType, int nShape) {

//   double LocDPhi[3];
// 
//   for(int dir = 0; dir < dim * dim; dir++) { Jac[dir] = 0.; }
// 
//   for(int s = 0; s < nShape; s++) {  // (s) shape derivative on
//     for(int dir = 0; dir < dim; dir++) {
//       LocDPhi[dir] = FirstDerivateOfLocalPhi(s, xyz, dim, dir, FamilyType);
//     }
// 
//     //  Jac[j * dim + i]_i,j= dx_j/deta_i
//     for(int j = 0; j < dim; j++) {  // field
//       for(int i = 0; i < dim; i++) {  // direction
//         Jac[j * dim + i] += ElemCoords[s + j * nShape] * LocDPhi[i];
//       }
//     }
//   }
// 
//   double det = 0.;
// 
//     det =   Jac[0]*(Jac[4] * Jac[8] - Jac[5] * Jac[7])
//             -Jac[3]*(Jac[1] * Jac[8] - Jac[2] * Jac[7])
//             +Jac[6]*(Jac[1] * Jac[5] - Jac[2] * Jac[4]);
// 
// //     for(int dir = 0; dir < dim; dir++) {
// //       int sign = 1 - 2 * (dir % 2);
// //       int idx1 = (dir + 1) % dim;
// //       int idx2 = (dir + 2) % dim;
// //       det += sign *JacT[dir * dim + dir]* (JacT[idx1 * dim + idx1] * JacT[idx2 * dim + idx2] -
// //                      JacT[idx1 * dim + idx2] * JacT[idx2 * dim + idx1]);
// //     }
// //   }
//   
// 
//   return det;
  return 0.;
}



template <> double MGFE<3,LAGRANGE_XFEM>::shape(const ElemType type,
                           const Order order,
                           const unsigned int i,
                           const double p[])
{
  return fe_lagrange_3D_shape<LAGRANGE_XFEM>(type, order, i, p);
}
// // TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
// // Basic TRI/TETRA
// // TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
// 
// // =======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Tri_2d_LinearPhi(
//   int nPhi,
//   double point[]
// ) {// ====================================================================================
//   double PhiVal;
//   int a = _CooTriEl[nPhi * 3];
//   int b = _CooTriEl[nPhi * 3 + 1];
//   int c = _CooTriEl[nPhi * 3 + 2];
//   PhiVal = a * point[0] + b * point[1] + c;
//   return PhiVal;
// }
// // =======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Tri_2d_QuadraticPhi(
//   int nPhi,
//   double point[]
// ) {// ====================================================================================
//   double PhiVal;
//   int m = nPhi / 3;
//   int l1 = (nPhi % 3);
//   int l2 = (l1 + 1) % 3;
//   double L1, L2;
//   L1 = Tri_2d_LinearPhi(l1, point);
//   L2 = Tri_2d_LinearPhi(l2, point);
//   PhiVal = L1 * ((1 - m) * (2. * L1 - 1) + m * 4. * L2);
//   return PhiVal;
// }
// // =======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Tri_2d_LinearDerPhi(int nPhi, double */*point[]*/, int dir) {
//   double PhiDer;
//   PhiDer = _CooTriEl[nPhi * 3 + dir];
//   return PhiDer;
// }
// // =======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Tri_2d_QuadraticDerPhi(int nPhi, double point[], int dir) {
//   double PhiDer;
//   int m = nPhi / 3;
//   int l1 = (nPhi % 3);
//   int l2 = (l1 + 1) % 3;
//   double L1, L2, dL1, dL2;
//   L1 = Tri_2d_LinearPhi(l1, point);
//   L2 = Tri_2d_LinearPhi(l2, point);
//   dL1 = Tri_2d_LinearDerPhi(l1, point, dir);
//   dL2 = Tri_2d_LinearDerPhi(l2, point, dir);
// 
//   PhiDer = dL1 * ((1 - m) * (2. * L1 - 1) + m * 4. * L2) + L1 * ((1 - m) * (2. * dL1) + m * 4. * dL2);
//   return PhiDer;
// }
// // =======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Tri_2d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2) {
//   double Phi2Der;
//   int m = nPhi / 3;
//   int l1 = (nPhi % 3);
//   int l2 = (l1 + 1) % 3;
//   double dL1_2, dL2_2, dL1_1, dL2_1;
// 
//   dL1_1 = Tri_2d_LinearDerPhi(l1, point, dir1);
//   dL2_1 = Tri_2d_LinearDerPhi(l2, point, dir1);
//   dL1_2 = Tri_2d_LinearDerPhi(l1, point, dir2);
//   dL2_2 = Tri_2d_LinearDerPhi(l2, point, dir2);
// 
//   Phi2Der =
//     dL1_1 * ((1 - m) * (2. * dL1_2) + m * 4. * dL2_2) + dL1_2 * ((1 - m) * (2. * dL1_1) + m * 4. * dL2_1);
//   return Phi2Der;
// }

// // QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
// // Basic QUAD/HEX
// // QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
// // =======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Edge_Quad_Phi(int PhiCoeff, double Coordinate) {
//   double PhiVal;
//   PhiVal = (1 - 0.5 * fabs(PhiCoeff)) * ((2 * fabs(PhiCoeff) - 1) * Coordinate * Coordinate +
//       PhiCoeff * Coordinate + (1 - fabs(PhiCoeff)));
//   return PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Edge_Quad_DPhi(int PhiCoeff, double Coordinate) {
//   double DPhiVal;
//   DPhiVal = (1 - 0.5 * fabs(PhiCoeff)) * (2. * (2 * fabs(PhiCoeff) - 1) * Coordinate + PhiCoeff);
//   return DPhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Edge_Quad_D2Phi(int PhiCoeff, double /*Coordinate*/) {
//   double D2PhiVal;
//   D2PhiVal = (1 - 0.5 * fabs(PhiCoeff)) * (2. * (2 * fabs(PhiCoeff) - 1));
//   return D2PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Edge_Lin_Phi(int PhiCoeff, double Coordinate) {
//   double PhiVal;
//   PhiVal = 0.5 * (1 + PhiCoeff * Coordinate);
//   return PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Edge_Lin_DPhi(int PhiCoeff, double /*Coordinate*/) {
//   double DPhiVal;
//   DPhiVal = 0.5 * PhiCoeff;
//   return DPhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Rec_Lin_Phi(int nPhi, double point[], int dimension) {
//   double PhiVal = 1.;
// 
//   for(int dir = 0; dir < dimension; dir++) {
//     int lambda_i;
// 
//     if(dimension == 1) {
//       lambda_i = _CooE3[nPhi];
//     } else if(dimension == 2) {
//       lambda_i = _CooQ9[nPhi + dir * _Q9Off];
//     } else if(dimension == 3) {
//       lambda_i = _CooH27[nPhi + dir * _H27Off];
//     } else {
//       printf("\033[1;31m MGFE1::Rec_Lin_Phi unkown phi for dimension %d \n \033[0m", dimension);
//     }
// 
//     PhiVal *= Edge_Lin_Phi(lambda_i, point[dir]);
//   }
// 
//   return PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Rec_Quad_Phi(int nPhi, double point[], int dimension) {
//   double PhiVal = 1.;
// 
//   for(int dir = 0; dir < dimension; dir++) {
//     int lambda_i;
// 
//     if(dimension == 1) {
//       lambda_i = _CooE3[nPhi];
//     } else if(dimension == 2) {
//       lambda_i = _CooQ9[nPhi + dir * _Q9Off];
//     } else if(dimension == 3) {
//       lambda_i = _CooH27[nPhi + dir * _H27Off];
//     } else {
//       printf("\033[1;31m MGFE1::Rec_Quad_Phi unkown phi for dimension %d \n \033[0m", dimension);
//     }
// 
//     PhiVal *= Edge_Quad_Phi(lambda_i, point[dir]);
//   }
// 
//   return PhiVal;
// }
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Rec_Lin_DPhi(int nPhi, double point[], int dimension, int DirDer) {
//   int lambda_i= _CooH27[nPhi + DirDer * _H27Off];
//  
// 
//   double PhiDer = Edge_Lin_DPhi(lambda_i, point[DirDer]);
// 
//   for(int dir = DirDer + 1; dir < dimension + DirDer; dir++) {
//     int direction = dir % dimension;
//     int lambda_i = _CooH27[nPhi + direction * _H27Off];
//     PhiDer *= Edge_Lin_Phi(lambda_i, point[direction]);
//   }
// 
//   return PhiDer;
// }
// 
// 
// 
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Rec_Quad_DPhi(int nPhi, double point[], int dimension, int DirDer) {
//   int lambda_i;
// 
//   if(dimension == 1) {
//     lambda_i = _CooE3[nPhi];
//   } else if(dimension == 2) {
//     lambda_i = _CooQ9[nPhi + DirDer * _Q9Off];
//   } else if(dimension == 3) {
//     lambda_i = _CooH27[nPhi + DirDer * _H27Off];
//   } else {
//     printf("\033[1;31m MGFE1::Rec_Quad_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//   }
// 
//   double PhiDer = Edge_Quad_DPhi(lambda_i, point[DirDer]);
// 
//   for(int dir = DirDer + 1; dir < dimension + DirDer; dir++) {
//     int lambda_i;
//     int direction = dir % dimension;
// 
//     if(dimension == 1) {
//       lambda_i = _CooE3[nPhi];
//     } else if(dimension == 2) {
//       lambda_i = _CooQ9[nPhi + direction * _Q9Off];
//     } else if(dimension == 3) {
//       lambda_i = _CooH27[nPhi + direction * _H27Off];
//     } else {
//       printf("\033[1;31m MGFE1::Rec_Quad_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//     }
// 
//     PhiDer *= Edge_Quad_Phi(lambda_i, point[direction]);
//   }
// 
//   return PhiDer;
// }
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Rec_Lin_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2) {
//   double PhiDer;
// 
//   if(DirDer1 == DirDer2) {
//     PhiDer = 0.;
//   } else {
//     PhiDer = 1;
// 
//     for(int dir = 0; dir < dimension; dir++) {
//       int lambda_i;
//       int direction = dir % dimension;
// 
//       if(dimension == 1) {
//         lambda_i = _CooE3[nPhi];
//       } else if(dimension == 2) {
//         lambda_i = _CooQ9[nPhi + direction * _Q9Off];
//       } else if(dimension == 3) {
//         lambda_i = _CooH27[nPhi + direction * _H27Off];
//       } else {
//         printf("\033[1;31m MGFE1::Rec_Lin_D2Phi unkown phi for dimension %d \n \033[0m", dimension);
//       }
// 
//       if(dir == DirDer1 || dir == DirDer2) {
//         PhiDer *= Edge_Lin_DPhi(lambda_i, point[dir]);
//       } else {
//         PhiDer *= Edge_Lin_Phi(lambda_i, point[dir]);
//       }
//     }
//   }
// 
//   return PhiDer;
// }
// // ======================================================================================
// template <> double  MGFE<3,LAGRANGE_XFEM>::Rec_Quad_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2) {
//   double PhiDer;
// 
//   if(DirDer1 == DirDer2) {
//     int lambda_i;
// 
//     if(dimension == 1) {
//       lambda_i = _CooE3[nPhi];
//     } else if(dimension == 2) {
//       lambda_i = _CooQ9[nPhi + DirDer1 * _Q9Off];
//     } else if(dimension == 3) {
//       lambda_i = _CooH27[nPhi + DirDer1 * _H27Off];
//     }
// 
//     PhiDer = Edge_Quad_D2Phi(lambda_i, point[DirDer1]);
// 
//     for(int dir = DirDer1 + 1; dir < dimension + DirDer1; dir++) {
//       int lambda_i;
//       int direction = dir % dimension;
// 
//       if(dimension == 1) {
//         lambda_i = _CooE3[nPhi];
//       } else if(dimension == 2) {
//         lambda_i = _CooQ9[nPhi + direction * _Q9Off];
//       } else if(dimension == 3) {
//         lambda_i = _CooH27[nPhi + direction * _H27Off];
//       } else {
//         printf("\033[1;31m MGFE1::Rec_Quad_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//       }
// 
//       PhiDer *= Edge_Quad_Phi(lambda_i, point[direction]);
//     }
//   } else {
//     PhiDer = 1;
// 
//     for(int dir = 0; dir < dimension; dir++) {
//       int lambda_i;
// 
//       if(dimension == 1) {
//         lambda_i = _CooE3[nPhi];
//       } else if(dimension == 2) {
//         lambda_i = _CooQ9[nPhi + dir * _Q9Off];
//       } else if(dimension == 3) {
//         lambda_i = _CooH27[nPhi + dir * _H27Off];
//       } else {
//         printf("\033[1;31m MGFE1::Rec_Quad_D2Phi unkown phi for dimension %d \n \033[0m", dimension);
//       }
// 
//       if(dir == DirDer1 || dir == DirDer2) {
//         PhiDer *= Edge_Quad_DPhi(lambda_i, point[dir]);
//       } else {
//         PhiDer *= Edge_Quad_Phi(lambda_i, point[dir]);
//       }
//     }
//   }
// 
//   return PhiDer;
// }


// ========================================
/// This function computes the 3D derivatives at
/// the nodal point ng:
template <> double  MGFE<3,LAGRANGE_XFEM>::Jac_nodes(  // 3D  jacobean ->
const int ng,          // nodal point <-
const double xyz[],          // coordinates <-
double InvJac[]        // Jacobean ->
) {                        // =====================================

  double x_xi = 0., x_eta = 0., x_zeta = 0.;
  double y_xi = 0., y_eta = 0., y_zeta = 0.;
  double z_xi = 0., z_eta = 0., z_zeta = 0.;
  int nshape = _NoShape;
  int offset = nshape * _NoGauss;

  for(int s = 0; s < nshape; s++) {
    int sng = s * _NoGauss + ng;
    double dphidx = _dphidxez_map_nodes[sng];
    double dphideta = _dphidxez_map_nodes[sng + offset];
    double dphidzeta = _dphidxez_map_nodes[sng + 2 * offset];

    x_xi += xyz[s] * dphidx;
    x_eta += xyz[s] * dphideta;
    x_zeta += xyz[s] * dphidzeta;
    y_xi += xyz[s + nshape] * dphidx;
    y_eta += xyz[s + nshape] * dphideta;
    y_zeta += xyz[s + nshape] * dphidzeta;
    z_xi += xyz[s + 2 * nshape] * dphidx;
    z_eta += xyz[s + 2 * nshape] * dphideta;
    z_zeta += xyz[s + 2 * nshape] * dphidzeta;
  }

  double det = x_xi * (y_eta * z_zeta - y_zeta * z_eta) - x_eta * (y_xi * z_zeta - y_zeta * z_xi) +
               x_zeta * (y_xi * z_eta - y_eta * z_xi);
  double invdet = 1. / det;
  InvJac[0] = (y_eta * z_zeta - y_zeta * z_eta) * invdet;
  InvJac[3] = -(x_eta * z_zeta - x_zeta * z_eta) * invdet;  //
  InvJac[6] = (y_zeta * x_eta - x_zeta * y_eta) * invdet;
  InvJac[1] = -(y_xi * z_zeta - y_zeta * z_xi) * invdet;  //
  InvJac[4] = (x_xi * z_zeta - x_zeta * z_xi) * invdet;
  InvJac[7] = -(y_zeta * x_xi - y_xi * x_zeta) * invdet;  //
  InvJac[2] = (z_eta * y_xi - y_eta * z_xi) * invdet;
  InvJac[5] = -(x_xi * z_eta - x_eta * z_xi) * invdet;  //
  InvJac[8] = (x_xi * y_eta - y_xi * x_eta) * invdet;

  return (det);
}


// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// // Basic piecewise continuous functions
// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// 
// /// /// This function generates the Lagrangian piecewise shape functions
// template <> void MGFE<3,LAGRANGE_XFEM>::init_pie_rec() {  // ================================
// 
//   //   CONTROLLARE PER IMPLEMENTAZIONE ELEMENTI MISTI
//   //
//   //
//   //
//   //               ********************************************
//   //                               S_0
//   //               ********************************************
//   //
//   //                               ___________
//   //                              |           |
//   //                              |     0     |
//   //                              |     x     |
//   //                              |           |
//   //                              |___________|
//   //
// 
// 
// 
// //Order nodes ---------------------------------------------------------------------
//   // Remark _CooH27[27 * 3] =
//   //       // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
//   //       -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1, 1,
//   //       1,  -1, 0,  1,  0,  -1, 0,  0,  1,  0,  -1, 0,  0,  /// xi
//   //       -1, -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1,
//   //       1,  1,  -1, 0,  1,  0,  0,  -1, 0,  1,  0,  0,  0,  /// eta
//   //       -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1, 0,  0,
//   //       0,  0,  1,  1,  1,  1,  -1, 0,  0,  0,  0,  1,  0  /// chi
//   //  QUAD9 nodes note _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)]
//   //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//   int hex_shift=0;
//   // ---------------------------------------------------------------------------------------
//   _fem_conn[0]=_NoElNodes-1;
// 
//   // dofs  ---------------------------------------------------------------------------
//   for(int j = 0; j < _NoElNodes; j++)
//     for(int dir = 0; dir < _dim; dir++) _xnodes[j+_NoElNodes*dir]= _CooH27[hex_shift+j+27*dir];
// 
//   // ---------------------------------------------------------------------------------------
// 
//   // nodes ----------------------------------------------------------------------------
//   double *pointNode=new double[_dim];
// 
//   for(int i = 0; i < _NoElNodes; i++) {  // points
//     for(int dir = 0; dir < _dim; dir++)  pointNode[dir] =   _xnodes[i+_NoElNodes*dir];
// 
//     for(int j = 0; j < _NoShape; j++) {  // shapse
//       //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//       _phi_map_nodes[j*_NoElNodes+i]   =1.;
// 
//       // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =0;     //
//       }
//     }
//   }
//   delete[]pointNode;
//   // ---------------------------------------------------------------------------------------
//   // Gauss integration points (rectangular coordinates) ---------------------------------------------------
//   // basic gaussian points coordinates
//   const double a = -sqrt(3. / 5.);  const double b = 0.;  const double c = -a;
//   const double x27[27*3] = {a, a, a, a, a, a, a, a, a, b, b, b, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c,
//       a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c,
//       a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c
//     };
//   // gaussian weights 1D
//   double weight1[3];  weight1[0] = 5. / 9.;  weight1[1] = 8. / 9.;  weight1[2] = 5. / 9.;
//   double *point=new double[_dim];  int iw[3];
// 
//   for(int i = 0; i < _NoGauss; i++) {
//     iw[0]=i%3;    iw[1]= (i/3)%3;    iw[2]= i/9;    _weight[i] =1.;
// 
//     for(int dir=0; dir<_dim; dir++) {
//       _weight[i] *= weight1[iw[dir]];        // gausssian weights
//       _xg[i+dir*_NoGauss] = x27[i+(3-_dim+dir)*27];
//       point[dir] =  _xg[i+dir*_NoGauss];  //gaussian points
//     }
// 
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       // phi[phi_1(xg1),phi_1(xg2),phi_1(xg3),...phi_2(xg1),phi_2(xg2),....]
//       phi[j*_NoGauss+i] = 1.; // shape functions
// 
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] = 0. ;
// 
//         for(int dir2 = 0; dir2 < _dim; dir2++)
//           _dphidxx_map[(j + (dir * _dim + dir2) *  _NoShape) * _NoGauss + i] =0.;
// 
// //           Rec_Quad_D2Phi(j, point, _dim, dir, dir2);  // d2/dxdy gaussian points
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   }
//   delete[]point;
//   return;
// }
// 
// /// /// This function generates the Lagrangian piecewise shape functions
// template <> void MGFE<3,LAGRANGE_XFEM>::init_pie_tri() {  // ================================
// 
//   /* //        ********************************************
//     //                               TRI 3
//     //     ********************************************
//     //                              2
//     //        |\
//     //        | \
//     //        |  \
//     //   r=l1 |   \ 1-s-r=l0
//     //        |    \
//     //        |     \
//     //        |______\
//     //        0  s=l2 1
//     */
//   // gaussian coordinates
//   // 1D --------------------------------
//   const double x[3] = {-sqrt(3. / 5.), 0., sqrt(3. / 5.)};
//   // 2D --------------------------------
//   const double alphabeta[4] = {1. - 2. * (2. / 7. + sqrt(15.) / 21.), 1. - 2. * (2. / 7. - sqrt(15.) / 21.),
//       2. / 7. + sqrt(15.) / 21., 2. / 7. - sqrt(15.) / 21.
//     };
//   const double r[7] = {1. / 3.,      alphabeta[2], alphabeta[0], alphabeta[2],
//       alphabeta[3], alphabeta[1], alphabeta[3]
//     };
//   const double s[7] = {1. / 3.,      alphabeta[2], alphabeta[2], alphabeta[0],
//       alphabeta[3], alphabeta[3], alphabeta[1]
//     };
// 
// 
// 
// //   int n_nodes=NDOF_FEM;
//   /*  if(_dim==1) {n_nodes=3;} //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//     if(_dim==2) {n_nodes=NDOF_FEMB;}//  QUAD9  _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)] */
//   const double x10[11*3] = {0,1,0,0,  .5,.5,0.,0.,0.,.5,  0.33333333333333333333333,
//       0,0,1,0,  0.,.5,.5,0.,.5,0.,  0.33333333333333333333333,
//       0,0,0,1,  0.,0.,0.,.5,.5,.5,  0.33333333333333333333333
//     };
// // 10 point = quad tetra  + center point
// 
// 
//   // ****LINEAR SHAPES AND DERIVATIVES****
//   // shape 2D in triangular coordinates              derivatives
//   //                                      |_______r_______|_______s_______
//   // phi0= 1-r-s                          |      -1.      |      -1.
//   // phi1= r                              |       1.      |       0.
//   // phi2= s                              |       0.      |       1.
//   // -----------------------------------------------------------------------
//   //                            1D ELEMENT
//   // -----------------------------------------------------------------------
// 
// 
//   // -----------------------------------------------------------------------
//   //                            2D ELEMENT
//   // -----------------------------------------------------------------------
//  
// 
//   }
// // =============================================================================
// template <> void MGFE<3,LAGRANGE_XFEM>::init_pie(
// ) {  // ========================================================================
//   if(_GeomType==1) init_pie_rec();  if(_GeomType==0) init_pie_tri();
//   return;
// }

// OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// Basic piecewise linear functions
// OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO





// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// // Basic Raviart Thomas non continuous functions
// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// //===============================================================================================
// template <> 
// std::array<double, 3 >  MGFE<3,LAGRANGE_XFEM>::Rec_RT_Phi(int iface, double x[],int dim) {
// 
//   
//   switch(iface) {
//   case 0:
//     return {0., 0.,-0.125*(1.-x[2])}; 
//     break;
//   case 1:
//     return {0.,-0.125*(1.-x[1]), 0.};    break;
//   case 2:
//     return {0.125*(x[0]+1.),0.,0.};    break;
//   case 3:
//     return {0.,0.125*(1.+x[1]), 0.};    break;
//   case 4:
//     return {-0.125*(1.-x[0]),0., 0.};    break;
//   case 5:
//     return {0.,0.,0.125*(1.+x[2])};    break;
//   default:  std::cout<< "MGFE1::RT_quad Error iface ";return {0.,0.,0.};exit(3);
//   }
//   return {0.,0.,0.};
// }
// 
// 
// //===============================================================================================
// 
// template <> 
// std::array<double,3>  MGFE<3,LAGRANGE_XFEM>::Rec_RT_DPhi(int iface, double x[], int dim, int dir) { 
//   if(dir==0) {
//     if(iface==2 || iface==4) {
//       return {0.125,0.,0.};
//     } else return {0.,0.,0.};
//   }
// 
//   if(dir==1) {
//     if(iface==1 || iface==3) {
//       return {0.,0.125,0.};
//     } else return {0.,0.,0.};
//   }
// 
//   if(dir==2) {
//     if(iface==0 || iface==5) {
//       return {0., 0.,0.125};
//     } else return {0.,0.,0.};
//   }
// }

// // =======================================================
// template <> void MGFE<3,LAGRANGE_XFEM>::init_lin_rec_rt() {
// 
//   //  ********************************************
//   //                 HEX RT 6
//   //  ********************************************
//   //    The \p Hex8 is an element in 3D composed of 8 nodes.
//   //    It is numbered like this:
//   //
//   //    HEXRT 6   ^
//   //             -|------
//   //           /: |     /|         zeta
//   //          / : o5 3 / |          ^   eta
//   //         /  :   o /  |          | /
//   //         --------   2|          |/
//   //     <--|4o .....|.o-|->        o---> xi
//   //        |  .1o   |  /
//   //        | . / o0 | /
//   //        |. /  |  |/
//   //         -----|--
//   //        0        1
//   //    (xi, eta, zeta) are the reference element coordinates associated with the given numbering.
//   //    ********************************************
//   
//   // ****LINEAR SHAPES AND DERIVATIVES****
//   // shape 2D                                            derivatives
//   //                                           |_______x_______|_______y_______
//   // phi0= (1.-x)*(1.-y)                       |     y.-1      |     x.-1
//   // phi1= x*(1.-y)                            |     1.-y      |      -x
//   // phi2= x*y                                 |       y       |       x
//   // phi3= y*(1.-x)                            |      -y       |     1.-x
// 
//   // shape 2D
//   // phi0 = {0,  -0.25*(1.-y)};
//   // phi1 = { 0.25*(1.+x), 0};
//   // phi2 = {0,   0.25*(1.+y)};
//   // phi1 = { -0.25*(1.-x), 0};
// 
// 
// //Order nodes ---------------------------------------------------------------------
//   // Remark _CooH27[27 * 3] =
//   //       // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
//   //       -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1, 1,
//   //       1,  -1, 0,  1,  0,  -1, 0,  0,  1,  0,  -1, 0,  0,  /// xi
//   //       -1, -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1,
//   //       1,  1,  -1, 0,  1,  0,  0,  -1, 0,  1,  0,  0,  0,  /// eta
//   //       -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1, 0,  0,
//   //       0,  0,  1,  1,  1,  1,  -1, 0,  0,  0,  0,  1,  0  /// chi
//   //  QUAD9 nodes note _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)]
//   //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//   int hex_shift=0;
//   int exx[6];
// 
//   
//     hex_shift= 0;    //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//     exx[0]=2;    exx[1]=1;    exx[2]=0;    exx[3]=1;    exx[4]=0;    exx[5]=2;
// 
//   int NDOF_FACE=6;
//  
// 
//   // ---------------------------------------------------------------------------------------
//   for(int j = 0; j < _NoShape; j++) _fem_conn[j]=j+_NoElNodes-NDOF_FACE-1;
// 
//   // dofs  ---------------------------------------------------------------------------
//   for(int j = 0; j < _NoElNodes; j++)
//     for(int dir = 0; dir < _dim; dir++) _xnodes[j+_NoElNodes*dir]= _CooH27[hex_shift+j+27*dir];
// 
//   // ---------------------------------------------------------------------------------------
// 
//   // nodes ----------------------------------------------------------------------------
//   double *pointNode=new double[_dim];
// 
//   for(int i = 0; i < _NoElNodes; i++) {  // points
//     for(int dir = 0; dir < _dim; dir++)  pointNode[dir] = _xnodes[i+_NoElNodes*dir];
// 
//     for(int j = 0; j < _NoShape; j++) {  // shapse
//       //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//       _phi_map_nodes[j*_NoElNodes+i]=Rec_RT_Phi(j, pointNode,_dim)[exx[j]];
// 
//       // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         // [N^x_1(x0),N^x_1(x1),...N^x_2(x0),N^x_2(x1),...N^x_1(x0),
//         //  N^y_1(x1),...N^y_2(x0),N^y_2(x1),...                    ]
//         _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =Rec_RT_Phi(j, pointNode,_dim)[dir];     //
// //         _dphidxez_map_nodes[(j + dir *_NoShape)*n_nodes+i] =  Rec_RT_DPhi(j, pointNode,_dim, dir);     // dphi/dxi
//       }
//     }
//   }
// 
//   delete[]pointNode;
//   // ---------------------------------------------------------------------------------------
// 
//   // Gauss integration points (rectangular coordinates) ---------------------------------------------------
//   // basic gaussian points coordinates
//   const double a = -sqrt(3. / 5.);
//   const double b = 0.;
//   const double c = -a;
//   const double x27[27*3] = {a, a, a, a, a, a, a, a, a, b, b, b, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c,
//       a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c,
//       a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c
//     };
//   // gaussian weights 1D
//   double weight1[3];
//   weight1[0] = 5. / 9.;
//   weight1[1] = 8. / 9.;
//   weight1[2] = 5. / 9.;
// 
//   double *point=new double[_dim];
//   int iw[3];
// 
//   for(int i = 0; i < _NoGauss; i++) {
//     iw[0]=i%3;    iw[1]= (i/3)%3;    iw[2]= i/9;    _weight[i] =1.;
// 
//     for(int dir=0; dir<_dim; dir++) {
//       _weight[i] *= weight1[iw[dir]];        // gausssian weights
//       _xg[i+dir*_NoGauss] = x27[i+(3-_dim+dir)*27];
//       point[dir] =  _xg[i+dir*_NoGauss];  //gaussian points
//     }
// 
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       // phi[phi_1(xg1),phi_1(xg2),phi_1(xg3),...phi_2(xg1),phi_2(xg2),....]
//       phi[j*_NoGauss+i] = Rec_RT_Phi(j, point,_dim)[exx[j]] ; // shape functions
// 
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//           // [N^x_1(xg0),N^x_1(xg1),...N^x_2(xg0),N^x_2(xg1),...N^x_1(xg0),
//         //  N^y_1(xg1),...N^y_2(xg0),N^y_2(xg1),...                    ]
//         _dphidxez_map[(j + dir *_NoShape)*_NoGauss+i] =Rec_RT_Phi(j, point,_dim)[dir];     //
// //         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] = Rec_RT_DPhi(j, point,_dim,dir)[exx[j]] ;
// 
//         for(int dir2 = 0; dir2 < _dim; dir2++)
//           _dphidxx_map[(j + (dir * _dim + dir2) *_NoShape) * _NoGauss + i] =0.;
// 
// //           Rec_Quad_D2Phi(j, point, _dim, dir, dir2);  // d2/dxdy gaussian points
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   }
// 
//   delete[]point;
//   //--------------------------------------------------------------------------------
// 
//   return;
// }


// // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
// template<> void MGFE<3,LAGRANGE_XFEM>::init_lin_rec(
// 
// ) {// LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
//   //  ********************************************
//   //                 HEX 8
//   //  ********************************************
//   //    The \p Hex8 is an element in 3D composed of 8 nodes.
//   //    It is numbered like this:
//   //
//   //    HEX8: 7        6
//   //            o--------z
//   //           /:       /|         zeta
//   //          / :      / |          ^   eta (into page)
//   //       4 /  :   5 /  |          | /
//   //        o--------o   |          |/
//   //        |   o....|...o 2        o---> xi
//   //        |  .3    |  /
//   //        | .      | /
//   //        |.       |/
//   //        o--------o
//   //        0        1
//   //    (xi, eta, zeta) are the reference element coordinates associated with the given numbering.
//   //    ********************************************
//   
//   //   // ****LINEAR SHAPES AND DERIVATIVES****
//   //   // shape 2D                                            derivatives
//   //   //                                           |_______x_______|_______y_______
//   //   // phi0= (1.-x)*(1.-y)                       |     y.-1      |     x.-1
//   //   // phi1= x*(1.-y)                            |     1.-y      |      -x
//   //   // phi2= x*y                                 |       y       |       x
//   //   // phi3= y*(1.-x)                            |      -y       |     1.-x
//   //
//   //Order nodes ---------------------------------------------------------------------
//   // Remark _CooH27[27 * 3] =
//   //       // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
//   //       -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1,
//   //     -1, 1,  1,  -1, 0,  1,  0,  -1, 0,  0,  1,  0,  -1, 0,  0,  /// xi
//   //       -1, -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1,
//   //       1,  1,  -1, 0,  1,  0,  0,  -1, 0,  1,  0,  0,  0,  /// eta
//   //       -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1, 0,  0,
//   //       0,  0,  1,  1,  1,  1,  -1, 0,  0,  0,  0,  1,  0  /// chi
//   //  QUAD9 nodes note _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)]
//   //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//   int hex_shift=0;
//   // ---------------------------------------------------------------------------------------
//   // dofs  ---------------------------------------------------------------------------
//   for(int j = 0; j <  _NoShape; j++) _fem_conn[j]=j;
// 
//   // nodes ----------------------------------------------------------------------------
//   for(int j = 0; j < _NoElNodes-1; j++) for(int dir=0; dir<_dim; dir++) _xnodes[j+_NoElNodes*dir]=_CooH27[hex_shift+j+27*dir];
// 
//   for(int dir = 0; dir < _dim; dir++)_xnodes[_NoElNodes-1+_NoElNodes*dir]=0.;  // central point
// 
//   // ---------------------------------------------------------------------------------------
// 
//   // nodes ----------------------------------------------------------------------------
//   for(int j = 0; j < _NoElNodes; j++)for(int dir = 0; dir < _dim; dir++) _xnodes[j+_NoElNodes*dir]= _CooH27[hex_shift+j+27*dir];
// 
//   for(int dir = 0; dir < _dim; dir++) _xnodes[_NoElNodes-1+_NoElNodes*dir]=0.;
// 
//   double *pointNode=new double[_dim];
// 
//   for(int i = 0; i < _NoElNodes; i++) {  // points
//     for(int dir = 0; dir < _dim; dir++)  pointNode[dir] = _xnodes[i+_NoElNodes*dir];
// 
//     for(int j = 0; j < _NoShape; j++) {  // shapse
//       //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//       _phi_map_nodes[j*_NoElNodes+i]   = Rec_Lin_Phi(j, pointNode,_dim);
// 
//       // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] = Rec_Lin_DPhi(j,pointNode,_dim,dir);     // dphi/dxi
//       }
//     }
//   }
// 
//   delete[]pointNode; // ---------------------------------------------------------------------------------------
// 
//   // Gauss integration points (rectangular coordinates) ---------------------------------------------------
//   // basic gaussian points coordinates
//   const double a = -sqrt(3. / 5.);
//   const double b = 0.;
//   const double c = -a;
//   const double x27[27*3] = {a, a, a, a, a, a, a, a, a, b, b, b, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c,
//       a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c,
//       a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c
//     };
//   // gaussian weights 1D
//   double weight1[3];
//   weight1[0] = 5. / 9.;
//   weight1[1] = 8. / 9.;
//   weight1[2] = 5. / 9.;
// 
//   double *point=new double[_dim];
//   int iw[3];
// 
//   for(int i = 0; i < _NoGauss; i++) {
//     iw[0]=i%3;
//     iw[1]= (i/3)%3;
//     iw[2]= i/9;
//     _weight[i] =1.;
// 
//     for(int dir=0; dir<_dim; dir++) {
//       _weight[i] *= weight1[iw[dir]];        // gausssian weights
//       _xg[i+dir*_NoGauss] = x27[i+(3-_dim+dir)*27];
//       point[dir] =  _xg[i+dir*_NoGauss];  //gaussian points
//     }
// 
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       // phi[phi_1(xg1),phi_1(xg2),phi_1(xg3),...phi_2(xg1),phi_2(xg2),....]
//       phi[j*_NoGauss+i] = Rec_Lin_Phi(j, point, _dim); // shape functions
// 
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map[(j+dir*_NoShape)*_NoGauss+i] = Rec_Lin_DPhi(j, point, _dim, dir);
// 
//         for(int dir2 = 0; dir2 < _dim; dir2++)
//           _dphidxx_map[(j+(dir*_dim+dir2)*_NoShape)*_NoGauss+i] =0.;
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   }
// 
//   delete[]point;
//   //--------------------------------------------------------------------------------
// 
// 
//   return;
// }
// // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
// template<> void MGFE<3,LAGRANGE_XFEM>::init_lin_tri(
// ) {// LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
//   /* //        ********************************************
//   //                               TRI 3
//   //     ********************************************
//   //                              2
//   //        |\
//   //        | \
//   //        |  \
//   //   r=l1 |   \ 1-s-r=l0
//   //        |    \
//   //        |     \
//   //        |______\
//   //        0  s=l2 1
//   */
//   // gaussian coordinates
//  
//   // 3D --------------------------------
//   const double p[14] = {0.31088591926330060980,
//       0.31088591926330060980,
//       1. - 3. * 0.31088591926330060980,
//       0.31088591926330060980,
//       0.092735250310891226402,
//       0.092735250310891226402,
//       1. - 3. * 0.092735250310891226402,
//       0.092735250310891226402,
//       0.5 - 0.045503704125649649492,
//       0.5 - 0.045503704125649649492,
//       0.045503704125649649492,
//       0.045503704125649649492,
//       0.5 - 0.045503704125649649492,
//       0.045503704125649649492
//     };
//   const double q[14] = {
//     0.31088591926330060980,  1. - 3. * 0.31088591926330060980, 0.31088591926330060980,
//     0.31088591926330060980,  0.092735250310891226402,          1. - 3. * 0.092735250310891226402,
//     0.092735250310891226402, 0.092735250310891226402,          0.5 - 0.045503704125649649492,
//     0.045503704125649649492, 0.045503704125649649492,          0.5 - 0.045503704125649649492,
//     0.045503704125649649492, 0.5 - 0.045503704125649649492
//   };
//   const double t[14] = {0.31088591926330060980,        0.31088591926330060980,
//       0.31088591926330060980,        1. - 3. * 0.31088591926330060980,
//       0.092735250310891226402,       0.092735250310891226402,
//       0.092735250310891226402,       1. - 3. * 0.092735250310891226402,
//       0.045503704125649649492,       0.045503704125649649492,
//       0.5 - 0.045503704125649649492, 0.045503704125649649492,
//       0.5 - 0.045503704125649649492, 0.5 - 0.045503704125649649492
//     };
// 
// 
// //   int n_nodes=NDOF_FEM;
//   /*  if(_dim==1) {n_nodes=3;} //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
//     if(_dim==2) {n_nodes=NDOF_FEMB;}//  QUAD9  _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)] */
//   const double x10[11*3] = {0,1,0,0,  .5,.5,0.,0.,0.,.5,  0.33333333333333333333333,
//       0,0,1,0,  0.,.5,.5,0.,.5,0.,  0.33333333333333333333333,
//       0,0,0,1,  0.,0.,0.,.5,.5,.5,  0.33333333333333333333333
//     };
// // 10 point = quad tetra  + center point
// 
// 
//   // ****LINEAR SHAPES AND DERIVATIVES****
//   // shape 2D in triangular coordinates              derivatives
//   //                                      |_______r_______|_______s_______
//   // phi0= 1-r-s                          |      -1.      |      -1.
//   // phi1= r                              |       1.      |       0.
//   // phi2= s                              |       0.      |       1.
//   // -----------------------------------------------------------------------
//   //                            1D ELEMENT
//   // -----------------------------------------------------------------------
//  
// 
//   // -----------------------------------------------------------------------
//   //                            3D ELEMENT
//   // -----------------------------------------------------------------------
//   // if(_dim == 3) {
// 
//     // dofs  ---------------------------------------------------------------------------
// 
//     for(int j = 0; j <  _NoShape; j++) _fem_conn[j]=j;
//     // nodes ----------------------------------------------------------------------------
//     for(int j = 0; j < _NoElNodes-1; j++)  for(int dir = 0; dir < _dim; dir++) _xnodes[j+_NoElNodes*dir]= x10[j+11*dir];
//     for(int dir = 0; dir < _dim; dir++) _xnodes[_NoElNodes-1+_NoElNodes*dir]=1./3.;
// 
//     double *pointNode=new double[_dim];
//     for(int i = 0; i <  _NoElNodes; i++) {  // points
//       for(int dir = 0; dir < _dim; dir++)  pointNode[dir] = _xnodes[i+_NoElNodes*dir];
// 
//       for(int j = 0; j < _NoShape; j++) {  // shapse
//         //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//         _phi_map_nodes[j*_NoElNodes+i]   = Tri_3d_QuadraticPhi(j, pointNode);
// 
//         // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//         for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//           _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =  Tri_3d_QuadraticDerPhi(j, pointNode, dir);     // dphi/dxi
//         }
//       }
//     // }
// 
//     delete[]pointNode;
// 
// 
// 
// 
// 
// 
//     // gaussian weights
//     const double weight_3[14] = {0.018781320953002641800,  0.018781320953002641800,  0.018781320953002641800,
//         0.018781320953002641800,  0.012248840519393658257,  0.012248840519393658257,
//         0.012248840519393658257,  0.012248840519393658257,  0.0070910034628469110730,
//         0.0070910034628469110730, 0.0070910034628469110730, 0.0070910034628469110730,
//         0.0070910034628469110730, 0.0070910034628469110730
//       };
// 
//     for(int i = 0; i < _NoGauss; i++) {
//       _weight[i] = weight_3[i];
//     }
// 
//     double point[3];
// 
//     for(int nPhi = 0; nPhi < _NoShape; nPhi++) {
//       for(int i = 0; i < _NoGauss; i++) {
//         point[0] = p[i];
//         point[1] = q[i];
//         point[2] = t[i];
//         phi[i + (nPhi)*_NoGauss] = Tri_3d_LinearPhi(nPhi, point);
// 
//         for(int dir1 = 0; dir1 < _dim; dir1++) {
//           _dphidxez_map[i + (nPhi + dir1 * _NoShape) * _NoGauss] = Tri_3d_LinearDerPhi(nPhi, point, dir1);
//         }
//       }
//     }
//   }
// 
//   return;
// }

// // ================================================================================================
// ///  This function generates the Lagrangian linear shape functions
// template <> void MGFE<3,LAGRANGE_XFEM>::init_lin() {  // ================================
// 
//   switch(_GeomType) {
//   case 1: // quad/hex ---------------------
//     if(_FamType==0)  init_lin_rec();
//     if(_FamType==1){      init_lin_rec_rt();    }
//     if(_FamType>1) {      std::cout<<" MGFE1::init_lin error: _FamType";      exit(1);    }   break;
// 
//   case 0:   // tri/tetra ---------------------
//     if(_FamType==0)  init_lin_tri();
//     else {      std::cout<<" MGFE1::init_lin error: _FamType";     exit(3);    }    break;
//   default:     std::cout<<" MGFE1::init_lin error: _GeomType";    exit(3);
//   }
//   return;
// }
// ===========================================================================================
// ===========================================================================================



// ================================================================================================
//    QUAD
// ============================================================================
// template<> void  MGFE<3,LAGRANGE_XFEM>::init_qua_rec(
// ) {//==========================================================================
// ///  The \p Hex27 is an element in 3D composed of 27 nodes.
// ///  It is numbered like this:
// ///
// ///    HEX27:     7              18             6
// ///               o--------------o--------------o
// ///              /:             /              /|
// ///             / :            /              / |
// ///            /  :           /              /  |
// ///         19/   :        25/            17/   |
// ///          o--------------o--------------o    |
// ///         /     :        /              /|    |
// ///        /    15o       /    23o       / |  14o
// ///       /       :      /              /  |   /|           zeta
// ///     4/        :   16/             5/   |  / |            ^   eta (into page)
// ///     o--------------o--------------o    | /  |            | /
// ///     |         :    |   26         |    |/   |            |/
// ///     |  24o    :    |    o         |  22o    |            o---> xi
// ///     |         :    |       10     |   /|    |
// ///     |        3o....|.........o....|../.|....o
// ///     |        .     |              | /  |   / 2
// ///     |       .    21|            13|/   |  /
// ///  12 o--------------o--------------o    | /
// ///     |     .        |              |    |/
// ///     |  11o         | 20o          |    o
// ///     |   .          |              |   / 9
// ///     |  .           |              |  /
// ///     | .            |              | /
// ///     |.             |              |/
// ///     o--------------o--------------o
// ///     0              8              1
// ///
// ///  (xi, eta, zeta) are the reference element coordinates associated with
// ///  the given numbering.
// //     ********************************************
// //                               QUAD 9
// //     ********************************************
// ///
// ///The \p Quad9 is an element in 3D composed of 27 nodes.
// /// It is numbered like this:
// ///
// /// Quad9
// ///  *       3 ______6_____ 2
// ///  *        |            |               ^   eta
// ///  *        |            |               |
// ///  *       7|      8     |5              |
// ///  *        |            |               o---> xi
// ///  *        |____________|
// ///  *       0       4      1
// /// * (xi, eta) are the reference element coordinates associated with
// /// * the given numbering.
// ///
// ///
// //     ********************************************
// //                               EDGE 3
// //     ********************************************
// ///
// ///  *        |____________|     o---> xi
// ///  *       0       2      1
//   /// * (xi) is the reference element coordinate associated with
// /// * the given numbering.
// ///
//   /// General formula for test function: N_i = (1-0.5*|xi_i|)*((2|xi_i|-1)*xi*xi + xi_i*xi + (1-|xi_i|)) for
// /// xi, eta, chi xi is the gaussian point and xi_i is given by _CooE3[i] for the 1D case Loop i is over the
// /// shape functions. Loop j is over the gaussian points
// 
//   //Order nodes ---------------------------------------------------------------------
//   // Remark _CooH27[27 * 3] =
//   //       // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
//   //       -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1, 1,
//   //       1,  -1, 0,  1,  0,  -1, 0,  0,  1,  0,  -1, 0,  0,  /// xi
//   //       -1, -1, 1,  1,  -1, -1, 1,  1,  -1, 0,  1,  0,  -1, -1,
//   //       1,  1,  -1, 0,  1,  0,  0,  -1, 0,  1,  0,  0,  0,  /// eta
//   //       -1, -1, -1, -1, 1,  1,  1,  1,  -1, -1, -1, -1, 0,  0,
//   //       0,  0,  1,  1,  1,  1,  -1, 0,  0,  0,  0,  1,  0  /// chi
//   //  QUAD9 nodes note _CooQ9[i]=[_CooH27[i+12](xi),_CooH27[27+i+12](eta)]
//   //  EDGE3  _CooE3[i]=[_CooH27[i+6]]
// 
//   int hex_shift=  0;
//     // Fem order
//   // =============================
//   
//   
//   
//   //  ------------------------------------------------------------------------
//   // gauss points (_NoGauss) ------------------------------------------------------
//   _NoGauss=0;  assert(_qrule); _NoGauss=_qrule->n_points();
//   if(_NoGauss==0) {   std::cout << "Quadrature not assigned"; abort();}
// 
//   // number of nodes (_NoElNodes) in the reference quad element
//    ElemType  elem_name=_fe_type.fetype; 
//    Order elem_order= _fe_type.order ;
//   // number of nodes (_NoElNodes) in the reference EDGE3 element
//   _NoElNodes= get_refspace_n_nodes(elem_name);//
//   _NoShape=_NoElNodes;  
//   
// //     TIME(std::clock_t start_time = std::clock();)     //  TC +++++++++++++++
// // dofs  ---------------------------------------------------------------------------
// // dofs -----------------------------------------------------------------------
//  for(int j = 0; j < _NoShape; j++) _fem_conn[j]=j;
//  double *xnodes=new double [_NoElNodes*_dim];
//  double *point=new double[_dim];
//     //
//     // ---------------------------------------------------------------------------------------
//     // nodes ----------------------------------------------------------------------------
//   get_refspace_nodes(elem_name,xnodes);// nodes[(x,y)_0,(x,y)_1.... ]
//     // _xnodes[0]=0.;    _xnodes[1]=1.;    _xnodes[2]=0.5;
//     for(int i = 0; i < _NoElNodes; i++) {  // points
//       for(int dir = 0; dir < _dim; dir++) { 
//         _xnodes[i+_NoElNodes*dir]=xnodes[dir+_dim*i]; 
//         point[dir] = _xnodes[i+_NoElNodes*dir];
//       }
//       
//       
//       
//   // for(int j = 0; j <  _NoShape; j++) _fem_conn[j]=j;
//   // 
//   // // nodes ----------------------------------------------------------------------------
//   // for(int j = 0; j < _NoElNodes-1; j++)
//   //   for(int dir = 0; dir < _dim; dir++) _xnodes[j+_NoElNodes*dir]= _CooH27[hex_shift+j+27*dir];
//   // 
//   // for(int dir = 0; dir < _dim; dir++)_xnodes[_NoElNodes-1+_NoElNodes*dir]=0.;
// 
// 
//   // double *pointNode=new double[_dim];
// 
//   // for(int i = 0; i <  _NoElNodes; i++) {  // points
//     // for(int dir = 0; dir < _dim; dir++)  pointNode[dir] = _xnodes[i+_NoElNodes*dir];
// 
//     for(int j = 0; j < _NoShape; j++) {  // shapse
//       //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//       _phi_map_nodes[j*_NoElNodes+i]   = fe_lagrange_3D_shape<LAGRANGE_XFEM>(elem_name , elem_order, j, point);// Rec_Quad_Phi(j, pointNode,_dim);
// 
//       // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(elem_name, elem_order, j,dir, point);// Rec_Quad_DPhi(j,pointNode,_dim,dir);     // dphi/dxi
//       }
//     }
//   }
// 
//   // delete[]pointNode;
//  delete[]xnodes;
// // // basic gaussian points coordinates  ---------------------------------------------------------------------
// //   double x1g[3];  x1g[0] =-sqrt(3. / 5.);  x1g[1] =0.;  x1g[2] =-x1g[0];  
// //   const double a = -sqrt(3. / 5.);  const double b = 0.;  const double c = -a;
// //   const double x27[27*3] = {a, a, a, a, a, a, a, a, a, b, b, b, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c,
// //       a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c, a, a, a, b, b, b, c, c, c,
// //       a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c, a, b, c
// //     };
//   // gaussian weights 1D
//   // double weight1[3];  weight1[0] = 5. / 9.;  weight1[1] = 8. / 9.;  weight1[2] = 5. / 9.;
//   // double *point=new double[_dim]; // int iw[3];
// 
//   for(int i = 0; i < _NoGauss; i++) {
//     _weight[i] =_qrule->w(i); // weight
//      _qrule->qp(point,i);     // ith-gaussian point ->point(i)
//     // iw[0]=i%3;    iw[1]= (i/3)%3;    iw[2]= i/9;    _weight[i] =1.;
// 
//     // for(int dir=0; dir<_dim; dir++) {
//     //   // _weight[i] *= weight1[iw[dir]];        // gausssian weights
//     //   _xg[i+dir*_NoGauss] = x27[i+(3-_dim+dir)*27];
//     //   point[dir] =  _xg[i+dir*_NoGauss];  //gaussian points
//     // }
// 
//     for(int j = 0; j < _NoShape; j++) {     // loop over test function id
//       // phi[phi_1(xg1),phi_1(xg2),phi_1(xg3),...phi_2(xg1),phi_2(xg2),....]
//       phi[j*_NoGauss+i] = fe_lagrange_3D_shape<LAGRANGE_XFEM>(elem_name , elem_order, j, point);//Rec_Quad_Phi(j, point, _dim); // shape functions
// 
//       for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//         _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] = fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(elem_name, elem_order, j,dir, point);//Rec_Quad_DPhi(j, point, _dim, dir);
// 
//         for(int dir2 = 0; dir2 < _dim; dir2++)
//           _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =fe_lagrange_3D_shape_second_deriv<LAGRANGE_XFEM>(elem_name, elem_order,j,dir +dir2, point);//       Rec_Quad_D2Phi(j, point, _dim, dir, dir2);  // d2/dxdy gaussian points
//       }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     }
//   } //--------------------------------------------------------------------------------
//   delete[]point;
//   return;
// }


// // ================================================================================
// template<> void  MGFE<3,LAGRANGE_XFEM>::init_qua_tri() {
//   /*           ********************************************
//   //                               TRI 6
//   //     ********************************************
//   //                              2
//   //        |\
//   //        | \
//   //        |  \
//   //         r=l1 |   \ 1-s-r=l0
//   //        |    \
//   //        |     \
//   //              |______\
//   //             0  s=l2 1
//   */
//   // gaussian coordinates
//   // 1D --------------------------------
//   const double x[3] = {-sqrt(3. / 5.), 0., sqrt(3. / 5.)};
// 
//   // 2D --------------------------------
//   const double alphabeta[4] = {1. - 2. * (2. / 7. + sqrt(15.) / 21.), 1. - 2. * (2. / 7. - sqrt(15.) / 21.),
//       2. / 7. + sqrt(15.) / 21., 2. / 7. - sqrt(15.) / 21.
//     };
//   const double r[7] = {1. / 3.,      alphabeta[2], alphabeta[0], alphabeta[2],
//       alphabeta[3], alphabeta[1], alphabeta[3]
//     };
//   const double s[7] = {1. / 3.,      alphabeta[2], alphabeta[2], alphabeta[0],
//       alphabeta[3], alphabeta[3], alphabeta[1]
//     };
// 
//   // 3D --------------------------------
//   const double p[14] = {0.31088591926330060980,
//       0.31088591926330060980,
//       1. - 3. * 0.31088591926330060980,
//       0.31088591926330060980,
//       0.092735250310891226402,
//       0.092735250310891226402,
//       1. - 3. * 0.092735250310891226402,
//       0.092735250310891226402,
//       0.5 - 0.045503704125649649492,
//       0.5 - 0.045503704125649649492,
//       0.045503704125649649492,
//       0.045503704125649649492,
//       0.5 - 0.045503704125649649492,
//       0.045503704125649649492
//     };
//   const double q[14] = {
//     0.31088591926330060980,  1. - 3. * 0.31088591926330060980, 0.31088591926330060980,
//     0.31088591926330060980,  0.092735250310891226402,          1. - 3. * 0.092735250310891226402,
//     0.092735250310891226402, 0.092735250310891226402,          0.5 - 0.045503704125649649492,
//     0.045503704125649649492, 0.045503704125649649492,          0.5 - 0.045503704125649649492,
//     0.045503704125649649492, 0.5 - 0.045503704125649649492
//   };
//   const double t[14] = {0.31088591926330060980,        0.31088591926330060980,
//       0.31088591926330060980,        1. - 3. * 0.31088591926330060980,
//       0.092735250310891226402,       0.092735250310891226402,
//       0.092735250310891226402,       1. - 3. * 0.092735250310891226402,
//       0.045503704125649649492,       0.045503704125649649492,
//       0.5 - 0.045503704125649649492, 0.045503704125649649492,
//       0.5 - 0.045503704125649649492, 0.5 - 0.045503704125649649492
//     };
// 
//   const double x10[11*3] = {0,1,0,0,  .5,.5,0.,0.,0.,.5,  0.33333333333333333333333,
//       0,0,1,0,  0.,.5,.5,0.,.5,0.,  0.33333333333333333333333,
//       0,0,0,1,  0.,0.,0.,.5,.5,.5,  0.33333333333333333333333
//     };
//   // ****QUADRATIC SHAPES AND DERIVATIVES****
//   // shape 2D in triangular coordinates                  derivatives
//   //                                           |_______r_______|_______s_______
//   // phi0= 2.*(1-r-s)*((1.-r-s)-0.5)           | -3.+4.s+4.r   |  -3.+4.*s+4.*r
//   // phi1= 2.*r*(r-0.5)                        |    4.*r-1     |      0.
//   // phi2= 2.*s*(s-0.5)                        |    0.         |   4.*s-1
//   // phi3= 4.*r*(1.-r-s)                       |  4.-4.*s-8.*r |   -4.*r
//   // phi4= 4.*r*s                              |   4.*s        |   4.*r
//   // phi5= 4.*s*(1.-r-s)                       |  -4.*s        |   4-4.*r-8.*s
// 
//   // -----------------------------------------------------------------------
//   //                            1D ELEMENT
//   // -----------------------------------------------------------------------
// 
// 
// 
//     // dofs  ---------------------------------------------------------------------------
// 
//     for(int j = 0; j <  _NoShape; j++) _fem_conn[j]=j;
//     // nodes ----------------------------------------------------------------------------
//     for(int j = 0; j < _NoElNodes-1; j++)  for(int dir = 0; dir < _dim; dir++) _xnodes[j+_NoElNodes*dir]= x10[j+11*dir];
// 
//     double *pointNode=new double[_dim];
//     for(int i = 0; i <  _NoElNodes; i++) {  // points
//       for(int dir = 0; dir < _dim; dir++)  pointNode[dir] = _xnodes[i+_NoElNodes*dir];
// 
//       for(int j = 0; j < _NoShape; j++) {  // shapse
//         //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
//         _phi_map_nodes[j*_NoElNodes+i]   = Tri_3d_QuadraticPhi(j, pointNode);
// 
//         // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
//         for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
//           _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =  Tri_3d_QuadraticDerPhi(j, pointNode, dir);     // dphi/dxi
//         }
//       }
//     }
// 
//     delete[]pointNode;
// 
// 
// 
// 
// 
// 
//     //  const double weight_3[5]={-16./120., 9./120., 9./120., 9./120., 9./120.};
//     const double weight_3[14] = {0.018781320953002641800,  0.018781320953002641800,  0.018781320953002641800,
//         0.018781320953002641800,  0.012248840519393658257,  0.012248840519393658257,
//         0.012248840519393658257,  0.012248840519393658257,  0.0070910034628469110730,
//         0.0070910034628469110730, 0.0070910034628469110730, 0.0070910034628469110730,
//         0.0070910034628469110730, 0.0070910034628469110730
//       };
// 
//     for(int i =0; i<_NoGauss; i++) {
//       _weight[i]=weight_3[i];
//     }
// 
//     // shape functions
// 
//     for(int nPhi = 0; nPhi < _NoShape; nPhi++) {
//       for(int i = 0; i < _NoGauss; i++) {
//         double point[3];
//         point[0] = p[i];
//         point[1] = q[i];
//         point[2] = t[i];
//         phi[i + (nPhi)*_NoGauss] = Tri_3d_QuadraticPhi(nPhi, point);
// 
//         for(int dir1 = 0; dir1 < _dim; dir1++) {
//           _dphidxez_map[i+(nPhi+dir1*_NoShape)*_NoGauss] = Tri_3d_QuadraticDerPhi(nPhi, point, dir1);
// 
//           for(int dir2=0; dir2<_dim; dir2++)
//             _dphidxx_map[i+(nPhi+(dir2+dir1*_dim)*_NoShape)*_NoGauss] =
//               Tri_3d_QuadraticDer2Phi(nPhi,point,dir1,dir2);
//         }
//       }
//     }
//   
// 
//   return;
// }


// // ================================================================================================
// /// This function generates the Lagrangian quad shape functions
// template<> void  MGFE<3,LAGRANGE_XFEM>::init_qua() {  // ================================
//   if(_GeomType==1) init_qua_rec();
//   if(_GeomType==0) init_qua_tri();
//   return;
// }

// ================================================================================================
// ================================================================================================
// ================================================================================================
//===============================================================================================

template <> double MGFE<3,LAGRANGE_XFEM>::ComputeInverseMatrix_vol(double Matrix[], double InvMatrix[]) {
  // Matrix  = 0  3   6
  //           1  4   7
  //           2  5   8
  double det  = Matrix[0]*(Matrix[4]*Matrix[8]-Matrix[5]*Matrix[7])+
    Matrix[1]*(Matrix[6]*Matrix[5]-Matrix[3]*Matrix[8])+
    Matrix[2]*(Matrix[3]*Matrix[7]-Matrix[4]*Matrix[6]);
  double invdet=1./det;
  // _InvMatrix  = 0  3   6
  //               1  4   7
  //               2  5   8
  InvMatrix[0]= (Matrix[8]*Matrix[4]-Matrix[7]*Matrix[5])*invdet;
  InvMatrix[1]= (Matrix[6]*Matrix[5]-Matrix[8]*Matrix[3])*invdet;
  InvMatrix[2]= (Matrix[7]*Matrix[3]-Matrix[6]*Matrix[4])*invdet;
  InvMatrix[3]= (Matrix[7]*Matrix[2]-Matrix[8]*Matrix[1])*invdet;
  InvMatrix[4]= (Matrix[8]*Matrix[0]-Matrix[6]*Matrix[2])*invdet;
  InvMatrix[5]= (Matrix[6]*Matrix[1]-Matrix[7]*Matrix[0])*invdet;
  InvMatrix[6]= (Matrix[5]*Matrix[1]-Matrix[4]*Matrix[2])*invdet;
  InvMatrix[7]= (Matrix[3]*Matrix[2]-Matrix[5]*Matrix[0])*invdet;
  InvMatrix[8]= (Matrix[4]*Matrix[0]-Matrix[3]*Matrix[1])*invdet;
//       for (int kr=0; kr<_nTdim; kr++) {
//          const  int kr1= (kr+1)%3; const  int kr2= (kr+2)%3;
//          for (int kc=0; kc<_nTdim; kc++) {
//            const  int kc1= (kc+1)%3; const  int kc2= (kc+2)%3;
//            InvMatrix[kc+kr*_nTdim] = (Matrix[3*kr1+kc1]*Matrix[3*kr2+kc2]-Matrix[3*kr1+kc2]*Matrix[3*kr2+kc1])*invdet;
//          }
//        }
  return det;
}
//===============================================================================================
template <> double MGFE<3,LAGRANGE_XFEM>::ComputeInverseMatrix_sur(double Matrix[], double InvMatrix[]) {
  double det = (Matrix[0] * Matrix[3] - Matrix[1] * Matrix[2]);
  double idet = 1. / det;
  InvMatrix[0] = Matrix[3] * idet;   // dxi dx
  InvMatrix[1] = -Matrix[2] * idet;  // deta dx
  InvMatrix[2] = -Matrix[1] * idet;  // dxi dy
  InvMatrix[3] = Matrix[0] * idet;   // deta dy
  return det;
}
//===============================================================================================
template <> double MGFE<3,LAGRANGE_XFEM>::ComputeInverseMatrix(double Matrix[], double InvMatrix[], int dim) {
  double det = 0.;

  if(dim == _dim)   det = ComputeInverseMatrix_vol(Matrix,InvMatrix);
  if(dim == _dim-1) det = ComputeInverseMatrix_sur(Matrix,InvMatrix);
  
  return det;
}

//===============================================================================================
template <> double MGFE<3,LAGRANGE_XFEM>::Tj_interp_lin(double *u,int j, double *xi,int ndof,int idim) {
  double sum=0.; /*int ndof=8;*/

  for(int k=0; k<ndof; k++) {
    double val=u[k+j];

    for(int kdim=0; kdim<idim; kdim++) {
      val *= (1.+ _CooH27[k+kdim*27]*xi[kdim]);
    }

    sum += val;
  }

  return sum/ndof;
}
//===============================================================================================
template <> double MGFE<3,LAGRANGE_XFEM>::dTjdxi_lin(double *u,int j, double *xi,int i,int ndof,int idim) {
  double sum=0.; /*int ndof=8;*/

  for(int k=0; k<ndof; k++) {
    double val=u[k+j];

    for(int kdim=0; kdim<i; kdim++) {
      val *= (1.+_CooH27[k+kdim*27]*xi[kdim]);
    }

    val *= _CooH27[k+i*27];

    for(int kdim=i+1; kdim<idim; kdim++) {
      val *= (1.+_CooH27[k+kdim*27]*xi[kdim]);
    }

    sum += val;
  }

  return sum/ndof;
}


// ===========================================================================================




template<> double  MGFE<3,LAGRANGE_XFEM>::JacOnGivenCanCoords(
const int dim,const double ElemCoords[], double CanCoords[], double InvJac[], int FamilyType, int nShape) {
  // double* LocDPhi = new double[dim];
  // double* CoordsDer = new double[dim * dim];
  // 
  // for(int dir = 0; dir < dim * dim; dir++) { CoordsDer[dir] = 0.; }
  // 
  // for(int s = 0; s < nShape; s++) {
  //   for(int dir = 0; dir < dim; dir++) {
  //     LocDPhi[dir] = FirstDerivateOfLocalPhi(s, CanCoords, dim, dir, FamilyType);
  //   }
  // 
  //   for(int dir1 = 0; dir1 < dim; dir1++)
  //     for(int dir2 = 0; dir2 < dim; dir2++) {
  //       CoordsDer[dir1 * dim + dir2] += ElemCoords[s + dir1 * nShape] * LocDPhi[dir2];
  //     }
  // }
  // 
  // double det = 0.;
  //   for(int dir = 0; dir < dim; dir++) {
  //     int sign = 1 - 2 * (dir % 2);
  //     int idx1 = (dir + 1) % dim;
  //     int idx2 = (dir + 2) % dim;
  //     det += sign * CoordsDer[dir * dim + dir]*(CoordsDer[idx1 * dim + idx1] * CoordsDer[idx2 * dim + idx2] -
  //            CoordsDer[idx1 * dim + idx2] * CoordsDer[idx2 * dim + idx1]);
  //   }
  // double idet = 1. / det;
  //   for(int row = 0; row < dim; row++) {   // LOOP OVER ROWS - row1<row2
  //     int sign_row = 1 - 2 * (row % 2);
  //     int row1 = ((row + 1) % dim < (row + 2) % dim) ? (row + 1) % dim : (row + 2) % dim;
  //     int row2 = ((row + 1) % dim > (row + 2) % dim) ? (row + 1) % dim : (row + 2) % dim;
  // 
  //     for(int col = 0; col < dim; col++) {   // LOOP OVER COLUMNS  - col1<col2
  //       int col1 = ((col + 1) % dim < (col + 2) % dim) ? (col + 1) % dim : (col + 2) % dim;
  //       int col2 = ((col + 1) % dim > (col + 2) % dim) ? (col + 1) % dim : (col + 2) % dim;
  //       int sign_col = 1 - 2 * (col % 2);
  //       InvJac[row * dim + col] = sign_col * sign_row * idet *
  //                                 (CoordsDer[row1 * dim + col1] * CoordsDer[row2 * dim + col2] -
  //                                  CoordsDer[row1 * dim + col2] * CoordsDer[row2 * dim + col1]);
  //     }
  //   }
  // 
  // delete[] LocDPhi;
  // delete[] CoordsDer;
  // return det;
  return 0.;
}


// ================================================
// JAC

// ================================================

//=============================================================
template <>
void MGFE<3,LAGRANGE_XFEM>::normal_g(
const double xx[],  // all surface coordinates <-
const double x_c[], ///< central point
double* normal_g   // normal ->
) const {          // ======================================


  
  int NDOF_FEMB=9;
  int PT2=3;
if(_GeomType==0){   NDOF_FEMB=6; PT2=2;}
  // coordinates
  double xx3D[3 * NDOF_FEMB];

  for(int i = 0; i < 3 * NDOF_FEMB; i++) { xx3D[i] = xx[i]; }
  double tg01[3];  // tangent  line
  double tg03[3];  // tangent plane

  //   the cross product of the two tangent vectors
  for(int i = 0; i < 3; i++) {
    tg01[i] = xx3D[1 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB];
    tg03[i] = xx3D[PT2 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB];

  }
  //   _mgutils.cross(tg01,tg03,normal_g);
  normal_g[0] = tg01[1] * tg03[2] - tg01[2] * tg03[1];
  normal_g[1] = tg01[2] * tg03[0] - tg01[0] * tg03[2];
  normal_g[2] = tg01[0] * tg03[1] - tg01[1] * tg03[0];
 if(normal_g[0] * (x_c[0] - xx3D[0]) + normal_g[1] * (x_c[1] - xx3D[0 + NDOF_FEMB])  + normal_g[2] * (x_c[2] - xx3D[0 +2* NDOF_FEMB]) > 0.) {
                 normal_g[0] *= -1.; tg01[0] *= -1.;   normal_g[1] *= -1.; tg01[1] *= -1.;   normal_g[2] *= -1.; tg01[2] *= -1.;
    //    std::cout << " Normal inverted ! ------------------------------------  \n";
  }
  // normalization -----------------------
  double mm = 0.;
  for(int idim = 0; idim < 3; idim++) { mm += normal_g[idim] * normal_g[idim]; }
  mm = std::sqrt(mm);
  for(int idim = 0; idim <3; idim++) { normal_g[idim] /= mm; }
  return;
}


//=============================================================
template <>
void MGFE<3,LAGRANGE_XFEM>::normal_g(
const double xx[],  // all surface coordinates <-
const double x_c[], ///< central point
double* normal_g,   // normal ->
int &sign
) const {          // ======================================


  
  int NDOF_FEMB=9;
  int PT2=3;
if(_GeomType==0){   NDOF_FEMB=6; PT2=2;}
  // coordinates
  double xx3D[3 * NDOF_FEMB];

  for(int i = 0; i < 3 * NDOF_FEMB; i++) { xx3D[i] = xx[i]; }
  double tg01[3];  // tangent  line
  double tg03[3];  // tangent plane

  //   the cross product of the two tangent vectors
  for(int i = 0; i < 3; i++) {
    tg01[i] = xx3D[1 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB];
    tg03[i] = xx3D[PT2 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB];

  }
  //   _mgutils.cross(tg01,tg03,normal_g);
  normal_g[0] = tg01[1] * tg03[2] - tg01[2] * tg03[1];
  normal_g[1] = tg01[2] * tg03[0] - tg01[0] * tg03[2];
  normal_g[2] = tg01[0] * tg03[1] - tg01[1] * tg03[0];
 if(normal_g[0] * (x_c[0] - xx3D[0]) + normal_g[1] * (x_c[1] - xx3D[0 + NDOF_FEMB])  + normal_g[2] * (x_c[2] - xx3D[0 +2* NDOF_FEMB]) > 0.) {
                 normal_g[0] *= -1.; tg01[0] *= -1.;   normal_g[1] *= -1.; tg01[1] *= -1.;   normal_g[2] *= -1.; tg01[2] *= -1.;
    //    std::cout << " Normal inverted ! ------------------------------------  \n";
  }
  // normalization -----------------------
  double mm = 0.;
  for(int idim = 0; idim < 3; idim++) { mm += normal_g[idim] * normal_g[idim]; }
  mm = std::sqrt(mm);
  for(int idim = 0; idim <3; idim++) { normal_g[idim] /= sign*mm; }
  
  
  
  
  
  return;
}


// ======================================
/// This function computes the normal at the gauss point
template <>
void MGFE<3,LAGRANGE_XFEM>::normal_g(
const double* xx,  // all surface coordinates <-
double* normal_g   // normal ->
) const {          // ======================================
  int NDOF_FEMB=9;
  int PT2=3;
if(_GeomType==0){   NDOF_FEMB=6; PT2=2;}
  // coordinates
  double xx3D[3 * NDOF_FEMB];

  for(int i = 0; i < 3 * NDOF_FEMB; i++) { xx3D[i] = xx[i]; }
  double tg01[3];  // tangent  line
  double tg03[3];  // tangent plane

  //   the cross product of the two tangent vectors
  for(int i = 0; i < 3; i++) {
    tg01[i] = xx3D[1 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB];
    tg03[i] = xx3D[PT2 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB];

  }
  //   _mgutils.cross(tg01,tg03,normal_g);
  normal_g[0] = tg01[1] * tg03[2] - tg01[2] * tg03[1];
  normal_g[1] = tg01[2] * tg03[0] - tg01[0] * tg03[2];
  normal_g[2] = tg01[0] * tg03[1] - tg01[1] * tg03[0];

  // normalization -----------------------
  double mm = 0.;
  for(int idim = 0; idim < 3; idim++) { mm += normal_g[idim] * normal_g[idim]; }
  mm = std::sqrt(mm);
  for(int idim = 0; idim <3; idim++) { normal_g[idim] /= mm; }
  
  

  return;
}
// ===================================================================
template <>
void  MGFE<3,LAGRANGE_XFEM>::Oxy_face(
const double* xx,   ///< coordinates
const double x_c[], ///< central point
double* normal_g,   ///<  i_n,i_t1,i_t2 face Oxyz
double* tg01,
double* /*tg03*/
) const {            // ======================================
  
  int NDOF_FEMB=9; int PT2=3; // quad9 geometry
  if(_GeomType==0){   NDOF_FEMB=6; PT2=2;}  // tet10 geometry
  double xx3D[3 * NDOF_FEMB]; // coordinates
  for(int i = 0; i < _dim * NDOF_FEMB; i++) { xx3D[i] = xx[i]; }

  double tg03[3];  // tangent plane
  //   the cross product of the two tangent vectors
  for(int i = 0; i < 3; i++) {
    tg01[i] = xx3D[1 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB];
    tg03[i] = xx3D[PT2 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB];
  }
  //   _mgutils.cross(tg01,tg03,normal_g);
  normal_g[0] = tg01[1] * tg03[2] - tg01[2] * tg03[1];
  normal_g[1] = tg01[2] * tg03[0] - tg01[0] * tg03[2];
  normal_g[2] = tg01[0] * tg03[1] - tg01[1] * tg03[0];
  // if the sign is not correct then reverse it
  if(normal_g[0] * (x_c[0] - xx3D[0])+ normal_g[1] * (x_c[1] - xx3D[0 + NDOF_FEMB]) +
      normal_g[2] * (x_c[2] - xx3D[0 + 2 * NDOF_FEMB]) >  0) {
     for(int i = 0; i < 3; i++) {  normal_g[i] *= -1;     tg01[i] *= -1;}
  }
  // normalization -----------------------
  double mn = 0.; double mtg = 0.;
  for(int idim = 0; idim < _dim; idim++) {mn += normal_g[idim]*normal_g[idim]; mtg +=tg01[idim]*tg01[idim];}
  mn = sqrt(mn); mtg=sqrt(mtg); for(int idim=0;idim <3;idim++) {normal_g[idim] /= mn;tg01[idim] /=mtg; }

  return;
}




//===================================================================
//===================================================================
// GET  shape dshape ddshape
//===================================================================
//===================================================================


//===================================================================
template <> void MGFE<3,LAGRANGE_XFEM>::get_dphi_on_given_node(
  const int dim, 
  double ElemCoords[], 
  double CanPos[], 
  double dphi[]
) {//===================================================================
  // double* InvJac = new double[dim * dim];
  // double* gradphi_g = new double[dim];
  // int nShape = _NoShape;
  // JacOnGivenCanCoords(dim, ElemCoords, CanPos, InvJac, _GeomType, nShape);
  // 
  // for(int eln = 0; eln < nShape; eln++) {
  //   for(int idim = 0; idim < dim; idim++) {
  //     gradphi_g[idim] = FirstDerivateOfLocalPhi(eln, CanPos, dim, idim, _GeomType);
  //   }
  // 
  //   for(int idim = 0; idim < dim; idim++) {
  //     double sum = 0.;
  // 
  //     for(int jdim = 0; jdim < dim; jdim++) { sum += InvJac[jdim + idim * dim] * gradphi_g[jdim]; }
  // 
  //     dphi[eln + idim * nShape] = sum;
  //   }
  // }
  // 
  // delete[] gradphi_g;
  // delete[] InvJac;
  return;
}



// ==================================================================
/// This function computes the shape values at the gauss point qp
template <> void MGFE<3,LAGRANGE_XFEM>::get_phi_gl_g(
const int qp,    // gaussian point <-
double phi[]     // shape functions ->
) {                  // =================================================
  for(int ish = 0; ish < _NoShape; ish++) phi[ish] = phi[ish * _NoGauss + qp];
  return;
}

// ==================================================================
/// This function computes the shape values at the gauss point qp
template <> void MGFE<3,LAGRANGE_XFEM>::get_phi_gl_g(
const int qp,             // gaussian point   <-
std::vector<double>& phi  // shape functions  ->
) {                           // =====================================
  for(int ish = 0; ish < _NoShape; ish++) phi[ish] = phi[ish * _NoGauss + qp];

  return;
}

// =================================================================
/// Shape functions derivatives dphi[id+i* el_nnodes] =Ti
/// Tx Ty Tz
///  tensor  order  at the gauss point qp, node id:
///  we have  dphi[Tx(0),Tx(1),..Tx(id), Ty(0),Ty(1),..Ty(id), Tz(0),Tz(1),..Tz(id) ]
template <> void MGFE<3,LAGRANGE_XFEM>::get_dphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double dphi[]           // global derivatives ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  double gradphi_g[3];  // temp grad phi

  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + qp;
    for(int idim = 0; idim < _dim; idim++) gradphi_g[idim]=_dphidxez_map[lqp+idim*goffset];
    for(int idim = 0; idim < _dim; idim++) {
      double sum=0.; for(int jdim=0;jdim<_dim;jdim++) sum +=InvJac[jdim+idim*_dim]*gradphi_g[jdim];
      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

// =============================================================
/// Shape functions 2nd derivatives at the gauss point qp
/// ddphi[id*_dim*_dim+i*_dim+j]= Tij
///    Txx Txy Txz
///    Tyx Tyy Tyz
///    Tzx Tzy Tzz
///
///     ddphi[Txx(0),Txy(0),Txz(0),Tyx(0),Tyy(0),Tyz(0),Tzx(0),Tzy(0),Tzz(0),
///                 ....................................
///                 Txx(id),Txy(id),Txz,Tyx(id),Tyy,Tyz,Tzx(id),Tzy(id),Tzz(id)]
///   tensor  order  at the gauss point qp, node inode
template <> void MGFE<3,LAGRANGE_XFEM>::get_ddphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double ddphi[]          // global derivatives ->
) {                         // =========================================
int NDOF_FEM=27;
if(_GeomType==0)  NDOF_FEM=6;
  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset
  double ddphi_loc[NDOF_FEM * 3 * 3];

  for(int eln = 0; eln < el_nnodes; eln++) {   // LOOP OVER NODES ===================================
    const int shift = eln * _dim * _dim;       // offset for node derivatives
    const int lqp = eln * el_ngauss + qp;

    for(int idim = 0; idim < _dim * _dim; idim++) {
      ddphi_loc[shift + idim] = _dphidxx_map[lqp + idim * goffset];
    }

    double Hess_tmp[3 * 3], JacTf[3 * 3];

    for(int init = 0; init < _dim; init++)
      for(int jnit = 0; jnit < _dim; jnit++) {
        Hess_tmp[init * _dim + jnit] = 0.;
        ddphi[shift + init * _dim + jnit] = 0.;
        JacTf[init * _dim + jnit] = InvJac[jnit * _dim + init];
      }

    for(int ii = 0; ii < _dim; ii++)
      for(int jj = 0; jj < _dim; jj++)
        for(int ss = 0; ss < _dim; ss++) {
          Hess_tmp[ii * _dim + jj] += ddphi_loc[shift + ii * _dim + ss] * JacTf[ss * _dim + jj];
        }

    // calculting d2 phi / dxi dxj
    for(int ii = 0; ii < _dim; ii++)
      for(int jj = 0; jj < _dim; jj++)
        for(int ss = 0; ss < _dim; ss++) {
          ddphi[shift + ii * _dim + jj] += InvJac[ii * _dim + ss] * Hess_tmp[ss * _dim + jj];
        }
  }  // END LOOP OVER NODES =====================================================================

  // Pay attention, different order with respect to first order derivatives!
  //   delete[]ddphi_loc;
  return;
}
// =============================================================
/// Shape functions derivatives at the nodal points
template <> void MGFE<3,LAGRANGE_XFEM>::get_dphi_node(
const int /*kdim*/,         // dimension <-
const int node,         // nodal point <-
const double InvJac[],  // Jacobean <-
double dphi[]           // global derivatives at the nodal point ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset
  double gradphi_g[3];  // temp grad phi
  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + node;
    for(int idim = 0; idim < _dim; idim++) { gradphi_g[idim] = _dphidxez_map_nodes[lqp + idim * goffset]; }
    for(int idim = 0; idim < _dim; idim++) {
      double sum = 0.;   for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * gradphi_g[jdim]; }
      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}



// =============================================================
/// Shape functions derivatives at the gauss point qp
template <> void MGFE<3,LAGRANGE_XFEM>::get_dphi_gl_g(
const int /*kdim*/,         // dimension <-
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double dphi[],          // global derivatives ->
int sdim) {             // =========================================

  const int el_nnodes = _NoShape;   // # of shape functions
  const int el_ngauss = _NoGauss;  // # of guass points
  const int goffset = el_nnodes * el_ngauss;  // gauss offset
  double gradphi_g[3];  // temp grad phi
  for(int idim = 0; idim < 3; idim++) { gradphi_g[idim] = 0.; }
  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + qp;
    for(int idim = 0; idim < sdim; idim++) { gradphi_g[idim] = _dphidxez_map[lqp + idim * goffset]; }
    for(int idim = 0; idim < sdim; idim++) {
      double sum = 0.; for(int jdim = 0; jdim < sdim; jdim++) { sum += InvJac[jdim + idim * sdim] * gradphi_g[jdim]; }
      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

// =============================================================
/// Shape functions derivatives at the gauss point qp
template <> void  MGFE<3,LAGRANGE_XFEM>::get_dphi_gl_g(
const int /*kdim*/,            // dimension <-
const int qp,              // gaussian point <-
const double InvJac[],     // Jacobean
std::vector<double>& dphi  // global derivatives ->
) {                            // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of guass points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset
  double dphidxi_g[3];  // temp grad phi
  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + qp;
    for(int idim = 0; idim < 3; idim++) { dphidxi_g[idim] = _dphidxez_map[lqp + idim * goffset]; }
    for(int idim = 0; idim < 3; idim++) {
      double sum = 0.;  for(int jdim = 0; jdim <3; jdim++) { sum += InvJac[jdim + idim * _dim] * dphidxi_g[jdim]; }
      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}



//  ========================================================================
// SET parameter
// =========================================================================

// ===================================================================================
template <> void MGFE<3,LAGRANGE_XFEM>::set_par(
  std::string fileh5,
  int dim_in,
  int order_in,
  int shapes_quad,
  int geom_type_in,
  int fam_type_in
) {_dim=dim_in, _order= order_in,  _NoShape=shapes_quad,   _GeomType=geom_type_in,     _FamType=fam_type_in;
   femfileh5_=fileh5;
     ///< family type  0=Lagrangian(C^0) 1=Raviart Thomas (L^2-C^0)
  // Fem order
  // =============================
    //  ------------------------------------------------------------------------
  // gauss points (_NoGauss) ------------------------------------------------------
  _NoGauss=0;  assert(_qrule); _NoGauss=_qrule->n_points();
  if(_NoGauss==0) {   std::cout << "Quadrature not assigned"; abort();}

  // number of nodes (_NoElNodes) in the reference quad element
    // number of nodes (_NoElNodes) in the reference quad element
   ElemType  elem_name=_fe_type.fetype; 
   ElemType  geom_elem_name=get_ref_elemtype(elem_name); 
   Order elem_order= _fe_type.order ;
    _order= (int)elem_order;
  // number of nodes (_NoElNodes) in the reference EDGE3 element
    _NoElNodes= get_refspace_n_nodes( geom_elem_name);// mesh nodes
  _NoShape=get_refspace_n_nodes(elem_name); // element nodes 
   
  
   
   
  // //  ------------------------------------------------------------------------
  // // gauss points (_NoGauss) ------------------------------------------------------
  //  _NoGauss=27; 
  // 
  // // number of nodes (_NoElNodes) in the reference quad element
  // _NoElNodes=_NoShape;  

  //  -----------------------------------------------------------------------------------
  // x node  coordinates   (_xnodes[xp1,xp2,xp3,...yp1,yp2,yp3,....zp1,zp2,zp3,...] )----
  _xnodes=new double[ _NoElNodes*_dim];
  // shape functions at nodes (_phi_map_nodes[sh1(x1),sh1(x2),....  sh2(x1),sh2(x2),.... ])
  _phi_map_nodes = new double[_NoShape* _NoElNodes];
  // 1st derivatives at nodes x _dphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxez_map_nodes = new double[_dim * _NoShape * _NoElNodes];
  // element dofs  (reference to nodes) --------------------------------------------------
  _fem_conn=new int[_NoShape];
  // gaussian (g) coordinates (_xg[xg1,xg2,xg3,...yg1,yg2,yg3,....zg1,zg2,zg3,...] )-----
  _xg=new double[_NoGauss*_dim];
  // weights [wg1,wg2,wg3,.......]
  _weight = new double[_NoGauss];
  // shape functions phi[sh1(xg1),sh1(xg2),....  sh2(xg1),sh2(xg2),.... ]
  _phi_map = new double[ _NoGauss*_NoShape ];
  // 1st derivatives at g _dphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxez_map = new double[_dim* _NoShape* _NoGauss];
  // 2nd derivatives at g _ddphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxx_map = new double[_dim * _dim  * _NoShape * _NoGauss];
  //   -------------------------------------------------------------------------------------
  
  
   // _deg=0;// DDDDDDDDDDDDDDDDDDDDDDDDD
  
  
 
  
//     TIME(std::clock_t start_time = std::clock();)     //  TC +++++++++++++++
// dofs  ---------------------------------------------------------------------------
// dofs -----------------------------------------------------------------------
 for(int j = 0; j < _NoShape; j++) _fem_conn[j]=j;
 double *xnodes=new double [_NoElNodes*_dim];
 double *point=new double[_dim];
    //
    // ---------------------------------------------------------------------------------------
    // nodes ----------------------------------------------------------------------------
  get_refspace_nodes(geom_elem_name,xnodes);// nodes[(x,y)_0,(x,y)_1.... ]
    // _xnodes[0]=0.;    _xnodes[1]=1.;    _xnodes[2]=0.5;
    for(int i = 0; i < _NoElNodes; i++) {  // points
      for(int dir = 0; dir < _dim; dir++) { 
        _xnodes[i+_NoElNodes*dir]=xnodes[dir+_dim*i]; 
        point[dir] = _xnodes[i+_NoElNodes*dir];
      }
      
    for(int j = 0; j < _NoShape; j++) {  // shapse
      //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
      _phi_map_nodes[j*_NoElNodes+i]   = fe_lagrange_3D_shape<LAGRANGE_XFEM>(elem_name , elem_order, j, point);// Rec_Quad_Phi(j, pointNode,_dim);

      // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
      for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
        _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(elem_name, elem_order, j,dir, point);// Rec_Quad_DPhi(j,pointNode,_dim,dir);     // dphi/dxi
      }
    }
  }

  // delete[]pointNode;
 delete[]xnodes;


  for(int i = 0; i < _NoGauss; i++) {
    _weight[i] =_qrule->w(i); // weight
     _qrule->qp(point,i);     // ith-gaussian point ->point(i)
  
    for(int j = 0; j < _NoShape; j++) {     // loop over test function id
      // phi[phi_1(xg1),phi_1(xg2),phi_1(xg3),...phi_2(xg1),phi_2(xg2),....]
      phi[j*_NoGauss+i] = fe_lagrange_3D_shape<LAGRANGE_XFEM>(elem_name , elem_order, j, point);//Rec_Quad_Phi(j, point, _dim); // shape functions

      for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
        _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] = fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(elem_name, elem_order, j,dir, point);//Rec_Quad_DPhi(j, point, _dim, dir);

        for(int dir2 = 0; dir2 < _dim; dir2++)
          _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =fe_lagrange_3D_shape_second_deriv<LAGRANGE_XFEM>(elem_name, elem_order,j,dir +dir2, point);//       Rec_Quad_D2Phi(j, point, _dim, dir, dir2);  // d2/dxdy gaussian points
      }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
  } //--------------------------------------------------------------------------------
  delete[]point;  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
//   if(_FamType==0) {  // 0= Lagrangian  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
//     // lagrangian -> n_point= n_shape =_order DDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
//     //            -> pol degree=_order/2 (int div)
//     _deg=0;// DDDDDDDDDDDDDDDDDDDDDDDDD
// 
//     switch(_order) {  // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
//     case 0:
//       init_pie();  break;
//     case 1:
//       init_lin();  break;
//     case 2:
//       init_qua();  break;
//     default:
//       std::cout<<"MGFE:"<<" Dim "<<_dim<<" case "<<_order<<" Not implemented";  exit(3);
//     }// OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// 
//   }//  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
// 
//   if(_FamType==1) {  //1= Raviart FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
//     if(_GeomType==1) {  //HEX
//       _deg=1;   if(_dim==1)  _deg=0;
//       switch(_order) {
//       case 0:
//         init_pie();   break;
//       case 1:
//         init_lin_rec_rt();   break;
//       default:
//         std::cout<<"MGFE::MGFE:"<<" Dim "<<_dim<<" case "<<_order<<"RT Not implemented"; exit(3);
//       }
//     }
//     if(_GeomType==0) {  // TETRA GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
//       std::cout<< " Not implemented";
//       exit(0);
//     }//  GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
//   }//  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
//   
//   if(_FamType==2) {  //2= discontinuous FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
//     std::cout<< " Not implemented";
//     exit(0);
//   } //  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

#ifdef PRINT_INFO
  std::cout << "\n MGFE: Dimension: " << _dim ;
  std::cout << "\n MGFE: gaussian points: " << _NoGauss ;
  std::cout << " MGFE: number of shape functions " << _NoShape ;
  std::cout << " MGFE: polynomial order " << _order ;
  std::cout << " \n";
#endif
  return;
}

template <>  bool MGFE<3,LAGRANGE_XFEM>::phys_to_ref(double xref[],double xphys[],const double Xel[], int max_it, double tol) {
  
  
  std::array<double,3> xi = {0.0, 0.0, 0.0};

    auto norm2 = [](const std::array<double,3>& v){
        return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    };

    for (int it = 0; it < max_it; ++it) {
        std::array<double,3> xmap = {0.0, 0.0, 0.0};
        double J[3][3] = {
            {0.0,0.0,0.0},
            {0.0,0.0,0.0},
            {0.0,0.0,0.0}
        };

        double p[3] = {xi[0], xi[1], xi[2]};

        for (int i = 0; i < _NoShape; ++i) {
            // === Shape functions del brick ===
            // const double Ni     = Rec_Quad_Phi (i, p, /*_dim=*/3);
            // const double dN_dxi = Rec_Quad_DPhi(i, p, /*_dim=*/3, /*dir=*/0);
            // const double dN_deta= Rec_Quad_DPhi(i, p, /*_dim=*/3, /*dir=*/1);
            // const double dN_dze = Rec_Quad_DPhi(i, p, /*_dim=*/3, /*dir=*/2);
            const double Ni     = fe_lagrange_3D_shape<LAGRANGE_XFEM>(HEX27,SECOND,i, p);
            const double dN_dxi = fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(HEX27,SECOND, i,0, p);
            const double dN_deta= fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(HEX27,SECOND, i,1, p);
            const double dN_dze = fe_lagrange_3D_shape_deriv<LAGRANGE_XFEM>(HEX27,SECOND, i,2, p);
            
        

            
          

            const double Xi = Xel[i + _NoShape*0];
            const double Yi = Xel[i + _NoShape*1];
            const double Zi = Xel[i + _NoShape*2];

            // mappa x(xi)
            xmap[0] += Ni * Xi;
            xmap[1] += Ni * Yi;
            xmap[2] += Ni * Zi;

            // Jacobiano J = dx/d(ref)
            J[0][0] += dN_dxi  * Xi;   // dx/dxi
            J[0][1] += dN_deta * Xi;   // dx/deta
            J[0][2] += dN_dze  * Xi;   // dx/dzeta

            J[1][0] += dN_dxi  * Yi;   // dy/dxi
            J[1][1] += dN_deta * Yi;   // dy/deta
            J[1][2] += dN_dze  * Yi;   // dy/dzeta

            J[2][0] += dN_dxi  * Zi;   // dz/dxi
            J[2][1] += dN_deta * Zi;   // dz/deta
            J[2][2] += dN_dze  * Zi;   // dz/dzeta
        }

        // residuo r = x(xi) - xphys
        std::array<double,3> r = {
            xmap[0] - xphys[0],
            xmap[1] - xphys[1],
            xmap[2] - xphys[2]
        };

        const double r2 = norm2(r);
        const double scale = 1.0 + xphys[0]*xphys[0] + xphys[1]*xphys[1] + xphys[2]*xphys[2];
        if (r2 <= tol*tol*scale) return {true};

        // Risolvi J * d = r (Newton step: xi_{k+1} = xi_k - d)
        const double a = J[0][0], b = J[0][1], c = J[0][2];
        const double d = J[1][0], e = J[1][1], f = J[1][2];
        const double g = J[2][0], h = J[2][1], k = J[2][2];

        const double det =
              a*(e*k - f*h)
            - b*(d*k - f*g)
            + c*(d*h - e*g);

        if (std::abs(det) < 1e-30) {
            return {false};
        }

        // inv(J) = adj(J)/det
        const double inv00 = (e*k - f*h) / det;
        const double inv01 = (c*h - b*k) / det;
        const double inv02 = (b*f - c*e) / det;

        const double inv10 = (f*g - d*k) / det;
        const double inv11 = (a*k - c*g) / det;
        const double inv12 = (c*d - a*f) / det;

        const double inv20 = (d*h - e*g) / det;
        const double inv21 = (b*g - a*h) / det;
        const double inv22 = (a*e - b*d) / det;

        std::array<double,3> delta;
        delta[0] = inv00*r[0] + inv01*r[1] + inv02*r[2];
        delta[1] = inv10*r[0] + inv11*r[1] + inv12*r[2];
        delta[2] = inv20*r[0] + inv21*r[1] + inv22*r[2];

        // line search (backtracking)
        double alpha = 1.0;
        double xi_try[3];  xi_try[0] = xref[0];xi_try[1] = xref[1];xi_try[2] = xref[2];

        for (int ls = 0; ls < 10; ++ls) {
            for(int kdim=0; kdim<3;kdim++) xi_try[kdim] = xref[kdim] - alpha*delta[kdim];

            double p2[3] = {xi_try[0], xi_try[1], xi_try[2]};
            std::array<double,3> x2 = {0.0, 0.0, 0.0};

            for (int i = 0; i < _NoShape; ++i) {
                // const double Ni = Rec_Quad_Phi(i, p2, /*_dim=*/3);
                const double Ni =fe_lagrange_3D_shape<LAGRANGE_XFEM>(HEX27,SECOND,i, p2);
                x2[0] += Ni * Xel[i + _NoShape*0];
                x2[1] += Ni * Xel[i + _NoShape*1];
                x2[2] += Ni * Xel[i + _NoShape*2];
            }

            std::array<double,3> r_try = {
                x2[0] - xphys[0],
                x2[1] - xphys[1],
                x2[2] - xphys[2]
            };

            if (norm2(r_try) < r2) break;
            alpha *= 0.5;
        }

        xref = xi_try;
    }

    return {false};
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
