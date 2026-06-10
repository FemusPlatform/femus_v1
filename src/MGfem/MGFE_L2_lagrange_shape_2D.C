// // // The libMesh Finite Element Library.
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
// std include ----------------->
#include <cmath>
#include <array>
#include <iomanip>
// Femus configuration include ----------------->
#include "Printinfo_conf.h"
// Femus includes ------------------------>
#include "MGFE_L2_FE_template.h"
#include "MGFE_L0_C.h"
#include "MGFE_L2_lagrange_shape_1D.h"
#include "MGquadrature_L1_gauss.h"
#include "MGenum_quadrature_type.h"

// ================================================================================================
// ================================================================================================
// // Auxiliary (Anonymous namespace) functions shared by LAGRANGE and
// // L2_LAGRANGE implementations. Implementations appear at the bottom
// // of this file.
// ==============================================================================================
// ================================================================================================
template <FEFamily T> double fe_lagrange_2D_shape
(const ElemType type_of_elem,const Order order, const unsigned int phi_i, const double  point[]); 
template <FEFamily T> double fe_lagrange_2D_shape_deriv
(const ElemType type_of_elem, const Order order, const unsigned int i, const unsigned int j, const double  p[]);
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <FEFamily T> double fe_lagrange_2D_shape_second_deriv
 (const ElemType type_of_elem,const Order order,const unsigned int phi_i,const unsigned int dxj,const double point[]);
// #endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
// ================================================================================================
// ================================================================================================

// =================================================================================================
// Constructor-Destructor  
// =================================================================================================
// // FE class members
template <> MGFE<2,LAGRANGE>::MGFE (const FEType & fet) :   MGFEbase (2,fet),
   // _dim(2),             ///< space dimension
   //  _order(2),         ///< pol approx (2=quad  1=lin 0=k)
   //  _NoShape(9),      ///< number of shapes
   //  _GeomType(1),  ///< geom type 0=hex/quad  1=tetra/tri
   //  _FamType(0), 
   last_side(INVALID_ELEM),
   last_edge(INVALID_ELEM)
 {// ==============================================================================================
   // Sanity check.  Make sure the Family specified in the template instantiation matches the one in the FEType object
   assert (LAGRANGE == this->get_family());
 }
  
  
// static std::unique_ptr<MGFE<2,LAGRANGE>> build (const unsigned int dim, const FEType & type){}
// ==============================================================================================
template <> MGFE<2,LAGRANGE>::~MGFE(){}




 template <>  const int MGFE<2,LAGRANGE>::get_dim()       {return _dim;}  ///< Dimension (3 or 2) used in the 3D FEM (1< _dim < MDIM_FEM)
 template <>  const int MGFE<2,LAGRANGE>::get_order()     {return _order;}///< Order of the shape functions 0=cpnst 1=linear 2=quad
 template <>  const int MGFE<2,LAGRANGE>::get_NoShape()   {return _NoShape;}  ///< Number of shape functions of the fem
 template <>  const int MGFE<2,LAGRANGE>::get_NoElNodes() {return _NoElNodes;}  ///< Number of nodes in the reference element (always ref quad)
 template <>  const int MGFE<2,LAGRANGE>::get_GeomType()  {return _GeomType;}   ///< Geom type 1=HEXAHEDRAL 0=TETRAHEDRAL
 template <>  const int MGFE<2,LAGRANGE>::get_FamType()   {return _FamType;}    ///< Family  (geom 1: Lagrangian =0  RT =1)
  
   // Gaussian points -----------------------------------------------------------------------------
 template <>  const int MGFE<2,LAGRANGE>::get_NoGauss(){ return _NoGauss;}  ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
  template <> const double  MGFE<2,LAGRANGE>::get_xg(int i){          return  _xg[i];}          ///< gaussian points
  template <> const double MGFE<2,LAGRANGE>::get_weight(int i){      return _weight[i];}            ///< Weight
  template <> const double & MGFE<2,LAGRANGE>::get_phi(){     return  *phi;}           ///< Shape functions
  template <> const double & MGFE<2,LAGRANGE>::get_dphidxez_map(){return  *_dphidxez_map;}      ///< Shape derivative functions in gaussian points
  template <> const double & MGFE<2,LAGRANGE>::get_dphidxx_map(){ return  *_dphidxx_map;}       ///< Second order shape derivatives in Gaussian points
  
   // Dofs points ----------------------------------------------------------------------------------
 template <>  const double &  MGFE<2,LAGRANGE>::get_xnodes(){return   *_xnodes;}              ///< xnode (hex or tetra) in (xi,eta,seta)
 template <>  const int &  MGFE<2,LAGRANGE>::get_fem_conn(){ return *_fem_conn;}              ///< dof node coords (respect to _xnodes)
  // // template <> const int     MGFE<2,LAGRANGE>::get_deg(){return _deg;}                  ///< Polinomial degree
  
   // Nodes (NDOF_FEM) -----------------------------------------------------------------------------
 template <>  const double & MGFE<2,LAGRANGE>::get_phi_map_nodes(){return   *_phi_map_nodes;   }          ///< Shape functions
 template <>  const double & MGFE<2,LAGRANGE>::get_dphidxez_map_nodes(){return   *_dphidxez_map_nodes;} 
 
// ===========================================================================
 template <> void MGFE<2,LAGRANGE>::attach_quadrature_rule (
   QBase *q
 ){ // ===========================================================================
   assert(q);
   this->_qrule = q;
   // make sure we don't cache results from a previous quadrature rule
   this->_elem_type = INVALID_ELEM;
   return;
 }
// =========================================================================== 
 template <> unsigned int MGFE<2,LAGRANGE>::n_shape_functions (
) const {// ===========================================================================
    return MGFE<2,LAGRANGE>::n_dofs (this->_elem_type, this->_fe_type.order + this->_p_level);
  }
  
// ===========================================================================
template <> int MGFE<2,LAGRANGE>::GetFamilyType(
  int elem_dof, 
  int dim
) {//=========================================================================
    int FamilyType=1; //lagrange27
    if(elem_dof == 3 || elem_dof == 6) {  FamilyType = 0;   }
    else if(elem_dof == 4 || elem_dof == 9) {   FamilyType = 1;  }
    else  printf( "MGFE1::GetFamilyType Unknown FamilyType for element with dimension %d and %d number of dof \n",
      dim, elem_dof);

  return FamilyType;
} 
 

//=================================================================================================
/// This function writes shape and derivative values at the gaussian points
 template <> void MGFE<2,LAGRANGE>::write_c(
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
// ===========================================================================   
 template <> void MGFE<2,LAGRANGE>::write(
const std::string& name  // file <-
) {  // ======================================================================
  std::ofstream in(name.c_str());
  this->write_c(in);
}  
   
// ===========================================================================
template <> void MGFE<2,LAGRANGE>::print_rec_xdmf(
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


// =========================================================================================
//  SHAPE   FUNCIONS,DERIVS,SECOND_DERIVS
// =========================================================================================


// =========================================================
// This function computes the shape functions  phi(i) (ElemType,Order)
// at the point p[] (xi,eta)
template <> double MGFE<2,LAGRANGE>::shape(
  const ElemType     type,  //  type of element (HEX# TRI# ...)
  const Order        order, //  Polynomial Order (FIRST, SECOND, ...) 
  const unsigned int i,     //  phi(i)
  const double       p[]    //  point[xi,eta] 
){ // ======================================================
  return fe_lagrange_2D_shape<LAGRANGE>(type, order, i, p);
}
//   LIBMESH ********************************************************
// template <> double MGFE<2,L2_LAGRANGE>::shape(const ElemType type, const Order order,
//                               const unsigned int i,  const double & p)
// {  return fe_lagrange_2D_shape<L2_LAGRANGE>(type, order, i, p); }
// -------------
// template <> double MGFE<2,LAGRANGE>::shape(const Elem * elem,   const Order order,  const unsigned int i,
//                            const Point & p,   const bool add_p_level) {
//   libmesh_assert(elem);
// call the orientation-independent shape functions
//   return fe_lagrange_2D_shape<LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), i, p); }
// -----------
// template <> double MGFE<2,L2_LAGRANGE>::shape(const Elem * elem,   const Order order,   const unsigned int i,
//                               const Point & p,    const bool add_p_level)
// {  libmesh_assert(elem);   // call the orientation-independent shape functions
//   return fe_lagrange_2D_shape<L2_LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), i, p);}
// ------------
// template <> double MGFE<2,LAGRANGE>::shape(const FEType fet,  const Elem * elem,    const unsigned int i,
//                            const Point & p, const bool add_p_level)
// {  libmesh_assert(elem);  return fe_lagrange_2D_shape<LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), i, p); }
// -----------
// template <> double MGFE<2,L2_LAGRANGE>::shape(const FEType fet,  const Elem * elem,  const unsigned int i,
//                               const Point & p,  const bool add_p_level) {
//   libmesh_assert(elem);  return fe_lagrange_2D_shape<L2_LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), i, p);}
// 


// =======================================================================
template <> double MGFE<2,LAGRANGE>::shape_deriv(
  const ElemType type,  
  const Order order,
  const unsigned int i,  
  const unsigned int j, 
  const double  p[]
) {// ====================================================================
  return fe_lagrange_2D_shape_deriv<LAGRANGE>(type, order, i, j, p); 
}
//   LIBMESH ********************************************************
// ----------------
// template <>double MGFE<2,L2_LAGRANGE>::shape_deriv(const ElemType type,   const Order order,   const unsigned int i,     const unsigned int j,
//                                     const Point & p)
// {return fe_lagrange_2D_shape_deriv<L2_LAGRANGE>(type, order, i, j, p); }
// -----------------------
// template <>
// double MGFE<2,LAGRANGE>::shape_deriv(const Elem * elem,   const Order order,        const unsigned int i,
//                                  const unsigned int j,   const Point & p,    const bool add_p_level)
// {  libmesh_assert(elem);  // call the orientation-independent shape functions
//   return fe_lagrange_2D_shape_deriv<LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), i, j, p);}
// ----------------
// template <> double MGFE<2,L2_LAGRANGE>::shape_deriv(const Elem * elem,  const Order order,
//                                     const unsigned int i, const unsigned int j,   const Point & p,  const bool add_p_level)
// {  libmesh_assert(elem);  // call the orientation-independent shape functions
//   return fe_lagrange_2D_shape_deriv<L2_LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), i, j, p);}
//  ---------------
// template <> double MGFE<2,LAGRANGE>::shape_deriv(const FEType fet,   const Elem * elem,    const unsigned int i,
//                                  const unsigned int j, const Point & p,    const bool add_p_level)
// {   libmesh_assert(elem);   return fe_lagrange_2D_shape_deriv<LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), i, j, p); }
// -----------------
// template <>  double MGFE<2,L2_LAGRANGE>::shape_deriv(const FEType fet,   const Elem * elem,   const unsigned int i,
//                                     const unsigned int j,    const Point & p,                   const bool add_p_level)
// {   libmesh_assert(elem);  return fe_lagrange_2D_shape_deriv<L2_LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), i, j, p); }
// 



// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// ================================================================
template <> double MGFE<2,LAGRANGE>::shape_second_deriv(
  const ElemType type,
  const Order order,
  const unsigned int i,
  const unsigned int j,
  const double  p[]
) {// ======================================================
   return fe_lagrange_2D_shape_second_deriv<LAGRANGE>(type, order, i, j, p); 
}    
//   LIBMESH ********************************************************  
// ---------------- 
// template <> double MGFE<2,L2_LAGRANGE>::shape_second_deriv(const ElemType type,  const Order order,
//                                            const unsigned int i,const unsigned int j,  const Point & p)
// {  return fe_lagrange_2D_shape_second_deriv<L2_LAGRANGE>(type, order, i, j, p); }
// ----------------
// template <> double MGFE<2,LAGRANGE>::shape_second_deriv(const Elem * elem,  const Order order,
//                                         const unsigned int i,   const unsigned int j,   const Point & p,  const bool add_p_level)
// {   libmesh_assert(elem);   // call the orientation-independent shape functions
//   return fe_lagrange_2D_shape_second_deriv<LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), i, j, p);}
// --------------------
// template <> double MGFE<2,L2_LAGRANGE>::shape_second_deriv(const Elem * elem, const Order order, const unsigned int i,
//  const unsigned int j, const Point & p, const bool add_p_level) { libmesh_assert(elem);
//   // call the orientation-independent shape functions
//   return fe_lagrange_2D_shape_second_deriv<L2_LAGRANGE>(elem->type(), order + add_p_level*elem->p_level(), i, j, p);
// }
// template <> double MGFE<2,LAGRANGE>::shape_second_deriv(const FEType fet, const Elem * elem, const unsigned int i,
//                                         const unsigned int j, const Point & p, const bool add_p_level){
//   libmesh_assert(elem); return fe_lagrange_2D_shape_second_deriv<LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), i, j, p);}
// template <> double MGFE<2,L2_LAGRANGE>::shape_second_deriv(const FEType fet,
//  const Elem * elem, const unsigned int i, const unsigned int j,  const Point & p,   const bool add_p_level) {
//   libmesh_assert(elem);  return fe_lagrange_2D_shape_second_deriv<L2_LAGRANGE>(elem->type(), fet.order + add_p_level*elem->p_level(), i, j, p);}
// #endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
// } // namespace libMesh
// 
// 
// 


// ================================================================================================
template <> 
double MGFE<2,LAGRANGE>::ComputeInverseMatrix_vol(double Matrix[], double InvMatrix[]
  
) {// ================================================================================================
  double det = (Matrix[0] * Matrix[3] - Matrix[1] * Matrix[2]);
  double idet = 1. / det;
  InvMatrix[0] = Matrix[3] * idet;   // dxi dx
  InvMatrix[1] = -Matrix[2] * idet;  // deta dx
  InvMatrix[2] = -Matrix[1] * idet;  // dxi dy
  InvMatrix[3] = Matrix[0] * idet;   // deta dy
  return det;
}
//===============================================================================================
template <> 
double MGFE<2,LAGRANGE>::ComputeInverseMatrix_sur(double Matrix[], double InvMatrix[]) {
  double det = (Matrix[0]);
  double idet = 1. / det;
  InvMatrix[0] =  idet;   // dxi dx
  return det;
}

//===============================================================================================
template <> 
double MGFE<2,LAGRANGE>::ComputeInverseMatrix(double Matrix[], double InvMatrix[], int dim) {
  double det = 0.;
  
  if(dim == _dim)   det = ComputeInverseMatrix_vol(Matrix,InvMatrix);
  if(dim == _dim-1) det = ComputeInverseMatrix_sur(Matrix,InvMatrix);
  
  return det;
}



//===============================================================================================
template <> 
double MGFE<2,LAGRANGE>::Tj_interp_lin(double *u,int j, double *xi,int ndof,int idim) {
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
template <> 
double MGFE<2,LAGRANGE>::dTjdxi_lin(double *u,int j, double *xi,int i,int ndof,int idim) {
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


// ================================================
/// This function computes the normal at the gauss point
template <>
void MGFE<2,LAGRANGE>::normal_g(
const double* xx,  // all surface coordinates <-
double* normal_g   // normal ->
) const {          // ======================================
 const int NDOF_FEMB=3;
  // coordinates
  double xx3D[3 * NDOF_FEMB];
  for(int i = 0; i < _dim * NDOF_FEMB; i++) { xx3D[i] = xx[i]; }
  double tg01[3];  // tangent  line
    //  2D -----------------------------------------------

  //  The surface elements are such that when you go from the 1st to
  //   the 2nd point, the outward normal is to the RIGHT
  for(int i = 0; i < 2; i++) { tg01[i] = xx3D[1 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB]; }
  //  Rotation matrix (0 1; -1 0), 90deg clockwise
  normal_g[0] = 0. * tg01[0] + 1. * tg01[1];
  normal_g[1] = -1. * tg01[0] + 0. * tg01[1];
  // normalization -----------------------
  double mm = 0.;  for(int idim=0;idim<2;idim++) {mm += normal_g[idim]*normal_g[idim];}
  mm = std::sqrt(mm);  for(int idim=0;idim<2;idim++) {normal_g[idim] /= mm;}

  return;
}

//=============================================================
template <>
void MGFE<2,LAGRANGE>::normal_g(
const double xx[],  // all surface coordinates <-
const double x_c[], ///< central point
double* normal_g   // normal ->
) const {          // ======================================
const int NDOF_FEMB=3;
  // coordinates
double xx3D[3 * NDOF_FEMB];double tg01[3]; tg01[2]=0.;
for(int i = 0; i <2*NDOF_FEMB; i++) { xx3D[i] = xx[i]; }
  //  The surface elements are such that when you go from the 1st to
  //   the 2nd point, the outward normal is to the RIGHT
  for(int i = 0; i < 2; i++) { tg01[i] = xx3D[1 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB]; }
  //  Rotation matrix (0 1; -1 0), 90deg clockwise
  normal_g[0] = 0. * tg01[0] + 1. * tg01[1];
  normal_g[1] = -1. * tg01[0] + 0. * tg01[1];
  // if the sign is not correct then reverse it
  if(normal_g[0] * (x_c[0] - xx3D[0]) + normal_g[1] * (x_c[1] - xx3D[0 + NDOF_FEMB]) > 0.) {
                 normal_g[0] *= -1.; tg01[0] *= -1.;   normal_g[1] *= -1.; tg01[1] *= -1.;
    //    std::cout << " Normal inverted ! ------------------------------------  \n";
  }
  // normalization -----------------------
  double mn = 0.; double mtg = 0.;
  for(int idim=0;idim <2;idim++) {mn += normal_g[idim]*normal_g[idim]; mtg +=tg01[idim]*tg01[idim];}
  mn = sqrt(mn); mtg = sqrt(mtg);  for(int idim =0;idim<2;idim++) {normal_g[idim] /= mn;  tg01[idim] /=mtg; }

  return;
}
//=============================================================
template <>
void MGFE<2,LAGRANGE>::normal_g(
const double xx[],  // all surface coordinates <-
const double x_c[], ///< central point
double* normal_g,   // normal ->
int &sign
) const {          // ======================================
const int NDOF_FEMB=3;
  // coordinates
double xx3D[3 * NDOF_FEMB];double tg01[3]; tg01[2]=0.;
for(int i = 0; i <2*NDOF_FEMB; i++) { xx3D[i] = xx[i]; }
  //  The surface elements are such that when you go from the 1st to
  //   the 2nd point, the outward normal is to the RIGHT
  for(int i = 0; i < 2; i++) { tg01[i] = xx3D[1 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB]; }
  //  Rotation matrix (0 1; -1 0), 90deg clockwise
  normal_g[0] = 0. * tg01[0] + 1. * tg01[1];
  normal_g[1] = -1. * tg01[0] + 0. * tg01[1];
  // if the sign is not correct then reverse it
  if(normal_g[0] * (x_c[0] - xx3D[0]) + normal_g[1] * (x_c[1] - xx3D[0 + NDOF_FEMB]) > 0.) {
                 normal_g[0] *= -1.; tg01[0] *= -1.;   normal_g[1] *= -1.; tg01[1] *= -1.;
    //    std::cout << " Normal inverted ! ------------------------------------  \n";
  }
  // normalization -----------------------
  double mn = 0.; double mtg = 0.;
  for(int idim=0;idim <2;idim++) {mn += normal_g[idim]*normal_g[idim]; mtg +=tg01[idim]*tg01[idim];}
  mn = sqrt(mn); mtg = sqrt(mtg);  for(int idim =0;idim<2;idim++) {normal_g[idim] /= sign*mn;  tg01[idim] /=sign*mtg; }

  return;
}

//=============================================================
template <> void  MGFE<2,LAGRANGE>::Oxy_face(
const double* xx,   ///< coordinates
const double x_c[], ///< central point
double* normal_g,   ///<  i_n,i_t1,i_t2 face Oxyz
double* tg01,
double* /*tg03*/
) const {            // ======================================
const int NDOF_FEMB=3;
  // coordinates
double xx3D[3 * NDOF_FEMB];
for(int i = 0; i <2*NDOF_FEMB; i++) { xx3D[i] = xx[i]; }

//   double tg01[MGFE1::MDIM_FEM];  // tangent  line


  //  The surface elements are such that when you go from the 1st to
  //   the 2nd point, the outward normal is to the RIGHT
  for(int i = 0; i < 2; i++) { tg01[i] = xx3D[1 + i * NDOF_FEMB] - xx3D[0 + i * NDOF_FEMB]; }

  //  Rotation matrix (0 1; -1 0), 90deg clockwise
  normal_g[0] = 0. * tg01[0] + 1. * tg01[1];
  normal_g[1] = -1. * tg01[0] + 0. * tg01[1];
  // if the sign is not correct then reverse it

  if(normal_g[0] * (x_c[0] - xx3D[0]) + normal_g[1] * (x_c[1] - xx3D[0 + NDOF_FEMB]) > 0.) {
    normal_g[0] *= -1.; tg01[0] *= -1.;   normal_g[1] *= -1.; tg01[1] *= -1.;
    //    std::cout << " Normal inverted ! ------------------------------------  \n";
  }
  // normalization -----------------------
  double mn = 0.; double mtg = 0.;
  for(int idim = 0; idim < _dim; idim++) {mn += normal_g[idim] * normal_g[idim]; mtg +=tg01[idim]*tg01[idim];}
  mn = sqrt(mn); mtg = sqrt(mtg);  for(int idim =0;idim<2;idim++) {normal_g[idim] /= mn;  tg01[idim] /=mtg; }

  return;
}

// =================================================================================
template <> void MGFE<2,LAGRANGE>::set_par(
  std::string fileh5,
  int dim_in,
  int order_in,
  int shapes_quad,
  int geom_type_in,
  int fam_type_in
) {//================================================================================
  _dim=dim_in; 
  // _order= order_in;  /*_NoShape=shapes_quad, */  
  _GeomType=geom_type_in;// 0 ->TRI
  
  // 1=HEX27
  _FamType=0;//Lagrangian 
   femfileh5_=fileh5;
     ///< family type  0=Lagrangian(C^0) 1=Raviart Thomas (L^2-C^0)
  // Fem order
  // =============================
  //  ------------------------------------------------------------------------
  // gauss points (_NoGauss) ------------------------------------------------------
  _NoGauss=0;  assert(_qrule); _NoGauss=_qrule->n_points();
  if(_NoGauss==0) {   std::cout << "Quadrature not assigned"; abort();}

  // number of nodes (_NoElNodes) in the reference quad element
   ElemType  elem_name=_fe_type.fetype; 
   ElemType  geom_elem_name=get_ref_elemtype(elem_name); 
   Order elem_order= _fe_type.order ;
   _order= (int)elem_order; 
  // number of nodes (_NoElNodes) in the reference EDGE3 element
  _NoElNodes= get_refspace_n_nodes(geom_elem_name);// mesh nodes
  _NoShape=get_refspace_n_nodes(elem_name); // element nodes 
  

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
  // if(_FamType==0) {  // 0= Lagrangian  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    // lagrangian -> n_point= n_shape =_order DDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    //            -> pol degree=_order/2 (int div)
    // _deg=0;// DDDDDDDDDDDDDDDDDDDDDDDDD

 // dofs -----------------------------------------------------------------------
 for(int j = 0; j < _NoShape; j++) _fem_conn[j]=j;
 double *xnodes=new double [_NoElNodes*_dim];
 double *point=new double[_dim];
    // ---------------------------------------------------------------------------------------
    // nodes ----------------------------------------------------------------------------
  get_refspace_nodes(elem_name,xnodes);// nodes[(x,y)_0,(x,y)_1.... ]
    // _xnodes[0]=0.;    _xnodes[1]=1.;    _xnodes[2]=0.5;
    for(int i = 0; i < _NoElNodes; i++) {  // points
      for(int dir = 0; dir < _dim; dir++) { 
        _xnodes[i+_NoElNodes*dir]=xnodes[dir+_dim*i]; 
        point[dir] = _xnodes[i+_NoElNodes*dir];
      }
    for(int j = 0; j < _NoShape; j++) {  // shapse
      //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
      
      _phi_map_nodes[j*_NoElNodes+i]=fe_lagrange_2D_shape<LAGRANGE>(elem_name , elem_order, j, point);//Rec_Quad_Phi(j, point,_dim);
     
      // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
  
    for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
        _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =fe_lagrange_2D_shape_deriv<LAGRANGE>(elem_name, elem_order, j,dir, point);// Rec_Quad_DPhi(j,point,_dim,dir);     // dphi/dxi
      }
    }
  }
   delete[]xnodes;
// gaussian points --------------------------------------------------------->

 for(int i = 0; i < _NoGauss; i++) {
    _weight[i] =_qrule->w(i); // weight
     _qrule->qp(point,i);     // ith-gaussian point ->point(i)
   
    for(int j = 0; j < _NoShape; j++) {     // loop over test function id
      _phi_map[j*_NoGauss+i] =fe_lagrange_2D_shape<LAGRANGE>(elem_name , elem_order, j, point);
        int flag_second_deriv=0;
      for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
        _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =fe_lagrange_2D_shape_deriv<LAGRANGE>(elem_name, elem_order, j,dir, point);
            // Order flag_sec_deriv= 0->d^2()/dxi^2; 1->d^2()/dxideta;2->d^2()/deta^2
        for(int dir2 = 0; dir2 <=dir; dir2++){//old ->lag_second_deriv ->dir*_dim+dim2
          _dphidxx_map[(j + flag_second_deriv * _NoShape) * _NoGauss + i] =
          fe_lagrange_2D_shape_second_deriv<LAGRANGE>(elem_name, elem_order,j,flag_second_deriv, point);
          flag_second_deriv++;
        }
      }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
  } 
  delete []point;

#ifdef PRINT_INFO
  std::cout << "\n MGFE: elem_name(Dim): "           <<elem_name<<"(" <<_dim <<")" ;
  std::cout << "\n       gaussian points: "          << _NoGauss ;
  std::cout << "\n       number of shape functions " << _NoShape ;
  std::cout << "\n       polynomial order "          << elem_order;
  std::cout << " \n";  
  
#endif
  return;
}

// ========================================
// ========================================
//   JAC
// ========================================
// ========================================
template <> double   MGFE<2,LAGRANGE>::compute_fe_face_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz_bd[] // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)]  3D point
 ){
  //    ElemType  elem_name=_fe_type.fetype; 
  //    ElemType  face_name=get_ref_face_elemtype(elem_name); 
  // //  Order elem_order= _fe_type.order ;
  // //   _order= (int)elem_order;
  // // // number of nodes (_NoElNodes) in the reference EDGE3 element
  // // _NoElNodes= get_refspace_n_nodes( geom_elem_name);// mesh nodes
  //  int NoShape_face=get_refspace_n_nodes(face_name); // element nodes 
  
  for(int s=0;s<_NoShape;s++)   {
      phi[s]=_phi_map[s* _NoGauss  +gp] ;// shape
      dphidxi[s]=_dphidxez_map[s*_NoGauss + gp];  dphi[s]=dphidxi[s];
      dphideta[s]=_dphidxez_map[(s+_NoShape)*_NoGauss + gp]; dphi[s+_NoShape]=dphideta[s];
      // dphidzeta[s]=_dphidxez_map[(s+2*_NoShape)*_NoGauss + gp]; dphi[s+2*_NoShape]=dphidzeta[s];
  }
  // if(!dphidzeta)  for(int l=0;l<_NoShape;l++) {dphidzeta[l]=_dphidxez_map[(l+1*_NoShape)*_NoGauss + gp];phi[l+2*_NoShape]=dphidzeta[l];}
    // Values to compute at gaussian points
  double dxdxi = 0.;  double dxdeta = 0.;  //
  double dydxi = 0.;  double dydeta = 0.;  // d(x,y,z)d(xi,eta)
  double dzdxi = 0.;  double dzdeta = 0.;  //
 double InvJac[4];
  //   int nshape=_NoShape[1];
  // int offset = _NoShape * _NoGauss;

  for(int s = 0; s < (int)_NoShape; s++) {
    // int sng = s * _NoGauss + gp;
    // double dphidxi = _dphidxez_map[ s * _NoGauss + gp];
    // double dphideta = _dphidxez_map[ s * _NoGauss + gp + _NoShape * _NoGauss];

    dxdxi +=  xyz_bd[s] * dphidxi[s];
    dxdeta += xyz_bd[s] * dphideta[s];
    dydxi +=  xyz_bd[s +_NoElNodes] * dphidxi[s];
    dydeta += xyz_bd[s +_NoElNodes] * dphideta[s];
    dzdxi +=  xyz_bd[s +2*_NoElNodes] * dphidxi[s];
    dzdeta += xyz_bd[s +2*_NoElNodes] * dphideta[s];
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
  
   for(int l=0;l<_NoShape;l++) {
  dphidx[l]=InvJac[0] *dphidxi[l]+ InvJac[1] *dphideta[l];
  dphidy[l]=InvJac[2] *dphidxi[l]+ InvJac[3] *dphideta[l];
  }  
  return  0.;
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
template <> void MGFE<2,LAGRANGE>::get_ddphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double ddphi_g[]          // global derivatives ->
) {                         // =========================================
  
  //Hessian=JacI*H_loc*jacI^T
    // const int el_nnodes = _NoShape;             // # of shape functions
  // const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset =  _NoShape * _NoGauss;  // gauss offset
  double ddphi_loc[36];
  double Hess_tmp[4], JacTf[4];
  
  for(int eln = 0; eln < _NoShape; eln++) {   // LOOP OVER NODES ===================================
    const int shift = eln * _dim * _dim;       // offset for node derivatives
    const int lqp = eln * _NoGauss + qp;
     ddphi_loc[shift +0] = _dphidxx_map[(eln +0) * _NoGauss + qp];//d2d2xi
     ddphi_loc[shift +1] = _dphidxx_map[(eln +1*_NoShape) * _NoGauss + qp];//d2dxideta
     ddphi_loc[shift +2] = _dphidxx_map[(eln +1*_NoShape) * _NoGauss + qp];//d2dxideta
     ddphi_loc[shift +3] = _dphidxx_map[(eln +2*_NoShape) * _NoGauss + qp];//d2d2eta
    // for(int idim = 0; idim < _dim * _dim; idim++) {  
    //   // ddphi_loc[shift + idim] = _dphidxx_map[lqp + idim * goffset];
    //    ddphi_loc[shift + idim] = _dphidxx_map[(eln +_dim*_NoShape) * _NoGauss + gp];
    // }//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    for(int init = 0; init < _dim; init++)
      for(int jnit = 0; jnit < _dim; jnit++) {
        Hess_tmp[init * _dim + jnit] = 0.;
        ddphi_g[shift + init * _dim + jnit] = 0.;
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
          ddphi_g[shift + ii * _dim + jj] += InvJac[ii * _dim + ss] * Hess_tmp[ss * _dim + jj];
        }
  }  // END LOOP OVER NODES =====================================================================

  // Pay attention, different order with respect to first order derivatives!
  //   delete[]ddphi_loc;
  
  

  //   delete[]ddphi_loc;
  return;
}
// ========================================
/// This function computes the  derivatives at
/// the gaussian point ng:
template <> double MGFE<2,LAGRANGE>::Jac(    //  jacobean ->
const int ng,    // gaussian point <-
const double xyz[],      // coordinates  <-  xyz[x(_NoElNodes),y(_NoElNodes)]
double InvJac[]  // Inverted Jacobean ->
) {                  // =====================================
  
    double x_xi = 0., x_eta = 0., y_xi = 0., y_eta = 0.;
    for(int l=0;l<_NoShape;l++)   {
      phi[l]=_phi_map[l* _NoGauss  +ng] ;// shape
      dphidxi[l]=_dphidxez_map[l*_NoGauss + ng];  
      dphideta[l]=_dphidxez_map[(l+_NoShape)*_NoGauss + ng];
      x_xi  += xyz[l]*dphidxi[l];
      x_eta += xyz[l]*dphideta[l];
      y_xi  += xyz[l+_NoElNodes]*dphidxi[l];  // -> _NoElNodes
      y_eta += xyz[l+_NoElNodes]*dphideta[l];
   }
  
  
 // int NDOF_FEM=9;
 // if(_GeomType==0) NDOF_FEM=6;
 //  double Jac = 0.;
 //  double x_xi = 0., x_eta = 0., y_xi = 0., y_eta = 0.;
 //  int nshape = _NoShape;
 //  int offset = _NoShape * _NoGauss;
 // 
 //  for(int s = 0; s < nshape; s++) {
 //    int sng = s * _NoGauss + ng;
 //    x_xi += x[s] * _dphidxez_map[sng];
 //    x_eta += x[s] * _dphidxez_map[sng + offset];
 //    y_xi += x[s + _NoShape] * _dphidxez_map[sng];
 //    y_eta += x[s +_NoShape] * _dphidxez_map[sng + offset];
 //  }

  double det = x_xi * y_eta - y_xi * x_eta;
  double idet = 1. / det;
  InvJac[0] = y_eta * idet;   // dxi dx
  InvJac[1] = -y_xi * idet;   // deta dx
  InvJac[2] = -x_eta * idet;  // dxi dy
  InvJac[3] = x_xi * idet;    // deta dy

  return (det);
}

//==========================================================================================================
template <> double   MGFE<2,LAGRANGE>::compute_fe_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz[], // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_NoElNodes)] 
  double InvJac[]
 ){
     // Compute the shape function values (and derivatives)
//   // at the Quadrature points qp. 
   
  //  ElemType  elem_name=_fe_type.fetype; 
  //  ElemType  geom_elem_name=get_ref_elemtype(elem_name); 
  //  Order elem_order= _fe_type.order ;
  //   _order= (int)elem_order;
  // // number of nodes (_NoElNodes) in the reference EDGE3 element
  // _NoElNodes= get_refspace_n_nodes( geom_elem_name);// mesh nodes
  // _NoShape=get_refspace_n_nodes(elem_name); // element nodes 
  // double InvJac[4];
   // double det = Jac(gp,xyz,InvJac); 
    
//   double x_xi = 0., x_eta = 0., y_xi = 0., y_eta = 0.;
//   
     for(int l=0;l<_NoShape;l++)   {
       phi[l]=_phi_map[l* _NoGauss  +gp] ;// shape
       dphidxi[l]=_dphidxez_map[l*_NoGauss + gp];  
       dphideta[l]=_dphidxez_map[(l+_NoShape)*_NoGauss + gp];
//       x_xi  += xyz[s]*dphidxi[s];
//       x_eta += xyz[s]*dphideta[s];
//       y_xi  += xyz[s+_NoElNodes]*dphidxi[s];  // -> _NoElNodes
//       y_eta += xyz[s+_NoElNodes]*dphideta[s];
       }
  // if(!dphidzeta)  for(int l=0;l<_NoShape;l++) {dphidzeta[l]=_dphidxez_map[(l+1*_NoShape)*_NoGauss + gp];phi[l+2*_NoShape]=dphidzeta[l];}
  

  // int nshape = _NoShape;
  // int offset = _NoShape * _NoGauss;

  // for(int s = 0; s < _NoShape; s++) {
  //   x_xi  += xyz[s]*dphidxi[s];
  //   x_eta += xyz[s]*dphideta[s];
  //   y_xi  += xyz[s+_NoShape]*dphidxi[s];
  //   y_eta += xyz[s+_NoShape]*dphideta[s];
  // }

  // double det = x_xi * y_eta - y_xi * x_eta;
  // double idet = 1. / det;
  // InvJac[0] = y_eta * idet;   // dxi dx
  // InvJac[1] = -y_xi * idet;   // deta dx
  // InvJac[2] = -x_eta * idet;  // dxi dy
  // InvJac[3] = x_xi * idet;    // deta dy
   
  // dphidxi dphideta from Jac()
   for(int l=0;l<_NoShape;l++) {
  dphidx[l]=InvJac[0] *dphidxi[l]+ InvJac[1] *dphideta[l];  dphi[l]=dphidx[l];
  dphidy[l]=InvJac[2] *dphidxi[l]+ InvJac[3] *dphideta[l];  dphi[l+_NoShape]=dphidy[l];
  
  // // secondorder derivatives // Order flag_sec_deriv= 0->d^2()/dxi^2; 1->d^2()/dxideta;2->d^2()/deta^2
  // // d2phi[l*_dim*_dim+dxi1*_dim+dxi2]=d^2phi/dxi1 dxi2
  //  d2phi[l*_dim*_dim+0*_dim+0]= _dphidxx_map[(l + 0 * _NoShape) * _NoGauss + gp];
  //   d2phi[l*_dim*_dim+0*_dim+0]=InvJac[0] *dphidxi[l]+ InvJac[1] *dphideta[l];  dphi[l]=dphidx[l];
  //  d2phi[l*_dim*_dim+1*_dim+0]= d2phi[l*_dim*_dim+0*_dim+1]= _dphidxx_map[(l + 1 * _NoShape) * _NoGauss + gp];
  //  d2phi[l*_dim*_dim+1*_dim+1]= _dphidxx_map[(l + 2 * _NoShape) * _NoGauss + gp];
  }  
  
  // 2nd derivative ------------------------------------------------------>
   /// Shape functions 2nd derivatives at the gauss point qp
/// d2phi[id*_dim*_dim+i*_dim+j]= Tij
///    Txx Txy Txz
///    Tyx Tyy Tyz
///    Tzx Tzy Tzz
  get_ddphi_gl_g(gp,InvJac,d2phi);          // global derivatives ->
  

  
  
  
  
  
  return 0;
}


// ========================================
/// This function computes the  derivatives at
/// the nodal point:
template <> double MGFE<2,LAGRANGE>::Jac_nodes( 
// double MGFE1::Jac_nodes(  //  jacobean ->
const int ng,        // nodal point <-
const double x[],          // coordinates  <-
double InvJac[]      // Jacobean ->
) {                      // =====================================
 int NDOF_FEM=9;
 if(_GeomType==0) NDOF_FEM=6;

  double x_xi = 0., x_eta = 0., y_xi = 0., y_eta = 0.;
  int nshape = _NoShape;
  int offset = _NoShape * _NoGauss;

  for(int s = 0; s < nshape; s++) {
    int sng = s * _NoGauss + ng;
    x_xi += x[s] * _dphidxez_map_nodes[sng];
    x_eta += x[s] * _dphidxez_map_nodes[sng + offset];
    y_xi += x[s + NDOF_FEM] * _dphidxez_map_nodes[sng];
    y_eta += x[s + NDOF_FEM] * _dphidxez_map_nodes[sng + offset];
  }

  //   std::cout<<x_xi<<"  "<<x_eta<<"  "<<y_xi<<"  "<<y_eta<<std::endl;
  //
  double det = x_xi * y_eta - y_xi * x_eta;
  double idet = 1. / det;
  InvJac[0] = y_eta * idet;   // dxi dx
  InvJac[1] = -y_xi * idet;   // deta dx
  InvJac[2] = -x_eta * idet;  // dxi dy
  InvJac[3] = x_xi * idet;    // deta dy

  return (det);
}

// ==========================================
template <> double MGFE<2,LAGRANGE>::JacSur(
// double MGFE1::JacSur(  // surface jacobean ->
const int ng,     // gaussian point
const double x[],       // coordinates
double InvJac[]   // Jacobean ->
)  {         // ================================
   int NDOF_FEMB=3;
 
  // Values to compute at gaussian points
  double dxdxi = 0.;  //
  double dydxi = 0.;  // d(x,y)d(xi,eta)
  const int nshape = _NoShape;

  //   int offset=_NoShape[0]*_NoGauss[0];
  for(int s = 0; s < nshape; s++) {
    int sng = s * _NoGauss + ng;
    double dphidxi = _dphidxez_map[sng];
    dxdxi += x[s] * dphidxi;
    dydxi += x[s + NDOF_FEMB] * dphidxi;
  }

  // surface weighted jacobean
  double det = sqrt(dxdxi * dxdxi + dydxi * dydxi);
  InvJac[0] = 1. / det;
  return (det);
}



// ========================================
/// This function computes the 3D derivatives at
/// the gaussian point ng. It computes the Jacobian (f=Jv+b  J=(nabla_v f)^T)
template <> double MGFE<2,LAGRANGE>::JacobianOnGauss(  // 3D  jacobean ->
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

  // if(dim==2) 
    det = (Jac[0]*Jac[3]-Jac[1]*Jac[2]);

  // if(dim==3) det = Jac[0]*(Jac[4]*Jac[8]-Jac[5]*Jac[7])
  //                    -Jac[3]*(Jac[1]*Jac[8]-Jac[2]*Jac[7])
  //                    +Jac[6]*(Jac[1]*Jac[5]-Jac[2]*Jac[4]);

  return (det);
}

// ========================================================================
template <> double MGFE<2,LAGRANGE>::JacOnGivenCanCoords(
const int dim,const double ElemCoords[], double CanCoords[], double InvJac[], int FamilyType, int nShape) {
//   double* LocDPhi = new double[dim];
//   double* CoordsDer = new double[dim * dim];
// 
//   for(int dir = 0; dir < dim * dim; dir++) { CoordsDer[dir] = 0.; }
// 
//   for(int s = 0; s < nShape; s++) {
//     for(int dir = 0; dir < dim; dir++) {
//       LocDPhi[dir] = FirstDerivateOfLocalPhi(s, CanCoords, dim, dir, FamilyType);
//     }
// 
//     for(int dir1 = 0; dir1 < dim; dir1++)
//       for(int dir2 = 0; dir2 < dim; dir2++) {
//         CoordsDer[dir1 * dim + dir2] += ElemCoords[s + dir1 * nShape] * LocDPhi[dir2];
//       }
//   }
// 
//   double det = 0.;
//  det = (CoordsDer[0] * CoordsDer[3] - CoordsDer[1] * CoordsDer[2]); 
// 
//   double idet = 1. / det;
//     InvJac[0] = CoordsDer[3] * idet;   // dxi dx
//     InvJac[1] = -CoordsDer[2] * idet;  // deta dx
//     InvJac[2] = -CoordsDer[1] * idet;  // dxi dy
//     InvJac[3] = CoordsDer[0] * idet;   // deta dy
//  
//     for(int row = 0; row < dim; row++) {   // LOOP OVER ROWS - row1<row2
//       int sign_row = 1 - 2 * (row % 2);
//       int row1 = ((row + 1) % dim < (row + 2) % dim) ? (row + 1) % dim : (row + 2) % dim;
//       int row2 = ((row + 1) % dim > (row + 2) % dim) ? (row + 1) % dim : (row + 2) % dim;
// 
//       for(int col = 0; col < dim; col++) {   // LOOP OVER COLUMNS  - col1<col2
//         int col1 = ((col + 1) % dim < (col + 2) % dim) ? (col + 1) % dim : (col + 2) % dim;
//         int col2 = ((col + 1) % dim > (col + 2) % dim) ? (col + 1) % dim : (col + 2) % dim;
//         int sign_col = 1 - 2 * (col % 2);
//         InvJac[row * dim + col] = sign_col * sign_row * idet *
//                                   (CoordsDer[row1 * dim + col1] * CoordsDer[row2 * dim + col2] -
//                                    CoordsDer[row1 * dim + col2] * CoordsDer[row2 * dim + col1]);
//       }
//     }
// 
//   delete[] LocDPhi;
//   delete[] CoordsDer;
//   return det;
  return 0.;
}
// ================================================================
template <> double MGFE<2,LAGRANGE>::JacobianOnPoint(
const int dim,const double xyz[], double ElemCoords[], double Jac[], int FamilyType, int nShape) {

//   double LocDPhi[2];
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
// //   if(dim == 3)  //{
// //     det =   Jac[0]*(Jac[4] * Jac[8] - Jac[5] * Jac[7])
// //             -Jac[3]*(Jac[1] * Jac[8] - Jac[2] * Jac[7])
// //             +Jac[6]*(Jac[1] * Jac[5] - Jac[2] * Jac[4]);
// // 
// // //     for(int dir = 0; dir < dim; dir++) {
// // //       int sign = 1 - 2 * (dir % 2);
// // //       int idx1 = (dir + 1) % dim;
// // //       int idx2 = (dir + 2) % dim;
// // //       det += sign *JacT[dir * dim + dir]* (JacT[idx1 * dim + idx1] * JacT[idx2 * dim + idx2] -
// // //                      JacT[idx1 * dim + idx2] * JacT[idx2 * dim + idx1]);
// // //     }
// // //   }
//   // if(dim == 2) 
//   { det = (Jac[0] * Jac[3] - Jac[1] * Jac[2]); }
// 
//   return det;
  return 0.;
}



// =========================================================
//  GET phi dphi ddphi
//==========================================================


//==============================================================================================
template <> void MGFE<2,LAGRANGE>::get_dphi_on_given_node(const int dim, double ElemCoords[], double CanPos[], double dphi[]) {
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
template <> void MGFE<2,LAGRANGE>::get_phi_gl_g(
const int qp,    // gaussian point <-
double phi_g[]     // shape functions ->
) {                  // =================================================
  for(int ish = 0; ish < _NoShape; ish++) phi_g[ish] = _phi_map[ish * _NoGauss + qp];
  return;
}

// ==================================================================
/// This function computes the shape values at the gauss point qp
template <> void MGFE<2,LAGRANGE>::get_phi_gl_g(
const int qp,             // gaussian point   <-
std::vector<double>& phi_g  // shape functions  ->
) {                           // =====================================
  for(int ish = 0; ish < _NoShape; ish++) phi_g[ish] = _phi_map[ish * _NoGauss + qp];

  return;
}

// =================================================================
/// Shape functions derivatives dphi[id+i* el_nnodes] =Ti
/// Tx Ty Tz
///  tensor  order  at the gauss point qp, node id:
///  we have  dphi[Tx(0),Tx(1),..Tx(id), Ty(0),Ty(1),..Ty(id), Tz(0),Tz(1),..Tz(id) ]
template <> void MGFE<2,LAGRANGE>::get_dphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double dphi_g[]           // global derivatives ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  double gradphi_g[2];  // temp grad phi

  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + qp;
    for(int idim = 0; idim < _dim; idim++) gradphi_g[idim]=_dphidxez_map[lqp+idim*goffset];

    for(int idim = 0; idim < _dim; idim++) {
      double sum=0.; for(int jdim=0;jdim<_dim;jdim++) sum +=InvJac[jdim+idim*_dim]*gradphi_g[jdim];
      dphi_g[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}


// =============================================================
/// Shape functions derivatives at the nodal points
template <> void MGFE<2,LAGRANGE>::get_dphi_node(
const int /*kdim*/,         // dimension <-
const int node,         // nodal point <-
const double InvJac[],  // Jacobean <-
double dphi_nodes[]           // global derivatives at the nodal point ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  double gradphi_g[2];  // temp grad phi

  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + node;

    for(int idim = 0; idim < _dim; idim++) {
      gradphi_g[idim] = _dphidxez_map_nodes[lqp + idim * goffset];
    }

    for(int idim = 0; idim < _dim; idim++) {
      double sum = 0.;

      for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * gradphi_g[jdim]; }

      dphi_nodes[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}


// =============================================================
/// Shape functions derivatives at the gauss point qp
template <> void MGFE<2,LAGRANGE>::get_dphi_gl_g(
const int /*kdim*/,         // dimension <-
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double dphi_g[],          // global derivatives ->
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
      double sum = 0.;

      for(int jdim = 0; jdim < sdim; jdim++) { sum += InvJac[jdim + idim * sdim] * gradphi_g[jdim]; }

      dphi_g[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

// =============================================================
/// Shape functions derivatives at the gauss point qp
template <> void MGFE<2,LAGRANGE>::get_dphi_gl_g(
const int /*kdim*/,            // dimension <-
const int qp,              // gaussian point <-
const double InvJac[],     // Jacobean
std::vector<double>& dphi_g  // global derivatives ->
) {                            // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of guass points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  double dphidxi_g[2];  // temp grad phi

  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + qp;

    for(int idim = 0; idim < _dim; idim++) { dphidxi_g[idim] = _dphidxez_map[lqp + idim * goffset]; }

    for(int idim = 0; idim < _dim; idim++) {
      double sum = 0.;

      for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * dphidxi_g[jdim]; }

      dphi_g[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}
//=========================================================================================================
template <>  bool MGFE<2,LAGRANGE>::phys_to_ref(double xref[],double xphys[],const double Xel[], int max_it, double tol) {
  // std::array<double,2> xi = {0.0, 0.0};

    auto norm2 = [](const std::array<double,2>& v){
        return v[0]*v[0] + v[1]*v[1];
    };

    for (int it = 0; it < max_it; ++it) {
        std::array<double,2> xmap = {0.0, 0.0};
        double J[2][2] = {{0.0,0.0},{0.0,0.0}};

        double p[2] = {xref[0], xref[1]};

        for (int i = 0; i < _NoShape; ++i) {
            // const double Ni = Rec_Quad_Phi(i, p, /*_dim=*/2);
             const double Ni = fe_lagrange_2D_shape<LAGRANGE>(QUAD9,SECOND,i, p);
            // const double dN_dxi  = Rec_Quad_DPhi(i, p, /*_dim=*/2, /*dir=*/0);
              const double dN_dxi  = fe_lagrange_2D_shape_deriv<LAGRANGE>(QUAD9,SECOND, i,0, p);
            // const double dN_deta = Rec_Quad_DPhi(i, p, /*_dim=*/2, /*dir=*/1);
            const double dN_deta = fe_lagrange_2D_shape_deriv<LAGRANGE>(QUAD9,SECOND, i,1, p);

            const double Xi = Xel[i + _NoShape*0];
            const double Yi = Xel[i + _NoShape*1];

            xmap[0] += Ni * Xi;
            xmap[1] += Ni * Yi;

            J[0][0] += dN_dxi  * Xi;  // dx/dxi
            J[0][1] += dN_deta * Xi;  // dx/deta
            J[1][0] += dN_dxi  * Yi;  // dy/dxi
            J[1][1] += dN_deta * Yi;  // dy/deta
        }

        std::array<double,2> r = {xmap[0] - xphys[0], xmap[1] - xphys[1]};

        const double r2 = norm2(r);
        const double scale = 1.0 + xphys[0]*xphys[0] + xphys[1]*xphys[1];
        if (r2 <= tol*tol*scale) return true;

        const double det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        if (std::abs(det) < 1e-30) {
            return false;
        }

        std::array<double,2> d;
        d[0] = ( r[0]*J[1][1] - r[1]*J[0][1]) / det;
        d[1] = ( J[0][0]*r[1] - J[1][0]*r[0]) / det;

        double alpha = 1.0;
        double xi_try[2];
        for (int ls = 0; ls < 10; ++ls) {
            xi_try[0] =  xref[0] - alpha*d[0];
            xi_try[1] =  xref[1] - alpha*d[1];


            double p2[2] = {xi_try[0], xi_try[1]};
            std::array<double,2> x2 = {0.0,0.0};
            for (int i = 0; i < _NoShape; ++i) {
                // const double Ni = Rec_Quad_Phi(i, p2, 2);
              const double Ni = fe_lagrange_2D_shape<LAGRANGE>(QUAD9,SECOND,i, p2);;
                x2[0] += Ni * Xel[i + _NoShape*0];
                x2[1] += Ni * Xel[i + _NoShape*1];
            }
            std::array<double,2> r_try = {x2[0]-xphys[0], x2[1]-xphys[1]};

            if (norm2(r_try) < r2) break;
            alpha *= 0.5;
        }

        xref = xi_try;
    }

    return false;
}


//=================================================================================
// // auxiliary function definitions
// using namespace libMesh;
// =========================================================================
template <FEFamily T> double fe_lagrange_2D_shape(
  const ElemType type,  // QUAD4 etc .... 
  const Order order,    // first second third ...
  const unsigned int i, //  shape(i)
  const double  p[]     // point[]={x,y}
){//==========================================================================
// #if LIBMESH_DIM > 1
  switch (order) {
    case CONSTANT: {     return 1.;    }
    case FIRST: { // linear Lagrange shape functions
        switch (type) {
          case QUAD4:    case QUADSHELL4:   case QUAD8:    case QUADSHELL8:   case QUAD9:  case QUADSHELL9:   {
              // Compute quad shape functions as a tensor-product
              const double xi  = p[0];  const double eta = p[1];
              assert (i< 4);
              //                                0  1  2  3
              static const unsigned int i0[] = {0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1};
              return (fe_lagrange_1D_linear_shape(i0[i], xi)* fe_lagrange_1D_linear_shape(i1[i], eta));
            }
          case TRI3:      case TRISHELL3:      case TRI6:     case TRI7:   {
              const double zeta1 = p[0];   const double zeta2 = p[1];
              const double zeta0 = 1. - zeta1 - zeta2;
              assert(i< 3);
              switch(i) {
                case 0: return zeta0;
                case 1: return zeta1;
                case 2: return zeta2;
                default: std::cerr<<"Invalid shape function index i = " << i;
                }
            }
          default:  std::cerr<<"ERROR: Unsupported 2D element type: " <<type;
          }
      }
    case SECOND:   {  // quadratic Lagrange shape functions   
      switch (type)      {
          case QUAD8:  case QUADSHELL8:      {
              const double xi  = p[0];   const double eta = p[1];     assert (i< 8);
              switch (i)    {
                case 0:  return .25*(1. - xi)*(1. - eta)*(-1. - xi - eta);
                case 1:  return .25*(1. + xi)*(1. - eta)*(-1. + xi - eta);
                case 2:  return .25*(1. + xi)*(1. + eta)*(-1. + xi + eta);
                case 3:   return .25*(1. - xi)*(1. + eta)*(-1. - xi + eta);
                case 4:  return .5*(1. - xi*xi)*(1. - eta);
                case 5:   return .5*(1. + xi)*(1. - eta*eta);
                case 6:   return .5*(1. - xi*xi)*(1. + eta);
                case 7:  return .5*(1. - xi)*(1. - eta*eta);
                default:   std::cerr<<"Invalid shape function index i = " << i;
                }
            }
          case QUAD4:
            assert(T == L2_LAGRANGE);std::cerr<< "High order on first order elements only supported for L2 families";
            (void(0));
          case QUAD9:   case QUADSHELL9: {
              // Compute quad shape functions as a tensor-product
              const double xi  = p[0];   const double eta = p[1];
              assert(i< 9);
              //                                0  1  2  3  4  5  6  7  8
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
              return (fe_lagrange_1D_quadratic_shape(i0[i], xi)*fe_lagrange_1D_quadratic_shape(i1[i], eta));
            }
          case TRI3:
            assert(T == L2_LAGRANGE);std::cerr<<  "High order on first order elements only supported for L2 families";
             (void(0));
          case TRI6:  case TRI7:  {
              const double zeta1 = p[0]; const double zeta2 = p[1];  const double zeta0 = 1. - zeta1 - zeta2;
              assert (i< 6);
              switch(i) {
                case 0: return 2.*zeta0*(zeta0-0.5);
                case 1: return 2.*zeta1*(zeta1-0.5);
                case 2: return 2.*zeta2*(zeta2-0.5);
                case 3: return 4.*zeta0*zeta1;
                case 4: return 4.*zeta1*zeta2;
                case 5: return 4.*zeta2*zeta0;
                default:  std::cerr<<"Invalid shape function index i = " << i;
                }
            }
          default:  std::cerr<<"ERROR: Unsupported 2D element type: " <<type;
          }
      }
      // "cubic" (one cubic bubble) Lagrange shape functions on TRI7
    case THIRD:  {
        switch (type)          {
          case TRI7:    {
              const double zeta1 = p[0];   const double zeta2 = p[1];
              const double zeta0 = 1. - zeta1 - zeta2;    const double bubble_27th = zeta0*zeta1*zeta2;
            assert (i< 7);
              switch(i)  {
                case 0:  return 2.*zeta0*(zeta0-0.5) + 3.*bubble_27th;
                case 1:  return 2.*zeta1*(zeta1-0.5) + 3.*bubble_27th;
                case 2:  return 2.*zeta2*(zeta2-0.5) + 3.*bubble_27th;
                case 3:  return 4.*zeta0*zeta1 - 12.*bubble_27th;
                case 4:  return 4.*zeta1*zeta2 - 12.*bubble_27th;
                case 5:  return 4.*zeta2*zeta0 - 12.*bubble_27th;
                case 6:  return 27.*bubble_27th;
                default:  std::cerr<<"Invalid shape function index i = " << i;
                }
            }
          default:  std::cerr<<"ERROR: Unsupported 2D element type: " <<type;
          }
      } // end case THIRD
    default:  std::cerr<<"ERROR: Unsupported 2D FE order: " << order; return 0;// unsupported order
    }
// #else // LIBMESH_DIM > 1
//   libmesh_ignore(type, order, i, p);
//   libmesh_not_implemented();
// #endif
}
// 
// 
// ===================================================================================== 
template <FEFamily T> double fe_lagrange_2D_shape_deriv(
  const ElemType type,
  const Order order,
  const unsigned int i,  // i-shape 
  const unsigned int j,  // d()/dx_j
  const double  p[]
){// ====================================================================================
// #if LIBMESH_DIM > 1
  assert (j< 2);
  switch (order)   { 
    case CONSTANT: {     return 0.;    }
    
    case FIRST: { // linear Lagrange shape functions // ------------------------------------------>
        switch (type)       {
          case QUAD4:  case QUADSHELL4:    case QUAD8:   case QUADSHELL8:   case QUAD9:   case QUADSHELL9:  {
              // Compute quad shape functions as a tensor-product
              const double xi  = p[0];  const double eta = p[1];   assert (i< 4);
              //                                0  1  2  3
              static const unsigned int i0[] = {0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1};
              switch (j) {
                case 0: // d()/dxi
                  return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)*fe_lagrange_1D_linear_shape(i1[i], eta));
                case 1: // d()/deta
                  return (fe_lagrange_1D_linear_shape(i0[i], xi)*fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta));
                default:  std::cerr<<"ERROR: Invalid derivative index j = " << j;abort(); return 0.;
                }
            }
          case TRI3:   case TRISHELL3:  case TRI6: case TRI7:{
              assert (i< 3);
              const double dzeta0dxi  = -1.;   const double dzeta1dxi  = 1.;  const double dzeta2dxi  = 0.;
              const double dzeta0deta = -1.;   const double dzeta1deta = 0.;    const double dzeta2deta = 1.;
              switch (j)  {
                case 0:  {   // d()/dxi          
                  switch(i)           {
                      case 0:     return dzeta0dxi;
                      case 1:     return dzeta1dxi;
                      case 2:     return dzeta2dxi;
                      default:   std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                case 1:   {  // d()/deta
                    switch(i)  {
                      case 0:  return dzeta0deta;
                      case 1:  return dzeta1deta;
                      case 2: return dzeta2deta;
                      default:   std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                default:   std::cerr<<"ERROR: Invalid derivative index j = " << j;abort(); return 0.;
                }
            }
          default:  std::cerr<<"ERROR: Unsupported 2D element type: " <<type; return 0.;
          }
      }
    
    case SECOND: { // quadratic Lagrange shape functions --------------------------------->
        switch (type)   {
          case QUAD8:  case QUADSHELL8:  {
              const double xi  = p[0];    const double eta = p[1];
              assert(i< 8);
              switch (j) {
                case 0: // d/dxi
                  switch (i) {
                    case 0:  return .25*(1. - eta)*((1. - xi)*(-1.) + (-1.)*(-1. - xi - eta));
                    case 1:  return .25*(1. - eta)*((1. + xi)*(1.) + (1.)*(-1. + xi - eta));
                    case 2:  return .25*(1. + eta)*((1. + xi)*(1.) +(1.)*(-1. + xi + eta));
                    case 3:  return .25*(1. + eta)*((1. - xi)*(-1.) + (-1.)*(-1. - xi + eta));
                    case 4: return .5*(-2.*xi)*(1. - eta);
                    case 5:  return .5*(1.)*(1. - eta*eta);
                    case 6:   return .5*(-2.*xi)*(1. + eta);
                    case 7:  return .5*(-1.)*(1. - eta*eta);
                    default:  std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                    }
                case 1:  // d/deta
                  switch (i)  {
                    case 0:  return .25*(1. - xi)*((1. - eta)*(-1.) + (-1.)*(-1. - xi - eta));
                    case 1:  return .25*(1. + xi)*((1. - eta)*(-1.) + (-1.)*(-1. + xi - eta));
                    case 2:  return .25*(1. + xi)*((1. + eta)*(1.) + (1.)*(-1. + xi + eta));
                    case 3:  return .25*(1. - xi)*((1. + eta)*(1.) + (1.)*(-1. - xi + eta));
                    case 4:  return .5*(1. - xi*xi)*(-1.);
                    case 5:  return .5*(1. + xi)*(-2.*eta);
                    case 6:  return .5*(1. - xi*xi)*(1.);
                    case 7:  return .5*(1. - xi)*(-2.*eta);
                    default:  std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                    }
                default:   std::cerr<<"ERROR: Invalid derivative index j = " << j;abort(); return 0.;
                }
            }
          case QUAD4: assert(T == L2_LAGRANGE);std::cerr<< "High order on first order elements only supported for L2 families";  (void(0));
          case QUAD9:  case QUADSHELL9:  {
              // Compute quad shape functions as a tensor-product
              const double xi  = p[0];   const double eta = p[1];   assert (i< 9);
              //                                0  1  2  3  4  5  6  7  8
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};

              switch (j) {
                case 0: // d()/dxi
                  return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*fe_lagrange_1D_quadratic_shape(i1[i], eta));
                case 1: // d()/deta
                  return (fe_lagrange_1D_quadratic_shape(i0[i], xi)*fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta));
                default: std::cerr<<"ERROR: Invalid derivative index j = " << j; abort();return 0.;
                }
            }
          case TRI3:  assert(T == L2_LAGRANGE);std::cerr<<"High order on first order elements only supported for L2 families"; (void(0));
          case TRI6:  case TRI7:  {
              assert (i< 6);
              const double zeta1 = p[0];  const double zeta2 = p[1]; const double zeta0 = 1. - zeta1 - zeta2;
              const double dzeta0dxi  = -1.;  const double dzeta1dxi  = 1.;  const double dzeta2dxi  = 0.;
              const double dzeta0deta = -1.; const double dzeta1deta = 0.; const double dzeta2deta = 1.;
              switch(j) {
                case 0: {
                    switch(i)  {
                      case 0: return (4.*zeta0-1.)*dzeta0dxi;
                      case 1: return (4.*zeta1-1.)*dzeta1dxi;
                      case 2: return (4.*zeta2-1.)*dzeta2dxi;
                      case 3: return 4.*zeta1*dzeta0dxi + 4.*zeta0*dzeta1dxi;
                      case 4: return 4.*zeta2*dzeta1dxi + 4.*zeta1*dzeta2dxi;
                      case 5: return 4.*zeta2*dzeta0dxi + 4*zeta0*dzeta2dxi;
                      default:  std::cerr<<"Invalid shape function index i = " << i; abort();return 0.;
                      }
                  }
                case 1: {
                    switch(i) {
                      case 0:  return (4.*zeta0-1.)*dzeta0deta;
                      case 1:  return (4.*zeta1-1.)*dzeta1deta;
                      case 2:  return (4.*zeta2-1.)*dzeta2deta;
                      case 3:  return 4.*zeta1*dzeta0deta + 4.*zeta0*dzeta1deta;
                      case 4:  return 4.*zeta2*dzeta1deta + 4.*zeta1*dzeta2deta;
                      case 5:  return 4.*zeta2*dzeta0deta + 4*zeta0*dzeta2deta;
                      default: std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                default:  std::cerr<<"ERROR: Invalid derivative index j = " << j;abort(); return 0.;
                }
            }
          default:  std::cerr<<"ERROR: Unsupported 2D element type: " <<type;abort(); return 0.;
          }
      }
     
    case THIRD:  { // "cubic" (one cubic bubble) Lagrange shape functions on TRI7 -------------------------->
        switch (type)   {
          case TRI7:  {
              assert (i< 7);
              const double zeta1 = p[0];   const double zeta2 = p[1];  const double zeta0 = 1. - zeta1 - zeta2;
              // const double bubble_27th = zeta0*zeta1*zeta2;
              const double dzeta0dxi  = -1.;   const double dzeta1dxi  = 1.;   const double dzeta2dxi  = 0.;
              const double dbubbledxi = zeta2 * (1. - 2.*zeta1 - zeta2);
              const double dzeta0deta = -1.;     const double dzeta1deta = 0.; const double dzeta2deta = 1.;
              const double dbubbledeta= zeta1 * (1. - zeta1 - 2.*zeta2);
              switch(j)  {
                case 0: {
                    switch(i)  {
                      case 0:  return (4.*zeta0-1.)*dzeta0dxi + 3.*dbubbledxi;
                      case 1:  return (4.*zeta1-1.)*dzeta1dxi + 3.*dbubbledxi;
                      case 2:  return (4.*zeta2-1.)*dzeta2dxi + 3.*dbubbledxi;
                      case 3:  return 4.*zeta1*dzeta0dxi + 4.*zeta0*dzeta1dxi - 12.*dbubbledxi;
                      case 4:  return 4.*zeta2*dzeta1dxi + 4.*zeta1*dzeta2dxi - 12.*dbubbledxi;
                      case 5:  return 4.*zeta2*dzeta0dxi + 4*zeta0*dzeta2dxi - 12.*dbubbledxi;
                      case 6:  return 27.*dbubbledxi;
                      default:  std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                case 1:  {
                    switch(i)   {
                      case 0:  return (4.*zeta0-1.)*dzeta0deta + 3.*dbubbledeta;
                      case 1:  return (4.*zeta1-1.)*dzeta1deta + 3.*dbubbledeta;
                      case 2:  return (4.*zeta2-1.)*dzeta2deta + 3.*dbubbledeta;
                      case 3:  return 4.*zeta1*dzeta0deta + 4.*zeta0*dzeta1deta - 12.*dbubbledeta;
                      case 4:  return 4.*zeta2*dzeta1deta + 4.*zeta1*dzeta2deta - 12.*dbubbledeta;
                      case 5:  return 4.*zeta2*dzeta0deta + 4*zeta0*dzeta2deta - 12.*dbubbledeta;
                      case 6:  return 27.*dbubbledeta;
                      default:  std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                default:   std::cerr<<"ERROR: Invalid derivative index j = " << j;abort(); return 0.;
                }
            }
          default:  std::cerr<<"ERROR: Unsupported 2D element type: " <<type; abort(); return 0.;// -------------------------->
          }
      } // end case THIRD
      // unsupported order
    default:  std::cerr<<"ERROR: Unsupported 2D FE order: " << order; abort(); return 0.;
    }
}
// 
// 
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
 template <FEFamily T> double fe_lagrange_2D_shape_second_deriv(
   const ElemType type,
   const Order order,
   const unsigned int i,
   const unsigned int j,
   const double p[])
{
// #if LIBMESH_DIM > 1
  // if(_dim<=1) {std::cerr << _dim ; std::abort();}
  // j = 0 ==> d^2 phi / dxi^2
  // j = 1 ==> d^2 phi / dxi deta
  // j = 2 ==> d^2 phi / deta^2
  assert (j< 3);
  switch (order)    {   // linear Lagrange shape functions
    case CONSTANT: {     return 0.;    }
    case FIRST:  {     switch (type)    {
          case QUAD4:   case QUADSHELL4:  case QUAD8:    case QUADSHELL8:     case QUAD9:    case QUADSHELL9:   {
              // Compute quad shape functions as a tensor-product
              const double xi  = p[0];              const double eta = p[1];  assert (i < 4);
              //                                0  1  2  3
              static const unsigned int i0[] = {0, 1, 1, 0};
              static const unsigned int i1[] = {0, 0, 1, 1};
              switch (j) {   // d^2() / dxi^2
                case 0:  return 0.;
                case 1: // d^2() / dxi deta
                  return (fe_lagrange_1D_linear_shape_deriv(i0[i], 0, xi)* fe_lagrange_1D_linear_shape_deriv(i1[i], 0, eta));
                case 2:  return 0.;// d^2() / deta^2
                default: std::cerr<<"ERROR: Invalid derivative index j = " << j;abort(); return 0.;
                }
            }
          case TRI3:  case TRISHELL3:  case TRI6: case TRI7:   {
              // All second derivatives for linear triangles are zero.
              return 0.;
            }
          default:  std::cerr<<"ERROR: Unsupported 2D element type: " ;return 0.;

          } // end switch (type)
      } // end case FIRST
      // quadratic Lagrange shape functions
    case SECOND:   {     switch (type)   {
          case QUAD8:    case QUADSHELL8:     {
              const double xi  = p[0];    const double eta = p[1];     assert (j< 3);
              switch (j)         {
                  // d^2z) / dxi^2
                case 0: {   switch (i)    {
                      case 0:   case 1:    return 0.5*(1.-eta);
                      case 2:   case 3:    return 0.5*(1.+eta);
                      case 4:   return eta - 1.;
                      case 5:   case 7:  return 0.0;
                      case 6:    return -1. - eta;
                      default:  std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                  // d^2() / dxi deta
                case 1: {    switch (i)  {
                      case 0: return 0.25*( 1. - 2.*xi - 2.*eta);
                      case 1:  return 0.25*(-1. - 2.*xi + 2.*eta);
                      case 2: return 0.25*( 1. + 2.*xi + 2.*eta);
                      case 3: return 0.25*(-1. + 2.*xi - 2.*eta);
                      case 4: return xi;
                      case 5: return -eta;
                      case 6: return -xi;
                      case 7: return eta;
                      default: std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }

                  // d^2() / deta^2
                case 2: {  switch (i) {
                      case 0:   case 3: return 0.5*(1.-xi);
                      case 1:   case 2: return 0.5*(1.+xi);
                      case 4:   case 6: return 0.0;
                      case 5:  return -1.0 - xi;
                      case 7:  return xi - 1.0;
                      default: std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                default:  std::cerr<<"ERROR: Invalid derivative index j = " << j;abort(); return 0.;
                } // end switch (j)
            } // end case QUAD8

          case QUAD4:
            assert(T == L2_LAGRANGE);std::cerr<<  "High order on first order elements only supported for L2 families"; (void(0));
          case QUAD9:       case QUADSHELL9:  {
              // Compute QUAD9 second derivatives as tensor product
              const double xi  = p[0]; const double eta = p[1];
              assert (i< 9);
              //                               0  1  2  3  4  5  6  7  8
              static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
              static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
              switch (j) {
                case 0:  // d^2() / dxi^2
                  return (fe_lagrange_1D_quadratic_shape_second_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape             (i1[i], eta));
                case 1:   // d^2() / dxi deta
                  return (fe_lagrange_1D_quadratic_shape_deriv(i0[i], 0, xi)*
                          fe_lagrange_1D_quadratic_shape_deriv(i1[i], 0, eta));
                case 2:   // d^2() / deta^2
                  return (fe_lagrange_1D_quadratic_shape             (i0[i], xi)*
                          fe_lagrange_1D_quadratic_shape_second_deriv(i1[i], 0, eta));

                default:   std::cerr<<"ERROR: Invalid derivative index j = " << j;
                }  // end switch (j)
            } // end case QUAD9

          case TRI3:  assert(T == L2_LAGRANGE);std::cerr<<  "High order on first order elements only supported for L2 families"; (void(0));
          case TRI6: case TRI7:      {
              const double dzeta0dxi  = -1.;   const double dzeta1dxi  = 1.;   const double dzeta2dxi  = 0.;
              const double dzeta0deta = -1.;    const double dzeta1deta = 0.;      const double dzeta2deta = 1.;   assert (j< 3);
              switch (j) {
                case 0:  {   switch (i)  {   // d^2() / dxi^2
                      case 0: return 4.*dzeta0dxi*dzeta0dxi;
                      case 1: return 4.*dzeta1dxi*dzeta1dxi;
                      case 2: return 4.*dzeta2dxi*dzeta2dxi;
                      case 3: return 8.*dzeta0dxi*dzeta1dxi;
                      case 4: return 8.*dzeta1dxi*dzeta2dxi;
                      case 5: return 8.*dzeta0dxi*dzeta2dxi;
                      default:  std::cerr<<"Invalid shape function index i = " << i;
                      }
                  }
                case 1:  { switch (i)  { // d^2() / dxi deta
                      case 0:return 4.*dzeta0dxi*dzeta0deta;
                      case 1:  return 4.*dzeta1dxi*dzeta1deta;
                      case 2:  return 4.*dzeta2dxi*dzeta2deta;
                      case 3:  return 4.*dzeta1deta*dzeta0dxi + 4.*dzeta0deta*dzeta1dxi;
                      case 4:  return 4.*dzeta2deta*dzeta1dxi + 4.*dzeta1deta*dzeta2dxi;
                      case 5:  return 4.*dzeta2deta*dzeta0dxi + 4.*dzeta0deta*dzeta2dxi;
                      default: std::cerr<<"Invalid shape function index i = " << i;
                      }
                  }
                case 2:                  {      switch (i)   {// d^2() / deta^2
                      case 0:  return 4.*dzeta0deta*dzeta0deta;
                      case 1: return 4.*dzeta1deta*dzeta1deta;
                      case 2: return 4.*dzeta2deta*dzeta2deta;
                      case 3: return 8.*dzeta0deta*dzeta1deta;
                      case 4: return 8.*dzeta1deta*dzeta2deta;
                      case 5: return 8.*dzeta0deta*dzeta2deta;
                      default:  std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                default:   std::cerr<<"ERROR: Invalid derivative index j = " << j;abort(); return 0.;
                } // end switch (j)
            }  // end case TRI6+TRI7
          default:  std::cerr<<"ERROR: Unsupported 2D element type: " ; return 0.;
          }
      } // end case SECOND
      // "cubic" (one cubic bubble) Lagrange shape functions on TRI7
    case THIRD:   {    switch (type)       {
          case TRI3:  assert(T == L2_LAGRANGE);std::cerr<<  "High order on first order elements only supported for L2 families"; (void(0));
          case TRI6:  case TRI7: {
              const double zeta1 = p[0];  const double zeta2 = p[1];
              // const double zeta0 = 1. - zeta1 - zeta2;
              const double dzeta0dxi  = -1.; const double dzeta1dxi  = 1.;   const double dzeta2dxi  = 0.;
              // const double dbubbledxi = zeta2 * (1. - 2.*zeta1 - zeta2);
              const double d2bubbledxi2 = -2. * zeta2;
              const double dzeta0deta = -1.; const double dzeta1deta = 0.;   const double dzeta2deta = 1.;
              // const double dbubbledeta= zeta1 * (1. - zeta1 - 2.*zeta2);
              const double d2bubbledeta2 = -2. * zeta1;
              const double d2bubbledxideta = (1. - 2.*zeta1 - 2.*zeta2);
              assert (j< 3);
              switch (j)  {
                case 0: { // d^2() / dxi^2
                    switch (i) {
                      case 0:  return 4.*dzeta0dxi*dzeta0dxi + 3.*d2bubbledxi2;
                      case 1:  return 4.*dzeta1dxi*dzeta1dxi + 3.*d2bubbledxi2;
                      case 2:  return 4.*dzeta2dxi*dzeta2dxi + 3.*d2bubbledxi2;
                      case 3:  return 8.*dzeta0dxi*dzeta1dxi - 12.*d2bubbledxi2;
                      case 4:  return 8.*dzeta1dxi*dzeta2dxi - 12.*d2bubbledxi2;
                      case 5:  return 8.*dzeta0dxi*dzeta2dxi - 12.*d2bubbledxi2;
                      case 6:  return 27.*d2bubbledxi2;
                      default:   std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                case 1: { // d^2() / dxi deta
                    switch (i)  {
                      case 0:   return 4.*dzeta0dxi*dzeta0deta + 3.*d2bubbledxideta;
                      case 1:   return 4.*dzeta1dxi*dzeta1deta + 3.*d2bubbledxideta;
                      case 2:   return 4.*dzeta2dxi*dzeta2deta + 3.*d2bubbledxideta;
                      case 3:   return 4.*dzeta1deta*dzeta0dxi + 4.*dzeta0deta*dzeta1dxi - 12.*d2bubbledxideta;
                      case 4:   return 4.*dzeta2deta*dzeta1dxi + 4.*dzeta1deta*dzeta2dxi - 12.*d2bubbledxideta;
                      case 5:   return 4.*dzeta2deta*dzeta0dxi + 4.*dzeta0deta*dzeta2dxi - 12.*d2bubbledxideta;
                      case 6:   return 27.*d2bubbledxideta;
                      default:  std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                case 2:  { // d^2() / deta^2
                    switch (i)  {
                      case 0:  return 4.*dzeta0deta*dzeta0deta + 3.*d2bubbledeta2;
                      case 1:  return 4.*dzeta1deta*dzeta1deta + 3.*d2bubbledeta2;
                      case 2:  return 4.*dzeta2deta*dzeta2deta + 3.*d2bubbledeta2;
                      case 3:  return 8.*dzeta0deta*dzeta1deta - 12.*d2bubbledeta2;
                      case 4:  return 8.*dzeta1deta*dzeta2deta - 12.*d2bubbledeta2;
                      case 5:  return 8.*dzeta0deta*dzeta2deta - 12.*d2bubbledeta2;
                      case 6:  return 27.*d2bubbledeta2;
                      default:  std::cerr<<"Invalid shape function index i = " << i;abort(); return 0.;
                      }
                  }
                default:   std::cerr<<"ERROR: Invalid derivative index j = " << j;abort(); return 0.;
                } // end switch (j)
            }  // end case TRI6+TRI7
          default:  std::cerr<<"ERROR: Unsupported 2D element type: " ; return 0.;
          }
      } // end case THIRD
      // unsupported order
    default:  std::cerr<<"ERROR: Unsupported 2D FE order: " << order; abort(); return 0.;
    } // end switch (order)
}
