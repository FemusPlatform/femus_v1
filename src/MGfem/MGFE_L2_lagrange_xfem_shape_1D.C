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

// std include ------------>
#include <cmath>
#include <array>
// // femus Local includes ----------->
#include "MGFE_L2_FE_template.h"
#include "MGFE_L0_C.h"
#include "MGFE_L2_lagrange_shape_1D.h"


// 
// LIBMESH_DEFAULT_VECTORIZED_FE(1,LAGRANGE_XFEM)
// LIBMESH_DEFAULT_VECTORIZED_FE(1,L2_LAGRANGE_XFEM)
// 
// // FE class members
template <>
 MGFE<1,LAGRANGE_XFEM>::MGFE (const FEType & fet) :
   MGFEbase (1,fet),
   // _dim(1),             ///< space dimension
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
 MGFE<1,LAGRANGE_XFEM>::~MGFE(){}
 
 
  template <>  const int MGFE<1,LAGRANGE_XFEM>::get_dim()       {return _dim;}  ///< Dimension (3 or 2) used in the 3D FEM (1< _dim < MDIM_FEM)
 template <>  const int MGFE<1,LAGRANGE_XFEM>::get_order()     {return _order;}///< Order of the shape functions 0=cpnst 1=linear 2=quad
 template <>  const int MGFE<1,LAGRANGE_XFEM>::get_NoShape()   {return _NoShape;}  ///< Number of shape functions of the fem
 template <>  const int MGFE<1,LAGRANGE_XFEM>::get_NoElNodes() {return _NoElNodes;}  ///< Number of nodes in the reference element (always ref quad)
 template <>  const int MGFE<1,LAGRANGE_XFEM>::get_GeomType()  {return _GeomType;}   ///< Geom type 1=HEXAHEDRAL 0=TETRAHEDRAL
 template <>  const int MGFE<1,LAGRANGE_XFEM>::get_FamType()   {return _FamType;}    ///< Family  (geom 1: Lagrangian =0  RT =1)
  
   // Gaussian points -----------------------------------------------------------------------------
 template <>  const int MGFE<1,LAGRANGE_XFEM>::get_NoGauss(){ return _NoGauss;}  ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
  template <> const double  MGFE<1,LAGRANGE_XFEM>::get_xg(int i){          return  _xg[i];}          ///< gaussian points
  template <> const double MGFE<1,LAGRANGE_XFEM>::get_weight(int i){      return _weight[i];}            ///< Weight
  template <> const double & MGFE<1,LAGRANGE_XFEM>::get_phi(){     return  *phi;}           ///< Shape functions
  template <> const double & MGFE<1,LAGRANGE_XFEM>::get_dphidxez_map(){return  *_dphidxez_map;}      ///< Shape derivative functions in gaussian points
  template <> const double & MGFE<1,LAGRANGE_XFEM>::get_dphidxx_map(){ return  *_dphidxx_map;}       ///< Second order shape derivatives in Gaussian points
   // Dofs points ----------------------------------------------------------------------------------
 template <>  const double &  MGFE<1,LAGRANGE_XFEM>::get_xnodes(){return   *_xnodes;}              ///< xnode (hex or tetra) in (xi,eta,seta)
 template <>  const int &  MGFE<1,LAGRANGE_XFEM>::get_fem_conn(){ return *_fem_conn;}              ///< dof node coords (respect to _xnodes)
  // template <> const int     MGFE<1,LAGRANGE_XFEM>::get_deg(){return _deg;}                  ///< Polinomial degree
  
   // Nodes (NDOF_FEM) -----------------------------------------------------------------------------
 template <>  const double & MGFE<1,LAGRANGE_XFEM>::get_phi_map_nodes(){return   *_phi_map_nodes;   }          ///< Shape functions
 template <>  const double & MGFE<1,LAGRANGE_XFEM>::get_dphidxez_map_nodes(){return   *_dphidxez_map_nodes;} 
 
 
 template <>    void MGFE<1,LAGRANGE_XFEM>::get_dphi_gl_g(const int qp, const double InvJac[], double dphi[]);
 template <>   void MGFE<1,LAGRANGE_XFEM>::get_ddphi_gl_g(const int qp, const double InvJac[], double ddphi[]);
template <>    void MGFE<1,LAGRANGE_XFEM>::get_dphi_gl_g(const int dim, const int qp, const double InvJac[], double dphi[], int sdim);
template <>    void MGFE<1,LAGRANGE_XFEM>::get_dphi_gl_g(const int dim, const int qp, const double InvJac[], std::vector<double>& dphi);
template <>    void MGFE<1,LAGRANGE_XFEM>::get_dphi_node(const int dim, const int qp, const double InvJac[], double dphi[]);

 template <>   void MGFE<1,LAGRANGE_XFEM>::get_dphi_on_given_node(const int dim, double ElemCoords[], double CanPos[], double dphi[]);
 template <>    void MGFE<1,LAGRANGE_XFEM>::get_dphi_on_given_nodeG(
    const int dim, double ElemCoords[], double CanPos[], double dphi[], int FamilyType){}
template <>    void MGFE<1,LAGRANGE_XFEM>::get_dphi_arb_node(std::vector<double> NodeCoord, const int order, double InvJac[], double dphi[]){}


  

 template <>   void MGFE<1,LAGRANGE_XFEM>::get_phi_gl_g(const int qp, double phi[]);
 template <>   void MGFE<1,LAGRANGE_XFEM>::get_phi_gl_g(const int qp, std::vector<double>& phi);
 template <>   void MGFE<1,LAGRANGE_XFEM>::get_phi_g_arb_el(const int dim, const int qp, double phi[], int FamilyType){}
  template <>  void MGFE<1,LAGRANGE_XFEM>::get_dphi_g_arb_el(const int dim, const int qp, const double InvJac[], double dphi[], int FamilyType){}
 
  template <>    void  MGFE<1,LAGRANGE_XFEM>::normal_g(double const*, double const*, double*) const{return;}
  template <>    void  MGFE<1,LAGRANGE_XFEM>::normal_g(double const*, double const*, double*, int &) const{return;}
  template <>    void MGFE<1,LAGRANGE_XFEM>::write(const std::string& filename){return ;}
 
 template <>    int MGFE<1,LAGRANGE_XFEM>::GetFamilyType(int elem_dof, int dim){return 1;}
 template <>    void MGFE<1,LAGRANGE_XFEM>::print_rec_xdmf(const std::string& filename,const std::string& namefem,const std::string& namelem){return ;}
 
 // ========================================================
 // Jacobian  ----------------------------------------------------
 // ========================================================
  template <> double MGFE<1,LAGRANGE_XFEM>::JacSur(const int ng,const   double x[],   double InvJac[]);
  template <> double MGFE<1,LAGRANGE_XFEM>::Jac(
     const int ng, 
     const double* x , 
     double* InvJac
  ){// ========================================================
      double x_xi = 0.;
      int nshape = _NoShape;
      for(int s = 0; s < nshape; s++)  x_xi += x[s] * _dphidxez_map[ s * _NoGauss + ng];
      double det=  x_xi;
      InvJac[0] = 1./(x_xi+1.e-20);  // deta dy

      return (det);
  }
  template <> double MGFE<1,LAGRANGE_XFEM>::JacG(const int,const double*, double*, int, int){return 0.;}
  template <> double MGFE<1,LAGRANGE_XFEM>::JacobianOnGauss(const int,const int,const double*, double*){return 0.;}
  template <> double MGFE<1,LAGRANGE_XFEM>::JacobianOnPoint(const int,const double*, double*, double*, int, int){return 0. ;}
 
  
  // quadrature -------------------------------------------------
 template <>
void MGFE<1,LAGRANGE_XFEM>::attach_quadrature_rule (QBase *q)
 {
assert(q);
   this->_qrule = q;
   // make sure we don't cache results from a previous quadrature rule
   this->_elem_type = INVALID_ELEM;
   return;
 }
 template <>
 unsigned int MGFE<1,LAGRANGE_XFEM>::n_shape_functions () const {
    return MGFE<1,LAGRANGE_XFEM>::n_dofs (this->_elem_type, this->_fe_type.order + this->_p_level);
  }
 
 
// ===========================================================================
//       Shape 
// ===========================================================================
template <> double MGFE<1,LAGRANGE_XFEM>::shape(
  const ElemType,
  const Order order,
  const unsigned int i,
  const double p[]
) {
   return fe_lagrange_1D_shape(order, i, p[0]);
 } // ===========================================================================
// 
// template <> Real  MGFE<1,L2_LAGRANGE_XFEM>::shape(const ElemType,   const Order order,const unsigned int i,const Point & p)
// {   return fe_lagrange_1D_shape(order, i, p(0)); }
//------------
// template <> Real  MGFE<1,LAGRANGE_XFEM>::shape(const Elem * elem, const Order order,  const unsigned int i,
//                            const Point & p,   const bool add_p_level) {
//   libmesh_assert(elem);   return fe_lagrange_1D_shape(order + add_p_level*elem->p_level(), i, p(0)); } 
// -------------------
// template <> Real  MGFE<1,LAGRANGE_XFEM>::shape(const FEType fet,  const Elem * elem,   const unsigned int i,
//                            const Point & p,   const bool add_p_level) {
//   libmesh_assert(elem);  return fe_lagrange_1D_shape(fet.order + add_p_level*elem->p_level(), i, p(0)); }
// ------------------
// template <> Real  MGFE<1,L2_LAGRANGE_XFEM>::shape(const Elem * elem,   const Order order,
//                               const unsigned int i,    const Point & p,   const bool add_p_level)
// { libmesh_assert(elem);  return fe_lagrange_1D_shape(order + add_p_level*elem->p_level(), i, p(0)); }
// -----------
// template <> Real  MGFE<1,L2_LAGRANGE_XFEM>::shape(const FEType fet,  const Elem * elem,   const unsigned int i,
//                               const Point & p,  const bool add_p_level)
// {  libmesh_assert(elem); return fe_lagrange_1D_shape(fet.order + add_p_level*elem->p_level(), i, p(0)); }
// 
// 
// ===========================================================================
//       Shape derivative (1st)
// ===========================================================================
// ==================================================================
 template <>
 double MGFE<1,LAGRANGE_XFEM>::shape_deriv(
   const ElemType, 
   const Order order, 
   const unsigned int i,
   const unsigned int j,
   const double  p[]
 ){ return fe_lagrange_1D_shape_deriv(order, i, j, p[0]);}

// template <> Real  MGFE<1,L2_LAGRANGE_XFEM>::shape_deriv(const ElemType,  const Order order,
//                                     const unsigned int i,  const unsigned int j,  const Point & p)
// {   return fe_lagrange_1D_shape_deriv(order, i, j, p(0)); }
// 
// 
// 
// template <>  Real  MGFE<1,LAGRANGE_XFEM>::shape_deriv(const Elem * elem,  const Order order,    const unsigned int i,
//                                  const unsigned int j, const Point & p,   const bool add_p_level)
// { libmesh_assert(elem);   return fe_lagrange_1D_shape_deriv(order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// template <> Real  MGFE<1,L2_LAGRANGE_XFEM>::shape_deriv(const Elem * elem, const Order order,  const unsigned int i,
//                                     const unsigned int j,  const Point & p, const bool add_p_level)
// {  libmesh_assert(elem);  return fe_lagrange_1D_shape_deriv(order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// template <> Real  MGFE<1,LAGRANGE_XFEM>::shape_deriv(const FEType fet,  const Elem * elem,   const unsigned int i,
//                                  const unsigned int j,  const Point & p,  const bool add_p_level) {
//   libmesh_assert(elem);  return fe_lagrange_1D_shape_deriv(fet.order + add_p_level*elem->p_level(), i, j, p(0)); }

// template <> Real  MGFE<1,L2_LAGRANGE_XFEM>::shape_deriv(const FEType fet, const Elem * elem,
//                                     const unsigned int i, const unsigned int j, const Point & p,
//                                     const bool add_p_level) {  libmesh_assert(elem);
// return fe_lagrange_1D_shape_deriv(fet.order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// ===========================================================================
//       Shape derivative (2nd)
// ===========================================================================
// ================================================================================================ 
template <>
double  MGFE<1,LAGRANGE_XFEM>::shape_second_deriv(
  const ElemType,
  const Order order,
  const unsigned int i,
  const unsigned int j,
  const double  p[]
) {// ================================================================================================
  return fe_lagrange_1D_shape_second_deriv(order, i, j, p[0]);
 }
//  
// template <> Real  MGFE<1,L2_LAGRANGE_XFEM>::shape_second_deriv(const ElemType, const Order order,   const unsigned int i,
//                                            const unsigned int j, const Point & p){
//   return fe_lagrange_1D_shape_second_deriv(order, i, j, p(0)); }

// template <>
// Real  MGFE<1,LAGRANGE_XFEM>::shape_second_deriv(const Elem * elem,    const Order order,
//                                         const unsigned int i,   const unsigned int j,
//                                         const Point & p,  const bool add_p_level)
// {   libmesh_assert(elem);  return fe_lagrange_1D_shape_second_deriv(order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// template <> Real  MGFE<1,L2_LAGRANGE_XFEM>::shape_second_deriv(const Elem * elem,     const Order order,
//                                            const unsigned int i,    const unsigned int j, const Point & p,
//                                            const bool add_p_level)
// {   libmesh_assert(elem);  return fe_lagrange_1D_shape_second_deriv(order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// template <> Real  MGFE<1,LAGRANGE_XFEM>::shape_second_deriv(const FEType fet,   const Elem * elem,      const unsigned int i,
//                                         const unsigned int j,           const Point & p,
//                                         const bool add_p_level) {
//   libmesh_assert(elem);   return fe_lagrange_1D_shape_second_deriv(fet.order + add_p_level*elem->p_level(), i, j, p(0));
// }

// template <>
// Real  MGFE<1,L2_LAGRANGE_XFEM>::shape_second_deriv(const FEType fet,     const Elem * elem,  const unsigned int i,
//                                            const unsigned int j,   const Point & p,  const bool add_p_level)
// {   libmesh_assert(elem);  return fe_lagrange_1D_shape_second_deriv(fet.order + add_p_level*elem->p_level(), i, j, p(0)); }
// 
// #endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
// } // namespace libMesh
// 

template <> double   MGFE<1,LAGRANGE_XFEM>:: compute_fe_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz[], // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
    double InvJac[]
 ){
  
  return 0.;
}
template <> double   MGFE<1,LAGRANGE_XFEM>:: compute_fe_face_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz[] // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
 ){
  
  return 0.;
}

template <>   double MGFE<1,LAGRANGE_XFEM>::JacSur(
  const int ng, 
  const double x[], 
  double InvJac[]
)   {// ========================================================
       double det=0.; double x_xi = 0.;
      int nshape = _NoShape;
       for(int  kdim = 0;  kdim < 2;  kdim++){ // FEM1dim ->JacSurf 2D -> 2
        x_xi = 0.;
      for(int s = 0; s < nshape; s++)  x_xi += x[s+ kdim *_NoGauss] * _dphidxez_map[ s * _NoGauss +   ng];
       det +=  x_xi* x_xi;
       }
      det= std::sqrt(det);
      InvJac[0] = 1./(det+1.e-20);  // deta dy

      return (det);
  }
  
  
  
// =======================================================
/// This function computes the  derivatives and Jacobian at the gaussian point ng
template <>   double MGFE<1,LAGRANGE_XFEM>::Jac_nodes(
const int ng,          // Gaussian point <-
 const double x[],            // coordinates  <-
double InvJac[]        // Jacobean
) {                        // ===================================================
  double x_xi = 0.;
  int nshape = _NoShape;
  //   int offset=_NoShape*_NoGauss;
  for(int s = 0; s < nshape; s++) {
    int sng = s * _NoGauss + ng;
    x_xi += x[s] * _dphidxez_map_nodes[sng];
  }
  double det = x_xi;
  double idet = 1. / det;
  InvJac[0] = x_xi * idet;  // deta dy
  return (det);
}





// ==============================================================
template<> double  MGFE<1,LAGRANGE_XFEM>::JacOnGivenCanCoords(
  const int dim, 
  const double ElemCoords[], 
  double CanCoords[], 
  double InvJac[], 
  int FamilyType, 
  int nShape
) {// ===========================================================
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
//   double det = CoordsDer[0];
// 
//  
// 
//   delete[] LocDPhi;
//   delete[] CoordsDer;
  // return det; 
  return 0.;
  
}

















// ======================================
/// This function computes the normal at the gauss point
template <>
void MGFE<1,LAGRANGE_XFEM>::normal_g(
const double* xx,  // all surface coordinates <-
double* normal_g   // normal ->
) const {          // ======================================
  int NDOF_FEMB=1;
   // coordinates
  double xx3D[3 * NDOF_FEMB];
  normal_g[0]=1.;
 

  return;
}

// =========================================================================================
template <>
void  MGFE<1,LAGRANGE_XFEM>::Oxy_face(
const double* xx,   ///< coordinates
const double x_c[], ///< central point
double* normal_g,   ///<  i_n,i_t1,i_t2 face Oxyz
double* tg01,
double* tg03
) const {            // ======================================
int NDOF_FEMB=1;
  // coordinates
  double xx3D[3 * NDOF_FEMB];
  normal_g[0]=1.;
  for(int i = 0; i < 1 * NDOF_FEMB; i++) { xx3D[i] = xx[i]; tg01[i] =0.;tg01[3] =0.;}
  // if the sign is not correct then reverse it
  if(normal_g[0] * (x_c[0] - xx3D[0])  >    0) {    normal_g[0] *= -1;    }

  return;
}

// =============================================================================
template <> void MGFE<1,LAGRANGE_XFEM>::set_par(
  std::string fileh5,
  int dim_in,
  int order_in,
  int shapes_quad,
  int geom_type_in,
  int fam_type_in
) {// =============================================================================
  _dim=dim_in, _order= order_in,  
  _NoShape=shapes_quad,   _GeomType=geom_type_in,     _FamType=fam_type_in;
   femfileh5_=fileh5;
     ///< family type  0=Lagrangian(C^0) 1=Raviart Thomas (L^2-C^0)
//   // Fem order
//   // =============================
//   //  ------------------------------------------------------------------------
//   // gauss points (_NoGauss) ------------------------------------------------------
//     _NoGauss=0;
//   assert(_qrule);
//  _NoGauss=_qrule->n_points();
//  if(_NoGauss==0) {
//    std::cout << "Quadrature not assigned"; abort();
// }
//    // _NoGauss=3; 
//   ElemType  elem_name=_fe_type.fetype; 
//   
//   // number of nodes (_NoElNodes) in the reference EDGE3 element
//   _NoElNodes= get_refspace_n_nodes(elem_name);// number of nodes on element  ElemType (it should be from Geom) todo3;  
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
    _NoElNodes= get_refspace_n_nodes(geom_elem_name);// mesh nodes
  _NoShape=get_refspace_n_nodes(elem_name); // element nodes 
 
  //  -----------------------------------------------------------------------------------
  // x node  coordinates   (_xnodes[xp1,xp2,xp3,...yp1,yp2,yp3,....zp1,zp2,zp3,...] )----
  _xnodes=new double[ _NoElNodes];
  // shape functions at nodes (_phi_map_nodes[sh1(x1),sh1(x2),....  sh2(x1),sh2(x2),.... ])
  _phi_map_nodes = new double[_NoShape* _NoElNodes];
  // 1st derivatives at nodes x _dphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxez_map_nodes = new double[ _NoShape * _NoElNodes];
  // element dofs  (reference to nodes) --------------------------------------------------
  _fem_conn=new int[_NoShape];
  // gaussian (g) coordinates (_xg[xg1,xg2,xg3,...yg1,yg2,yg3,....zg1,zg2,zg3,...] )-----
  _xg=new double[_NoGauss];
  // weights [wg1,wg2,wg3,.......]
  _weight = new double[_NoGauss];
  // shape functions _phi[sh1(xg1),sh1(xg2),....  sh2(xg1),sh2(xg2),.... ]
  // 1st derivatives at g _dphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxez_map = new double[ _NoShape* _NoGauss];
  // 2nd derivatives at g _ddphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxx_map = new double[_NoShape * _NoGauss];
  //   -------------------------------------------------------------------------------------
  
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
      _phi_map_nodes[j*_NoElNodes+i]=fe_lagrange_1D_shape( elem_order, j, point[0]);//Rec_Quad_Phi(j, point,_dim);

      // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
      for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
        _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] = fe_lagrange_1D_shape_deriv( elem_order, j,dir, point[0]);// Rec_Quad_DPhi(j,point,_dim,dir);     // dphi/dxi
      }
    }
  }
   delete[]xnodes;
// gaussian points --------------------------------------------------------->
 for(int i = 0; i < _NoGauss; i++) {
    _weight[i] =_qrule->w(i); // weight
     _qrule->qp(point,i);     // ith-gaussian point ->point(i)
    for(int j = 0; j < _NoShape; j++) {     // loop over test function id
      phi[j*_NoGauss+i] =fe_lagrange_1D_shape( elem_order, j, point[0]);
      for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
        _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_deriv( elem_order, j,dir, point[0]);
        for(int dir2 = 0; dir2 < _dim; dir2++)// dir +dir2=0->d^2phi/dxi^2; dir +dir2=2->d^2phi/deta^2; dir +dir2=1->d^2phi/deta dxi;
          _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =fe_lagrange_1D_shape_second_deriv(elem_order,j,dir +dir2, point[0]);
      }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
  } 
  delete []point;
 
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
  std::cout << "\n MGFE: elem_name(Dim): "           <<elem_name<<"("<< _dim <<")" ;
  std::cout << "\n       gaussian points: "          << _NoGauss ;
  std::cout << "\n       number of shape functions " << _NoShape ;
  std::cout << "\n       polynomial order "          << elem_order;
  std::cout << " \n";
#endif
  return;
}



// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// // Basic Raviart Thomas non continuous functions
// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// //===============================================================================================
// template <>
// std::array<double,3 >  MGFE<1,LAGRANGE_XFEM>::Rec_RT_Phi(int iface, double x[],int dim) {
// //===============================================================================================
//   switch(iface) {
//   case 0:
//     return {1.,0.,0.};    break;
//   case 1:
//     return {1.,0.,0.};    break;
//   default:
//     std::cout<< "MGFE1::RT_quad Error iface "; return  {0.,0.,0.};  exit(3);
//   }
//   return  {0.,0.,0.};
// }
// //===============================================================================================
// 
// 
// template <>
// std::array<double,3>  MGFE<1,LAGRANGE_XFEM>::Rec_RT_DPhi(int /*iface*/, double x[],int dim, int /*dir*/) {
//   return  {0.,0.,0.};  
// }



// OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// Basic special operators
// OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO



template <> double MGFE<1,LAGRANGE_XFEM>::ComputeInverseMatrix_vol(double Matrix[], double InvMatrix[]) {
  double det = Matrix[0] ;
  double idet = 1. / det;
  InvMatrix[0] =  idet;   // dxi dx
  return det;
}
template <> double MGFE<1,LAGRANGE_XFEM>::ComputeInverseMatrix_sur(double Matrix[], double InvMatrix[]) {
  double det = 1.;
  InvMatrix[0] = 1;   // dxi dx
  return det;
}


//===============================================================================================
template <> double MGFE<1,LAGRANGE_XFEM>::ComputeInverseMatrix(
  double Matrix[], 
  double InvMatrix[], 
  int dim) {
  double det = 0.;

  if(dim == 1)   det = ComputeInverseMatrix_vol(Matrix,InvMatrix);
  if(dim == 0) det = ComputeInverseMatrix_sur(Matrix,InvMatrix);
  
  return det;
}
//===============================================================================================







//===============================================================================================
template <>  double MGFE<1,LAGRANGE_XFEM>::Tj_interp_lin(double *u,int j, double *xi,int ndof,int idim) {
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
template <>  double MGFE<1,LAGRANGE_XFEM>::dTjdxi_lin(double *u,int j, double *xi,int i,int ndof,int idim) {
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
//========================================================================================
// =========================================================
//  GET phi dphi ddphi
//==========================================================
//========================================================================================
template <> void MGFE<1,LAGRANGE_XFEM>::get_dphi_on_given_node(
  const int dim, 
  double ElemCoords[], 
  double CanPos[], 
  double dphi[]
) {
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
}//======================================================================



// ==================================================================
/// This function computes the shape values at the gauss point qp
template <> void MGFE<1,LAGRANGE_XFEM>::get_phi_gl_g(
const int qp,    // gaussian point <-
double phi[]     // shape functions ->
) {                  // =================================================
  for(int ish = 0; ish < _NoShape; ish++) phi[ish] = phi[ish * _NoGauss + qp];
  return;
}

// ==================================================================
/// This function computes the shape values at the gauss point qp
template <> void MGFE<1,LAGRANGE_XFEM>::get_phi_gl_g(
const int qp,             // gaussian point   <-
std::vector<double>& phi  // shape functions  ->
) { // ==============================================================
  for(int ish = 0; ish < _NoShape; ish++) phi[ish] = phi[ish * _NoGauss + qp];

  return;
}

// // ==================================================================
// void MGFE1::get_phi_gl_g(
//   double point[],           // gaussian point <
//   double phi[],
//   int fam_type,
//   int sdim
// ) { // ==============================================================
// 
//   const int el_nnodes = _NoShape;   // # of shape functions
// 
//   for(int eln = 0; eln < el_nnodes; eln++) {
// 
//     phi[eln] = LocalPhi(eln, point, sdim, fam_type) ;
// 
//   }
// 
//   return;
// }
// void MGFE1::get_phi_g_arb_el(
// const int kdim,    // dimension <-
// const int qp,      // gaussian point <-
// double phi[],      // shape functions ->
// int FamilyType) {  // =================================================
//   for(int ish = 0; ish < _GNoShape[FamilyType][kdim - 1]; ish++) {
//     phi[ish] = _Gphi_map1[(kdim - 1) + FamilyType * 3][ish * _GNoGauss1[FamilyType][kdim - 1] + qp];
//   }
//   return;
// }
// =================================================================
/// Shape functions derivatives dphi[id+i* el_nnodes] =Ti
/// Tx Ty Tz
///  tensor  order  at the gauss point qp, node id:
///  we have  dphi[Tx(0),Tx(1),..Tx(id), Ty(0),Ty(1),..Ty(id), Tz(0),Tz(1),..Tz(id) ]
template <> void MGFE<1,LAGRANGE_XFEM>::get_dphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double dphi[]           // global derivatives ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  std::vector<double> gradphi_g(_dim);  // temp grad phi

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

// void MGFE1::get_dphi_g_arb_el(
// const int kdim,         // dimension <-
// const int qp,           // gaussian point <
// const double InvJac[],  // Jacobean
// double dphi[],          // global derivatives ->
// int FamilyType) {       // =========================================
//
//   const int el_nnodes = _GNoShape[FamilyType][kdim - 1];             // # of shape functions
//   const int el_ngauss = _GNoGauss1[FamilyType][kdim - 1];            // # of gauss points
//   const int goffset = el_nnodes * _GNoGauss1[FamilyType][kdim - 1];  // gauss offset
//
//   double gradphi_g[DIMENSION];  // temp grad phi
//
//   for(int eln = 0; eln < el_nnodes; eln++) {
//     int lqp = eln * el_ngauss + qp;
//     for(int idim = 0; idim < _dim; idim++) {
//       gradphi_g[idim] = _Gdphidxez_map1[kdim - 1 + 3 * FamilyType][lqp + idim * goffset];
//     }
//     for(int idim = 0; idim < _dim; idim++) {
//       double sum = 0.;
//       for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * gradphi_g[jdim]; }
//       dphi[eln + idim * el_nnodes] = sum;
//     }
//   }
//
//   return;
// }
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
template <> void MGFE<1,LAGRANGE_XFEM>::get_ddphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double ddphi[]          // global derivatives ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset
  // double ddphi_loc[NDOF_FEM * DIMENSION * DIMENSION];
   double ddphi_loc[3 * 1 * 1];

  for(int eln = 0; eln < el_nnodes; eln++) {   // LOOP OVER NODES ===================================
    const int shift = eln * _dim * _dim;       // offset for node derivatives
    const int lqp = eln * el_ngauss + qp;

    for(int idim = 0; idim < _dim * _dim; idim++) {
      ddphi_loc[shift + idim] = _dphidxx_map[lqp + idim * goffset];
    }

    // double Hess_tmp[DIMENSION * DIMENSION], JacTf[DIMENSION * DIMENSION];
      double Hess_tmp[1 * 1], JacTf[1 * 1];

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
template <> void MGFE<1,LAGRANGE_XFEM>::get_dphi_node(
const int /*kdim*/,         // dimension <-
const int node,         // nodal point <-
const double InvJac[],  // Jacobean <-
double dphi[]           // global derivatives at the nodal point ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  // double gradphi_g[DIMENSION];  // temp grad phi
  double gradphi_g[1];  // temp grad phi

  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + node;

    for(int idim = 0; idim < _dim; idim++) {
      gradphi_g[idim] = _dphidxez_map_nodes[lqp + idim * goffset];
    }

    for(int idim = 0; idim < _dim; idim++) {
      double sum = 0.;

      for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * gradphi_g[jdim]; }

      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

// // =============================================================
// /// Shape functions derivatives at the nodal points
// void MGFE1::get_dphi_arb_node(
// std::vector<double> NodeCoord,
// const int order,
// double InvJac[],
// double dphi[]) {  // =========================================
//   const int kdim = NodeCoord.size();
//
//   const int el_nnodes = _NoShape;             // # of shape functions
//   const int el_ngauss = _NoGauss1[kdim - 1];            // # of gauss points
//   const int goffset = el_nnodes * _NoGauss1[kdim - 1];  // gauss offset
//
//   double gradphi_g[DIMENSION];  // temp grad phi
//
//   int el_nodes = (order == 1) ? 4 : 9;
//   if(kdim == 3) { el_nodes *= (order == 1) ? 2 : 3; }
//
//   for(int eln = 0; eln < el_nodes; eln++) {
//     gradphi_g[0] = 1.;
//     gradphi_g[1] = 1.;
//     gradphi_g[kdim - 1] = 1.;
//     for(int idim = 0; idim < kdim; idim++) {   // derivative of eln-th phi ind idim direction
//       const int idim2 = (idim + 1) % kdim;
//       gradphi_g[idim] =
//       0.5 * _CooQ9[eln + idim * _Q9Off] * (0.5 * (1. + NodeCoord[idim2] * _CooQ9[eln + idim2 * _Q9Off]));
//     }
//     for(int idim = 0; idim < _dim; idim++) {
//       double sum = 0.;
//       for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * gradphi_g[jdim]; }
//       dphi[eln + idim * el_nodes] = sum;
//     }
//   }
//   return;
// }

// =============================================================
/// Shape functions derivatives at the gauss point qp
template <> void MGFE<1,LAGRANGE_XFEM>::get_dphi_gl_g(
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
      double sum = 0.;
      for(int jdim = 0; jdim < sdim; jdim++) { sum += InvJac[jdim + idim * sdim] * gradphi_g[jdim]; }
      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

// // =============================================================
// /// Shape functions derivatives at the gauss point qp
// template <> void MGFE<1,LAGRANGE_XFEM>::get_dphi_gl_g(
// const int /*kdim*/,         // dimension <-
// double point[],           // gaussian point <
// const double InvJac[],  // Jacobean
// double dphi[],          // global derivatives ->
// int sdim,
// int fam_type) {             // =========================================
// 
//   const int el_nnodes = _NoShape;   // # of shape functions
//   const int el_ngauss = _NoGauss;  // # of guass points
//   const int goffset = el_nnodes * el_ngauss;  // gauss offset
// 
//   double gradphi_g[3];  // temp grad phi
// 
//   for(int idim = 0; idim < 3; idim++) { gradphi_g[idim] = 0.; }
// 
//   for(int eln = 0; eln < el_nnodes; eln++) {
// 
//     for(int idim = 0; idim < sdim; idim++) { gradphi_g[idim] = FirstDerivateOfLocalPhi(eln, point, sdim, idim, fam_type); }
// 
//     for(int idim = 0; idim < sdim; idim++) {
//       double sum = 0.;
// 
//       for(int jdim = 0; jdim < sdim; jdim++) { sum += InvJac[jdim + idim * sdim] * gradphi_g[jdim]; }
// 
//       dphi[eln + idim * el_nnodes] = sum;
//     }
//   }
// 
//   return;
// }

// =============================================================
/// Shape functions derivatives at the gauss point qp
template <> void MGFE<1,LAGRANGE_XFEM>::get_dphi_gl_g(
const int /*kdim*/,            // dimension <-
const int qp,              // gaussian point <-
const double InvJac[],     // Jacobean
std::vector<double>& dphi  // global derivatives ->
) {                            // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of guass points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  double dphidxi_g[1];  // temp grad phi

  for(int eln = 0; eln < el_nnodes; eln++) {
    int lqp = eln * el_ngauss + qp;

    for(int idim = 0; idim < _dim; idim++) { dphidxi_g[idim] = _dphidxez_map[lqp + idim * goffset]; }

    for(int idim = 0; idim < _dim; idim++) {
      double sum = 0.;

      for(int jdim = 0; jdim < _dim; jdim++) { sum += InvJac[jdim + idim * _dim] * dphidxi_g[jdim]; }

      dphi[eln + idim * el_nnodes] = sum;
    }
  }

  return;
}

template <>  bool MGFE<1,LAGRANGE_XFEM>::phys_to_ref(double xref[],double xphys[],const double Xel[], int max_it, double tol) {return false;}


