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
// 
// // Local includes
 #include "MGFE_L2_FE_template.h"
 #include "MGFE_L0_C.h"
// #include "libmesh/elem.h"
// 
// 
#include<array>
#include<cmath>
// namespace libMesh
// {
// 
// 
// LIBMESH_DEFAULT_VECTORIZED_FE(0,LAGRANGE)
// LIBMESH_DEFAULT_VECTORIZED_FE(0,L2_LAGRANGE)
// 
template <>
 MGFE<0,LAGRANGE>::MGFE (const FEType & fet) :
   MGFEbase (0,fet),
   // _dim(0),             ///< space dimension
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
   assert (LAGRANGE == this->get_family());
 }
template <>
 MGFE<0,LAGRANGE>::~MGFE(){}

 template <>
 unsigned int MGFE<0,LAGRANGE>::n_shape_functions () const {
    return MGFE<0,LAGRANGE>::n_dofs (this->_elem_type, this->_fe_type.order + this->_p_level);
  }
  
   template <>  const int MGFE<0,LAGRANGE>::get_dim()       {return _dim;}  ///< Dimension (3 or 2) used in the 3D FEM (1< _dim < MDIM_FEM)
 template <>  const int MGFE<0,LAGRANGE>::get_order()     {return _order;}///< Order of the shape functions 0=cpnst 1=linear 2=quad
 template <>  const int MGFE<0,LAGRANGE>::get_NoShape()   {return _NoShape;}  ///< Number of shape functions of the fem
 template <>  const int MGFE<0,LAGRANGE>::get_NoElNodes() {return _NoElNodes;}  ///< Number of nodes in the reference element (always ref quad)
 template <>  const int MGFE<0,LAGRANGE>::get_GeomType()  {return _GeomType;}   ///< Geom type 1=HEXAHEDRAL 0=TETRAHEDRAL
 template <>  const int MGFE<0,LAGRANGE>::get_FamType()   {return _FamType;}    ///< Family  (geom 1: Lagrangian =0  RT =1)
  
   // Gaussian points -----------------------------------------------------------------------------
 template <>  const int MGFE<0,LAGRANGE>::get_NoGauss(){ return _NoGauss;}  ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
  template <> const double  MGFE<0,LAGRANGE>::get_xg(int i){          return  _xg[i];}          ///< gaussian points
  template <> const double MGFE<0,LAGRANGE>::get_weight(int i){      return _weight[i];}            ///< Weight
  template <> const double & MGFE<0,LAGRANGE>::get_phi(){     return  *phi;}           ///< Shape functions
  template <> const double & MGFE<0,LAGRANGE>::get_dphidxez_map(){return  *_dphidxez_map;}      ///< Shape derivative functions in gaussian points
  template <> const double & MGFE<0,LAGRANGE>::get_dphidxx_map(){ return  *_dphidxx_map;}       ///< Second order shape derivatives in Gaussian points
   template <>    void MGFE<0,LAGRANGE>::get_dphi_on_given_nodeG(int, double*, double*, double*, int){return ;}
    template <>    void MGFE<0,LAGRANGE>::get_dphi_g_arb_el(int, int, double const*, double*, int){return ;}
   // Nodes (NDOF_FEM) -----------------------------------------------------------------------------
 template <>  const double &  MGFE<0,LAGRANGE>::get_xnodes(){return   *_xnodes;}              ///< xnode (hex or tetra) in (xi,eta,seta)
 template <>  const int &  MGFE<0,LAGRANGE>::get_fem_conn(){ return *_fem_conn;}              ///< dof node coords (respect to _xnodes)
  // template <> const int     MGFE<0,LAGRANGE>::get_deg(){return _deg;}                  ///< Polinomial degree
 template <>  const double & MGFE<0,LAGRANGE>::get_phi_map_nodes(){return   *_phi_map_nodes;   }          ///< Shape functions
 template <>  const double & MGFE<0,LAGRANGE>::get_dphidxez_map_nodes(){return   *_dphidxez_map_nodes;} 
     template <>    void MGFE<0,LAGRANGE>::get_dphi_arb_node(std::vector<double> NodeCoord, const int order, double InvJac[], double dphi[]){return ;}
 
 
  // quadrature rule -------------------------------------------------
template <>
void MGFE<0,LAGRANGE>::attach_quadrature_rule (
  QBase *q){
assert(q);
   this->_qrule = q;
   // make sure we don't cache results from a previous quadrature rule
   this->_elem_type = INVALID_ELEM;
   return;
 }
 
  
   template <>    void MGFE<0,LAGRANGE>::print_rec_xdmf(const std::string& filename,const std::string& namefem,const std::string& namelem){return ;}
   template <>    void MGFE<0,LAGRANGE>::write(const std::string& filename){return ;}
   // template <>    void MGFE<0,LAGRANGE>::init_pie_tri(){return ;}
  
   template <>    int MGFE<0,LAGRANGE>::GetFamilyType(int elem_dof, int dim){return 1;}
  
  
 
 
  

   
      template <>    double MGFE<0,LAGRANGE>::ComputeInverseMatrix_vol(double*, double*){return 0.;}
 // ========================================================
 // Jacobian  ----------------------------------------------------
 // ========================================================  
   template<> double MGFE<0,LAGRANGE>::Jac   (const int ng,const double*, double*){return 0.;}
   template<> double MGFE<0,LAGRANGE>::JacSur(const int ng,const double x[], double InvJac[]) {return 0.;} 
   template<> double MGFE<0,LAGRANGE>::JacG  (const int,const double*, double*, int, int){return 0.;}
   template<> double MGFE<0,LAGRANGE>::Jac_nodes          (const int ng,const double x[], double InvJac[])  {return 0.;}
   template<> double MGFE<0,LAGRANGE>::JacobianOnGauss    (const int ng,const int,const double*, double*){return 0.;}
   template<> double MGFE<0,LAGRANGE>::JacobianOnPoint    (const int ng,const double*, double*, double*, int, int){return 0. ;}
   template<> double MGFE<0,LAGRANGE>::JacOnGivenCanCoords(const int dim,const double ElemCoords[], double CanCoords[], double InvJac[], int FamilyType, int nShape) {return 0.;}
 // ========================================================
 // Normal  ----------------------------------------------------
 // ======================================================== 
  template <>    void  MGFE<0,LAGRANGE>::normal_g(double const*, double const*, double*, int&) const{return;}
  template <>    void  MGFE<0,LAGRANGE>::normal_g(double const*, double const*, double*) const{return;}
   // template <>    double & MGFE<0,LAGRANGE>::normal_g(double const*, double const*, double*) const{return ;}
 // ========================================================
 //  RETURN SHAPE FUNCTIONS AT GAUSS POINTS --------------------------------------------
    ///< \param[in]  <dim>    Dimension
    ///< \param[in]  <qp>     Gaussian point
    ///< \param[out] <phi>    Shape function
  // ========================================================  

  template <> void MGFE<0,LAGRANGE>::get_phi_g_arb_el(const int dim, const int qp, double phi[], int FamilyType){return ;}
  
  template <>  bool MGFE<0,LAGRANGE>::phys_to_ref(double xref[],double xphys[],const double Xel[], int max_it, double tol) {return false;}
  
//  LAGRANGE
// template <>
//  Real MGFE<0,LAGRANGE>::shape(const Elem *,  const Order,  const unsigned int libmesh_dbg_var(i),   const Point &,   const bool)
//  { assert (i< 1);   return 1.;
// }
// template <>
// Real MGFE<0,LAGRANGE>::shape(const FEType,  const Elem *,    const unsigned int libmesh_dbg_var(i),
//                            const Point &,   const bool) {
//   libmesh_assert_less (i, 1);  return 1.;}
// 
template <> double  MGFE<0,LAGRANGE>::shape(
  const ElemType,
  const Order,
  const unsigned int i,
  const double p[]
){ assert(i< 1); return 1.;}

//  L2-LAGRANGE
// template <>
// Real MGFE<0,L2_LAGRANGE>::shape(const ElemType,      const Order,      const unsigned int libmesh_dbg_var(i),
//                               const Point &) {libmesh_assert_less (i, 1);   return 1.;
// }
// template <> Real MGFE<0,L2_LAGRANGE>::shape(const Elem *,  const Order,  const unsigned int libmesh_dbg_var(i), const Point &,  const bool)
// {  libmesh_assert_less (i, 1); return 1.;
// }

// template <> Real MGFE<0,L2_LAGRANGE>::shape(const FEType,
//                               const Elem *,
//                               const unsigned int libmesh_dbg_var(i),
//                               const Point &,
//                               const bool)
// {   libmesh_assert_less (i, 1);   return 1.; }
// 


// 
// template <>
// Real MGFE<0,L2_LAGRANGE>::shape_deriv(const ElemType,     const Order,    const unsigned int,
//                                     const unsigned int,      const Point &)
// {libmesh_error_msg("No spatial derivatives in 0D!"); return 0.; }
// 
// 
// 
// template <> Real MGFE<0,L2_LAGRANGE>::shape_deriv(const Elem *,     const Order,  const unsigned int,
//                                     const unsigned int, const Point &,  const bool)
// {
//   libmesh_error_msg("No spatial derivatives in 0D!");
//   return 0.;
// }
// 
// 
// 
template <> double MGFE<0,LAGRANGE>::shape_deriv(
  const ElemType,
  const Order,
  const unsigned int,
   const unsigned int,
  const double [])
{ std::cerr<<"No spatial derivatives in 0D!"; return 0.;
}
// 
// 
// 
// template <>
// Real MGFE<0,LAGRANGE>::shape_deriv(const Elem *,
//                                  const Order,
//                                  const unsigned int,
//                                  const unsigned int,
//                                  const Point &,
//                                  const bool)
// {
//   libmesh_error_msg("No spatial derivatives in 0D!");
//   return 0.;
// }
// 
// 
// 
// template <>
// Real MGFE<0,L2_LAGRANGE>::shape_deriv(const FEType,
//                                     const Elem *,
//                                     const unsigned int,
//                                     const unsigned int,
//                                     const Point &,
//                                     const bool)
// {
//   libmesh_error_msg("No spatial derivatives in 0D!");
//   return 0.;
// }
// 
// 
// template <>
// Real MGFE<0,LAGRANGE>::shape_deriv(const FEType,
//                                  const Elem *,
//                                  const unsigned int,
//                                  const unsigned int,
//                                  const Point &,
//                                  const bool)
// {
//   libmesh_error_msg("No spatial derivatives in 0D!");
//   return 0.;
// }
// 
// 
// 
// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// 
// template <>
// double MGFE<0,L2_LAGRANGE>::shape_second_deriv(const ElemType,
//                                            const Order,
//                                            const unsigned int,
//                                            const unsigned int,
//                                            const double &)
// {
//   std::cerr<<"No spatial derivatives in 0D!";
//   return 0.;
// }
// 
// 
// 
// template <>
// Real MGFE<0,L2_LAGRANGE>::shape_second_deriv(const Elem *,
//                                            const Order,
//                                            const unsigned int,
//                                            const unsigned int,
//                                            const Point &,
//                                            const bool)
// {
//   libmesh_error_msg("No spatial derivatives in 0D!");
//   return 0.;
// }
// 
 template <>  double MGFE<0,LAGRANGE>::shape_second_deriv(const ElemType,
                                         const Order,
                                         const unsigned int,
                                         const unsigned int,
                                         const double p[])
 {
   std::cerr<<"No spatial derivatives in 0D!";
   return 0.;
 }
// 
// 
// 
// template <>
// Real MGFE<0,LAGRANGE>::shape_second_deriv(const Elem *,
//                                         const Order,
//                                         const unsigned int,
//                                         const unsigned int,
//                                         const Point &,
//                                         const bool)
// {
//   libmesh_error_msg("No spatial derivatives in 0D!");
//   return 0.;
// }
// 
// 
// template <>
// Real MGFE<0,L2_LAGRANGE>::shape_second_deriv(const FEType,
//                                            const Elem *,
//                                            const unsigned int,
//                                            const unsigned int,
//                                            const Point &,
//                                            const bool)
// {
//   libmesh_error_msg("No spatial derivatives in 0D!");
//   return 0.;
// }
// 
// 
// template <>
// Real MGFE<0,LAGRANGE>::shape_second_deriv(const FEType,
//                                         const Elem *,
//                                         const unsigned int,
//                                         const unsigned int,
//                                         const Point &,
//                                         const bool)
// {
//   libmesh_error_msg("No spatial derivatives in 0D!");
//   return 0.;
// }
// 
// #endif
// 
// 
// } // namespace libMesh


// // TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
// // Basic TRI/TETRA
// // TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
// 
// // =======================================================================================
// template <>  double MGFE<0,LAGRANGE>::Tri_2d_LinearPhi(
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
// template <>  double MGFE<0,LAGRANGE>::Tri_2d_QuadraticPhi(
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
// template <>  double MGFE<0,LAGRANGE>::Tri_2d_LinearDerPhi(int nPhi, double */*point[]*/, int dir) {
//   double PhiDer;
//   PhiDer = _CooTriEl[nPhi * 3 + dir];
//   return PhiDer;
// }
// // =======================================================================================
// template <>  double MGFE<0,LAGRANGE>::Tri_2d_QuadraticDerPhi(int nPhi, double point[], int dir) {
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
// template <>  double MGFE<0,LAGRANGE>::Tri_2d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2) {
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
// // =======================================================================================
// template <>  double MGFE<0,LAGRANGE>::Tri_3d_LinearPhi(int nPhi, double point[]) {
//   double PhiVal;
//   int a = _CooTetra4[nPhi * 4];
//   int b = _CooTetra4[nPhi * 4 + 1];
//   int c = _CooTetra4[nPhi * 4 + 2];
//   int d = _CooTetra4[nPhi * 4 + 3];
//   PhiVal = a * point[0] + b * point[1] + c * point[2] + d;
//   return PhiVal;
// }
// // =======================================================================================
// template <>  double MGFE<0,LAGRANGE>::Tri_3d_LinearDerPhi(int nPhi, double */*point[]*/, int dir) {
//   double PhiDer;
//   PhiDer = _CooTetra4[nPhi * 4 + dir];
//   return PhiDer;
// }
// // =======================================================================================
// template <>  double MGFE<0,LAGRANGE>::Tri_3d_QuadraticPhi(int nPhi, double point[]) {
//   double PhiVal;
//   int l1 = _CooTetra10[nPhi * 4];
//   int l2 = _CooTetra10[nPhi * 4 + 1];
//   int a = _CooTetra10[nPhi * 4 + 2];
//   int b = _CooTetra10[nPhi * 4 + 3];
//   double L1, L2;
//   L1 = Tri_3d_LinearPhi(l1, point);
//   L2 = Tri_3d_LinearPhi(l2, point);
//   PhiVal = L1 * (a * L2 + b);
//   return PhiVal;
// }
// // =======================================================================================
// template <>  double MGFE<0,LAGRANGE>::Tri_3d_QuadraticDerPhi(int nPhi, double point[], int dir) {
//   double PhiDer;
//   int l1 = _CooTetra10[nPhi * 4];
//   int l2 = _CooTetra10[nPhi * 4 + 1];
//   int a = _CooTetra10[nPhi * 4 + 2];
//   int b = _CooTetra10[nPhi * 4 + 3];
//   double L1, L2, dL1, dL2;
//   L1 = Tri_3d_LinearPhi(l1, point);
//   L2 = Tri_3d_LinearPhi(l2, point);
//   dL1 = Tri_3d_LinearDerPhi(l1, point, dir);
//   dL2 = Tri_3d_LinearDerPhi(l2, point, dir);
//   PhiDer = dL1 * (a * L2 + b) + L1 * a * dL2;
//   return PhiDer;
// }
// // =======================================================================================
// template <>  double MGFE<0,LAGRANGE>::Tri_3d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2) {
//   double PhiDer;
//   int l1 = _CooTetra10[nPhi * 4];
//   int l2 = _CooTetra10[nPhi * 4 + 1];
//   int a = _CooTetra10[nPhi * 4 + 2];
//   double dL1_1, dL2_1, dL1_2, dL2_2;
//   dL1_1 = Tri_3d_LinearDerPhi(l1, point, dir1);
//   dL2_1 = Tri_3d_LinearDerPhi(l2, point, dir1);
//   dL1_2 = Tri_3d_LinearDerPhi(l1, point, dir2);
//   dL2_2 = Tri_3d_LinearDerPhi(l2, point, dir2);
//   PhiDer = a * (dL1_1 * dL2_2 + dL1_2 * dL2_1);
//   return PhiDer;
// }






// // QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
// // Basic QUAD/HEX
// // QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
// // =======================================================================================
// template <>
//  double MGFE<0,LAGRANGE>::Edge_Quad_Phi(int PhiCoeff, double Coordinate) {
//   double PhiVal;
//   PhiVal = (1 - 0.5 * std::fabs(PhiCoeff)) * ((2 * fabs(PhiCoeff) - 1) * Coordinate * Coordinate +
//       PhiCoeff * Coordinate + (1 - fabs(PhiCoeff)));
//   return PhiVal;
// }
// // ======================================================================================
// template <>
//  double MGFE<0,LAGRANGE>::Edge_Quad_DPhi(int PhiCoeff, double Coordinate) {
//   double DPhiVal;
//   DPhiVal = (1 - 0.5 * fabs(PhiCoeff)) * (2. * (2 * fabs(PhiCoeff) - 1) * Coordinate + PhiCoeff);
//   return DPhiVal;
// }
// // ======================================================================================
// template <>
//  double MGFE<0,LAGRANGE>::Edge_Quad_D2Phi(int PhiCoeff, double /*Coordinate*/) {
//   double D2PhiVal;
//   D2PhiVal = (1 - 0.5 * fabs(PhiCoeff)) * (2. * (2 * fabs(PhiCoeff) - 1));
//   return D2PhiVal;
// }
// // ======================================================================================
// template <>
//  double MGFE<0,LAGRANGE>::Edge_Lin_Phi(int PhiCoeff, double Coordinate) {
//   double PhiVal;
//   PhiVal = 0.5 * (1 + PhiCoeff * Coordinate);
//   return PhiVal;
// }
// // ======================================================================================
// template <>
//  double MGFE<0,LAGRANGE>::Edge_Lin_DPhi(int PhiCoeff, double /*Coordinate*/) {
//   double DPhiVal;
//   DPhiVal = 0.5 * PhiCoeff;
//   return DPhiVal;
// }
// // ======================================================================================
// template <>
//  double MGFE<0,LAGRANGE>::Rec_Lin_Phi(int nPhi, double point[], int dimension) {
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
// template <>
//  double MGFE<0,LAGRANGE>::Rec_Quad_Phi(int nPhi, double point[], int dimension) {
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
// template <>
//  double MGFE<0,LAGRANGE>::Rec_Lin_DPhi(int nPhi, double point[], int dimension, int DirDer) {
//   int lambda_i;
// 
//   if(dimension == 1) {
//     lambda_i = _CooE3[nPhi];
//   } else if(dimension == 2) {
//     lambda_i = _CooQ9[nPhi + DirDer * _Q9Off];
//   } else if(dimension == 3) {
//     lambda_i = _CooH27[nPhi + DirDer * _H27Off];
//   } else {
//     printf("\033[1;31m MGFE1::Rec_Lin_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//   }
// 
//   double PhiDer = Edge_Lin_DPhi(lambda_i, point[DirDer]);
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
//       printf("\033[1;31m MGFE1::Rec_Lin_DPhi unkown phi for dimension %d \n \033[0m", dimension);
//     }
// 
//     PhiDer *= Edge_Lin_Phi(lambda_i, point[direction]);
//   }
// 
//   return PhiDer;
// }
// // ======================================================================================
// template <>
//  double MGFE<0,LAGRANGE>::Rec_Quad_DPhi(int nPhi, double point[], int dimension, int DirDer) {
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
// template <>
//  double MGFE<0,LAGRANGE>::Rec_Lin_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2) {
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
// template <>
//  double MGFE<0,LAGRANGE>::Rec_Quad_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2) {
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
// 
// 
// 
// template <>
//  double MGFE<0,LAGRANGE>::FirstDerivateOfLocalPhi(int nPhi, double point[], int Dimension, int DirDer, int FamilyType) {
//   double PhiDer;
// 
//   if(FamilyType == 0) {   // TRIANGULAR ELEMENTS
//     if(_order == 1) {
//       if(Dimension == 1) {
//         PhiDer = Rec_Lin_DPhi(nPhi, point, Dimension, DirDer);
//       } else if(Dimension == 2) {
//         PhiDer = Tri_2d_LinearDerPhi(nPhi, point, DirDer);
//       } else if(Dimension == 3) {
//         PhiDer = Tri_3d_LinearDerPhi(nPhi, point, DirDer);
//       }
//     }
// 
//     if(_order == 2) {
//       if(Dimension == 1) {
//         PhiDer = Rec_Quad_DPhi(nPhi, point, Dimension, DirDer);
//       } else if(Dimension == 2) {
//         PhiDer = Tri_2d_QuadraticDerPhi(nPhi, point, DirDer);
//       } else if(Dimension == 3) {
//         PhiDer = Tri_3d_QuadraticDerPhi(nPhi, point, DirDer);
//       }
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
// 
// 
// template <> void MGFE<0,LAGRANGE>::init_pie_rec() { return;}
// template <> void MGFE<0,LAGRANGE>::init_lin_rec() { return;}
// template <> void MGFE<0,LAGRANGE>::init_lin_rec_rt() { return;}
// template <> void MGFE<0,LAGRANGE>::init_lin() { return;}
// template <> void MGFE<0,LAGRANGE>::init_pie() { return;}
// 
// 
// 
// // LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
// template <> void MGFE<0,LAGRANGE>::init_lin_tri(
// ) {// LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
// 
// 
//   return;
// }
// 
template <> double   MGFE<0,LAGRANGE>:: compute_fe_face_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz[] // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
 ){
return 0.;}

template <> double   MGFE<0,LAGRANGE>:: compute_fe_gauss_points(
   const int  gp, // gaussian point gp
   const double xyz[], // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
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
  
//   for(int l=0;l<_NoShape;l++)       {
//     phi[l]=_phi_map[l* _NoGauss  +gp] ;// shape
//    dphidxi[l]=_dphidxez_map[(l)*_NoGauss + gp];  phi[l]=dphidxi[l];
// }
//   double InvJac=0;
//   double x_xi = 0.; 
//   for(int s = 0; s < _NoShape; s++) {
//     x_xi += xyz[s] *  _dphidxez_map[s * _NoGauss + gp];
//   }
//   double det = x_xi; InvJac =1./det;
//   for(int l=0;l<_NoShape;l++)  {  dphidx[l]=InvJac* dphidxi[l];  }
  
  
   
  return 0.;
}

// 
// // template <> double MGFE<0,LAGRANGE>::ComputeInverseMatrix(double *,double *,int) { return 0.;}
template <> double MGFE<0,LAGRANGE>::ComputeInverseMatrix_sur(double *,double *) { return 0.;}
// // template <> double MGFE<0,LAGRANGE>::FirstDerivateOfLocalPhi(int , double *,int,int,int) { return 0.;}
// 

// 
// 
// 
// // ================================================================================================
// //    QUAD
// // ============================================================================
// template<> void  MGFE<0,LAGRANGE>::init_qua_rec(
// ) {//==========================================================================
// ///  The \p Hex27 is an element in 3D composed of 27 nodes.
// 
// 
//   return;
// }
// 
// 
// // ================================================================================
// template<> void  MGFE<0,LAGRANGE>::init_qua_tri() {
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
//  
//   return;
// }


// // ================================================================================================
// /// This function generates the Lagrangian quad shape functions
// template<> void  MGFE<0,LAGRANGE>::init_qua() {  // ================================
//   if(_GeomType==1) init_qua_rec();  
//   if(_GeomType==0) init_qua_tri();
//   return;
// }

// ================================================================================================
// ================================================================================================
// ================================================================================================



// ================================================
/// This function computes the normal at the gauss point
template <>
void MGFE<0,LAGRANGE>::normal_g(
const double* xx,  // all surface coordinates <-
double* normal_g   // normal ->
) const {          // ======================================


  for(int idim = 0; idim < _dim; idim++) { normal_g[idim] = 0.; }

  return;
}
//=============================================================
template <>
void  MGFE<0,LAGRANGE>::Oxy_face(
const double* xx,   ///< coordinates
const double x_c[], ///< central point
double* normal_g,   ///<  i_n,i_t1,i_t2 face Oxyz
double* tg01,
double* tg03
) const {            // ======================================

  for(int idim = 0; idim < _dim; idim++) {normal_g[idim] = 0.; tg01[idim]=0.;tg03[idim]=0.; }
  return;
}

// =========================================================================================================================
template <> void MGFE<0,LAGRANGE>::set_par(std::string fileh5,int dim_in,int order_in,int shapes_quad,int geom_type_in,int fam_type_in)  
{_dim=dim_in,  _NoShape=shapes_quad,   _GeomType=geom_type_in,     _FamType=fam_type_in;
   femfileh5_=fileh5;
     ///< family type  0=Lagrangian(C^0) 1=Raviart Thomas (L^2-C^0)
  // // Fem order
  // // =============================
  // //  ------------------------------------------------------------------------
  // // gauss points (_NoGauss) ------------------------------------------------------
  //  _NoGauss=1; 
  // 
  // // number of nodes (_NoElNodes) in the reference quad element
  // _NoElNodes=1;  
   _NoGauss=1; // assert(_qrule); _NoGauss=_qrule->n_points();
 // if(_NoGauss==0) {   std::cout << "Quadrature not assigned"; abort();}

  // number of nodes (_NoElNodes) in the reference quad element
   ElemType  elem_name=_fe_type.fetype; 
   // ElemType  geom_elem_name=get_ref_elemtype(elem_name); 
   Order elem_order= _fe_type.order ;
    _order= (int)elem_order;
  // number of nodes (_NoElNodes) in the reference EDGE3 element
    // _NoElNodes= get_refspace_n_nodes( geom_elem_name);// mesh nodes
  // _NoShape=get_refspace_n_nodes(elem_name); // element nodes 
    
  //  -----------------------------------------------------------------------------------
  // x node  coordinates   (_xnodes[xp1,xp2,xp3,...yp1,yp2,yp3,....zp1,zp2,zp3,...] )----
  _xnodes=new double[1];
  // shape functions at nodes (_phi_map_nodes[sh1(x1),sh1(x2),....  sh2(x1),sh2(x2),.... ])
  _phi_map_nodes = nullptr;//new double[_NoShape* _NoElNodes];
  // 1st derivatives at nodes x _dphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxez_map_nodes = nullptr;//new double[_dim * _NoShape * _NoElNodes];
  // element dofs  (reference to nodes) --------------------------------------------------
  _fem_conn=new int[1];
  // gaussian (g) coordinates (_xg[xg1,xg2,xg3,...yg1,yg2,yg3,....zg1,zg2,zg3,...] )-----
  _xg=new double[1];
  // weights [wg1,wg2,wg3,.......]
  _weight = new double[_NoGauss];
  // shape functions phi_map[sh1(xg1),sh1(xg2),....  sh2(xg1),sh2(xg2),.... ]
  _phi_map = nullptr;//new double[ _NoGauss*_NoShape ];
  // 1st derivatives at g _dphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxez_map = nullptr;//new double[_dim* _NoShape* _NoGauss];
  // 2nd derivatives at g _ddphi_map[dphi/de1,.....,dphi/de2,....,dphi/de3]
  _dphidxx_map = nullptr;//new double[_dim * _dim  * _NoShape * _NoGauss];
  //   -------------------------------------------------------------------------------------
  // dofs -----------------------------------------------------------------------
 _fem_conn[0]=0;
 // double *xnodes=new double [_NoElNodes*_dim];
 // double *point=new double[_dim];
    //
    // ---------------------------------------------------------------------------------------
    // nodes ----------------------------------------------------------------------------
  // get_refspace_nodes(geom_elem_name,xnodes);// nodes[(x,y)_0,(x,y)_1.... ]
    // _xnodes[0]=0.;    _xnodes[1]=1.;    _xnodes[2]=0.5;
    // for(int i = 0; i < _NoElNodes; i++) {  // points
    //   for(int dir = 0; dir < _dim; dir++) { 
        _xnodes[0]=0.; 
        // point[0] = _xnodes[0];
      // }

    // for(int j = 0; j < _NoShape; j++) {  // shapse
    //   //  _phi_map_nodes[phi1(xp1),phi1(xp2),....phi2(xp1),phi2(xp2),...  ]
    //   _phi_map_nodes[j*_NoElNodes+i]=1.;//Rec_Quad_Phi(j, point,_dim);
    // 
    //   // _dphidxez_map_nodes[dphi1(xp1)/de1,dphi1(xp2)/de1,....dphi2(xp1)/de1,dphi2(xp2)/e1,...  ]
    //   for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
    //     _dphidxez_map_nodes[(j + dir *_NoShape)*_NoElNodes+i] =0.;// Rec_Quad_DPhi(j,point,_dim,dir);     // dphi/dxi
    //   }
    // }
  // }
   // delete[]xnodes;
// gaussian points --------------------------------------------------------->
 // for(int i = 0; i < _NoGauss; i++) {
        
    _weight[0] =1.;//_qrule->w(i); // weight
    
     // _qrule->qp(point,i);     // ith-gaussian point ->point(i)
    // for(int j = 0; j < _NoShape; j++) {     // loop over test function id
      // phi[j*_NoGauss+i] =1.;
      // for(int dir = 0; dir < _dim; dir++) {  //  derivatives ++++++
      //   _dphidxez_map[(j + dir * _NoShape) * _NoGauss + i] =0.;
      //   for(int dir2 = 0; dir2 < _dim; dir2++)// dir +dir2=0->d^2phi/dxi^2; dir +dir2=2->d^2phi/deta^2; dir +dir2=1->d^2phi/deta dxi;
      //     _dphidxx_map[(j + (dir * _dim + dir2) * _NoShape) * _NoGauss + i] =0.;
      // }// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // }
  // } 
  // delete []point;
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


template <> void  MGFE<0,LAGRANGE>::get_dphi_on_given_node(const int dim, double ElemCoords[], double CanPos[], double dphi[]) {
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

// void MGFE1::get_dphi_on_given_nodeG(
// const int dim, double ElemCoords[], double CanPos[], double dphi[], int FamilyType) {
//   double* InvJac = new double[dim * dim];
//   double* gradphi_g = new double[dim];
//   int nShape = _GNoShape[FamilyType][dim - 1];
//   JacOnGivenCanCoords(dim, ElemCoords, CanPos, InvJac, FamilyType, nShape);
//
//   for(int eln = 0; eln < nShape; eln++) {
//     for(int idim = 0; idim < dim; idim++) {
//       gradphi_g[idim] = FirstDerivateOfLocalPhi(eln, CanPos, dim, idim, FamilyType);
//     }
//
//     for(int idim = 0; idim < dim; idim++) {
//       double sum = 0.;
//       for(int jdim = 0; jdim < dim; jdim++) { sum += InvJac[jdim + idim * dim] * gradphi_g[jdim]; }
//       dphi[eln + idim * nShape] = sum;
//     }
//   }
//
//   delete[] gradphi_g;
//   delete[] InvJac;
//   return;
// }

// ==================================================================
/// This function computes the shape values at the gauss point qp
template <> void  MGFE<0,LAGRANGE>::get_phi_gl_g(
const int qp,    // gaussian point <-
double phi[]     // shape functions ->
) {                  // =================================================
  for(int ish = 0; ish < _NoShape; ish++) phi[ish] = phi[ish * _NoGauss + qp];
  return;
}

// ==================================================================
/// This function computes the shape values at the gauss point qp
template <> void  MGFE<0,LAGRANGE>::get_phi_gl_g(
const int qp,             // gaussian point   <-
std::vector<double>& phi  // shape functions  ->
) {                           // =====================================
  for(int ish = 0; ish < _NoShape; ish++) phi[ish] = phi[ish * _NoGauss + qp];

  return;
}
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
template <> void  MGFE<0,LAGRANGE>::get_dphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double dphi[]           // global derivatives ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  double gradphi_g[1];  // temp grad phi

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
//   double gradphi_g[MGFE1::MDIM_FEM];  // temp grad phi
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
template <> void  MGFE<0,LAGRANGE>::get_ddphi_gl_g(
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double ddphi[]          // global derivatives ->
) {                         // =========================================
  int NDOF_FEM=1;
  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset
  double ddphi_loc[NDOF_FEM * 1 * 1];

  for(int eln = 0; eln < el_nnodes; eln++) {   // LOOP OVER NODES ===================================
    const int shift = eln * _dim * _dim;       // offset for node derivatives
    const int lqp = eln * el_ngauss + qp;

    for(int idim = 0; idim < _dim * _dim; idim++) {
      ddphi_loc[shift + idim] = _dphidxx_map[lqp + idim * goffset];
    }

    double Hess_tmp[1 * 1], JacTf[1*1];

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
template <> void  MGFE<0,LAGRANGE>::get_dphi_node(
const int /*kdim*/,         // dimension <-
const int node,         // nodal point <-
const double InvJac[],  // Jacobean <-
double dphi[]           // global derivatives at the nodal point ->
) {                         // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of gauss points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

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
//   double gradphi_g[MGFE1::MDIM_FEM];  // temp grad phi
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
template <> void  MGFE<0,LAGRANGE>::get_dphi_gl_g(
const int /*kdim*/,         // dimension <-
const int qp,           // gaussian point <
const double InvJac[],  // Jacobean
double dphi[],          // global derivatives ->
int sdim) {             // =========================================

  const int el_nnodes = _NoShape;   // # of shape functions
  const int el_ngauss = _NoGauss;  // # of guass points
  const int goffset = el_nnodes * el_ngauss;  // gauss offset

  double gradphi_g[1];  // temp grad phi

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

// =============================================================
/// Shape functions derivatives at the gauss point qp
template <> void  MGFE<0,LAGRANGE>::get_dphi_gl_g(
const int /*kdim*/,            // dimension <-
const int qp,              // gaussian point <-
const double InvJac[],     // Jacobean
std::vector<double>& dphi  // global derivatives ->
) {                            // =========================================

  const int el_nnodes = _NoShape;             // # of shape functions
  const int el_ngauss = _NoGauss;            // # of guass points
  const int goffset = el_nnodes * _NoGauss;  // gauss offset

  double dphidxi_g[1]; dphidxi_g[0]=0.;  // temp grad phi

  

  return;
}


// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// // Basic Raviart Thomas non continuous functions
// // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// //===============================================================================================
//  template<>
// std::array<double,3 >   MGFE<0,LAGRANGE>::Rec_RT_Phi(int iface, double x[],int dim) {
// 
// 
//   return {0.,0.,0.};
// }
// 
// //===============================================================================================
// 
// template<>
// std::array<double,3>   MGFE<0,LAGRANGE>::Rec_RT_DPhi(int /*iface*/, double x[],int dim,  int dir) {
//   return  {0.,0.,0.};  
// }
// template <> double   MGFE<0,LAGRANGE>::compute_fe_gauss_points(
//    const int  gp, // gaussian point gp
//    double xyz[] // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
//  ){
//   
//   return det 0.;
// }
// 
// template <> double   MGFE<0,LAGRANGE>::compute_fe_face_gauss_points(
//    const int  gp, // gaussian point gp
//    double xyz[] // x_elem[xyz(1),xyz(2), xyz(3),   xyz(_Nshape)] 
//  ){
//   
//   return det 0.;
// }
//===============================================================================================
template<>
double  MGFE<0,LAGRANGE>::ComputeInverseMatrix(double Matrix[], double InvMatrix[], int dim) {
  double det = 0.;

 
  
  return det;
}

//===============================================================================================
template<>
double  MGFE<0,LAGRANGE>::Tj_interp_lin(double *u,int j, double *xi,int ndof,int idim) {
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
template<>
double  MGFE<0,LAGRANGE>::dTjdxi_lin(double *u,int j, double *xi,int i,int ndof,int idim) {
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



