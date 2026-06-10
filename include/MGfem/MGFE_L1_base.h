#pragma once


// C++ includes
#include <cstddef>
#include <vector>
#include <memory>
#include <cassert>


#include "MGFE_L0_abstract.h"
// #include "MGFE_L0_II_enum_order.h"

/**
 * This class forms the foundation from which generic finite
 * elements may be derived.  In the current implementation the
 * templated derived class \p FE offers a wide variety of commonly
 * used finite element concepts.  Check there for details.
 *
 * Use the \p MGFEbase<OutputType>::build() method to create an
 * object of any of the derived classes which is compatible with
 * OutputType.
 *
 * \author Benjamin S. Kirk
 * \date 2002  
 */
// template <typename OutputType> 

class MGFEbase : public MGFE0 {
 
 
 public:
 // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // static constexpr int MDIM_FEM=3; // max geom dim FEM (volume(3))
    //   static constexpr int NDOF_FEM=27; // max dof FEM (volume(3) LAGRANGE quad=27)
    //     static constexpr int NDOF_FEMB=9; // max dof FEM (volume(2) LAGRANGE quad=9)
    //Input parameter ----------------------------------------------------------------------------
    std::string femfileh5_;
     int _dim;      ///< Dimension (3 or 2) used in the 3D FEM (1< _dim < MDIM_FEM)
     int _order;    ///< Order of the shape functions 0=cpnst 1=linear 2=quad
     int _NoShape;  ///< Number of shape functions of the fem
     int _NoElNodes;///< Number of nodes in the reference element (always ref quad)
     int _GeomType; ///< Geom type 1=HEXAHEDRAL 0=TETRAHEDRAL
     int _FamType;  ///< Family  (geom 1: Lagrangian =0  RT =1)

    // Gaussian points -----------------------------------------------------------------------------
    int     _NoGauss;  ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
    double*  _xg;                 ///< gaussian points
    double* _weight;              ///< Weight
    double* _phi_map;             ///< Shape functions
    double* _dphidxez_map;        ///< Shape derivative functions in gaussian points
    double* _dphidxx_map;         ///< Second order shape derivatives in Gaussian points

    // Nodes (NDOF_FEM) -----------------------------------------------------------------------------
    double *_xnodes;              ///< xnode (hex or tetra) in (xi,eta,seta)
    int *_fem_conn;              ///< dof node coords (respect to _xnodes)
    // int  _deg;                  ///< Polinomial degree
    double* _phi_map_nodes;             ///< Shape functions
    double* _dphidxez_map_nodes;  ///< Shape derivative functions in nodal points
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
protected:
  /** Constructor.  Optionally initializes required data structures.  Protected so that this base class
   * cannot be explicitly instantiated.  */
  MGFEbase (const unsigned int dim, const FEType & fet);

public:

  /* Destructor. */
  virtual ~MGFEbase();

  /* Builds a specific finite element type.  */
  static std::unique_ptr<MGFEbase> build (const unsigned int dim, const FEType & type);

  /* Convenient typedefs for gradients of output, hessians of output, and potentially-complex-valued versions of same. */
  // typedef OutputType                                                      OutputShape;
  // typedef typename TensorTools::IncrementRank<OutputShape>::type          OutputGradient;
  // typedef typename TensorTools::IncrementRank<OutputGradient>::type       OutputTensor;
  // typedef typename TensorTools::DecrementRank<OutputShape>::type          OutputDivergence;
  // typedef typename TensorTools::MakeNumber<OutputShape>::type             OutputNumber;
  // typedef typename TensorTools::IncrementRank<OutputNumber>::type         OutputNumberGradient;
  // typedef typename TensorTools::IncrementRank<OutputNumberGradient>::type OutputNumberTensor;
  // typedef typename TensorTools::DecrementRank<OutputNumber>::type         OutputNumberDivergence;
  virtual void set_par(std::string fileh5,int dim_in,int order_in,int shapes_quad,int geom_type_in,int fam_type_in)=0;
  
  
 // get function
 virtual   int GetFamilyType(int elem_dof, int dim)=0;  // the function returns 0 or 1
  //Input parameter ----------------------------------------------------------------------------
 virtual const int get_dim() =0     ;// {return _dim;}  ///< Dimension (3 or 2) used in the 3D FEM (1< _dim < MDIM_FEM)
 virtual const int get_order() =0   ;// {return _order;}///< Order of the shape functions 0=cpnst 1=linear 2=quad
 virtual const int get_NoShape() =0 ;// {return _NoShape;}  ///< Number of shape functions of the fem
 virtual const int get_NoElNodes()=0 ;//{return _NoElNodes;}  ///< Number of nodes in the reference element (always ref quad)
 virtual const int get_GeomType()=0;//  {return _GeomType;}   ///< Geom type 1=HEXAHEDRAL 0=TETRAHEDRAL
 virtual const int get_FamType()=0 ;//  {return _FamType;}    ///< Family  (geom 1: Lagrangian =0  RT =1)
  
   // Gaussian points -----------------------------------------------------------------------------
  virtual const int get_NoGauss()=0;//{ return _NoGauss;}  ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
  virtual  const double  get_xg(int gp)=0;//{          return  *_xg;}          ///< gaussian points
  virtual const double   get_weight(int gp)=0;//            ///< Weight
  virtual const double & get_xg()          {return  *_xg;}          ///< gaussian points
  virtual const double & get_weight()      {return  *_weight;} ///< gaussian weights
  virtual const double & get_phi_map()     {return  *_phi_map;}           ///< Shape functions
  virtual const double & get_dphidxez_map(){return  *_dphidxez_map;}      ///< Shape derivative functions in gaussian points
  virtual const double & get_dphidxx_map() {return  *_dphidxx_map;}       ///< Second order shape derivatives in Gaussian points
 
 // xfem ======================================================
   // inline double get_weight(int qp) {return _GaussWeights[qp];}
  virtual   inline double get_xi(int qp) {return 0.;}
  virtual   inline double get_eta(int qp) {return 0.;}
  virtual   inline double get_zeta(int qp) {return 0.;}
 
 // Nodes (NDOF_FEM) -----------------------------------------------------------------------------
 virtual  const double &  get_xnodes()=0;//{return   *_xnodes;}              ///< xnode (hex or tetra) in (xi,eta,seta)
 virtual const int    &  get_fem_conn()=0;//{ return *_fem_conn;}              ///< dof node coords (respect to _xnodes)
 // virtual const int       get_deg()=0;//{return _deg};                  ///< Polinomial degree
 virtual const double &   get_phi_map_nodes()=0;//{return   *_phi_map_nodes;   }          ///< Shape functions
 virtual  const double & get_dphidxez_map_nodes()=0;//{return   *_dphidxez_map_nodes;}  ///< Shape derivative functions in nodal points
  
  // virtual double Jac(const int ng, double x[], double InvJac[]);
  // virtual double& normal_g(const double* xx, const double x_c[],double* normal_g) const;          
// virtual      void init_qua()=0;  ///< Generates the Lagrangian quad shape functions
// virtual    void init_lin()=0;  ///< Generates the Lagrangian linear shape functions
// virtual    void init_pie()=0;  ///< Generates the Lagrangian piecewise shape functions
// virtual    void init_pie_rec()=0;
//  virtual   void init_pie_tri()=0;
// 
//  virtual   void init_qua_tri()=0;
// //   void init_qua_Gtri();
// virtual    void init_qua_rec()=0;
// //   void init_qua_Grec();
// 
//  virtual   void init_lin_tri()=0;
// //   void init_lin_Gtri();
//  virtual   void init_lin_rec()=0;
//  virtual   void init_lin_rec_rt()=0;
//   void init_lin_Grec();
    /// Clear data substructures -------------------------------------------------------------
// virtual    void clear()=0;
    // ~MGFE();

 
 
 

    // -------------------------------------------------------
    // Triangle - based elements:

// virtual    double Tri_2d_LinearPhi(int nPhi, double point[])=0;
//  virtual   double Tri_2d_LinearDerPhi(int nPhi, double point[], int dir)=0;
//  virtual   double Tri_2d_QuadraticPhi(int nPhi, double point[])=0;
//  virtual   double Tri_2d_QuadraticDerPhi(int nPhi, double point[], int dir)=0;
//  virtual   double Tri_2d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2)=0;
// 
//  virtual   double Tri_3d_LinearPhi(int nPhi, double point[])=0;
// virtual    double Tri_3d_LinearDerPhi(int nPhi, double point[], int dir)=0;
// virtual    double Tri_3d_QuadraticPhi(int nPhi, double point[])=0;
//  virtual   double Tri_3d_QuadraticDerPhi(int nPhi, double point[], int dir)=0;
// virtual    double Tri_3d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2)=0;
// 
//     // -------------------------------------------------------
//     // Quadrangle - based elements:
//  virtual   double Edge_Lin_Phi(int PhiCoeff, double Coordinate)=0;
//  virtual   double Edge_Lin_DPhi(int PhiCoeff, double Coordinate)=0;
//  virtual   double Edge_Quad_Phi(int PhiCoeff, double Coordinate)=0;
// virtual    double Edge_Quad_DPhi(int PhiCoeff, double Coordinate)=0;
//  virtual   double Edge_Quad_D2Phi(int PhiCoeff, double Coordinate)=0;
// 
// virtual    double Rec_Lin_Phi(int nPhi, double point[], int dimension)=0;
//  virtual   double Rec_Lin_DPhi(int nPhi, double point[], int dimension, int DirDer)=0;
//  virtual   double Rec_Lin_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2)=0;
//  virtual   double Rec_Quad_Phi(int nPhi, double point[], int dimension)=0;
//  virtual   double Rec_Quad_DPhi(int nPhi, double point[], int dimension, int DirDer)=0;
//  virtual   double Rec_Quad_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2)=0;
// 
//  virtual   double FirstDerivateOfLocalPhi(int nPhi, double point[], int dimension, int DirDer, int FamilyType)=0;
// 
// // Raviart Thomas function
// virtual    std::array<double,3> Rec_RT_Phi(int iface, double x[],int dim)=0;
//  virtual   std::array<double,3> Rec_RT_DPhi(int iface, double x[],int dim,int dir)=0;
//=================================================================================================
    /// \NAME RETURN SHAPE DERIVATIVE FUNCTIONS AT GAUSS POINTS ---------------------------
    ///< \param[in]  <dim>    Dimension
    ///< \param[in]  <qp>     Gaussian point
    ///< \param[in]  <InvJac> Jacobean
    ///< \param[out] <dphi>   Derivative

virtual    void get_dphi_gl_g(const int qp, const double InvJac[], double dphi[])=0;
virtual    void get_dphi_gl_g(const int dim, const int qp, const double InvJac[], double dphi[], int sdim)=0;
virtual    void get_dphi_gl_g(const int dim, const int qp, const double InvJac[], std::vector<double>& dphi)=0;
// =============================================================
/// Shape functions derivatives at the gauss point qp
virtual void get_dphi_gl_g(
 const int kdim,         // dimension <-
double point[],           // gaussian point <
const double InvJac[],  // Jacobean
double dphi[],          // global derivatives ->
int sdim,
int fam_type){}             // =========================================



virtual    void get_dphi_node(const int dim, const int qp, const double InvJac[], double dphi[])=0;
virtual   void get_ddphi_gl_g(const int qp, const double InvJac[], double ddphi[])=0;
virtual   void get_dphi_on_given_node(const int dim, double ElemCoords[], double CanPos[], double dphi[])=0;
 // virtual    void get_dphi_on_given_nodeG(
    // const int dim, double ElemCoords[], double CanPos[], double dphi[], int FamilyType)=0;
// virtual    void get_dphi_arb_node(std::vector<double> NodeCoord, const int order, double InvJac[], double dphi[])=0;


    //  RETURN SHAPE FUNCTIONS AT GAUSS POINTS --------------------------------------------
    ///< \param[in]  <dim>    Dimension
    ///< \param[in]  <qp>     Gaussian point
    ///< \param[out] <phi>    Shape function

 virtual   void get_phi_gl_g(const int qp, double phi[])=0;
 virtual   void get_phi_gl_g(const int qp, std::vector<double>& phi)=0;
 // virtual   void get_phi_g_arb_el(const int dim, const int qp, double phi[], int FamilyType)=0;
  // virtual  void get_dphi_g_arb_el(const int dim, const int qp, const double InvJac[], double dphi[], int FamilyType)=0;

 //==================================================================================================
    //  Jacobian functions -------------------------------------------------------------------
 virtual   double ComputeInverseMatrix(double Matrix[], double InvMatrix[], int Dimension)=0;
    // Jacobian and inverse Jacobian at gaussian points
 virtual   double JacobianOnGauss(const int dim,const int ng,const double ElemCoords[], double Jac[])=0;    ///< Jacobian (3D-2D)
 virtual   double JacobianOnPoint(const int dim,const double Point[],  double ElemCoords[], double Jac[], int FamilyType, int nShape)=0;
    //   double NablaCoordT2D(const int ng, double x[], double  Jac[])=0;
//   double NablaCoordT3D(const int ng,double xyz[],double Jac[])=0;                   // =====================================
  // virtual  double Jac3D(const int ng, double x[], double InvJac[])=0;                ///< Jacobian (3D)
   virtual double Jac(const int ng,const double x[], double InvJac[])=0;                  ///< Jacobian (3D-2D)
 // virtual   double Jac1D(const int ng, double x[], double InvJac[])=0;                ///< Jacobian (1D)
 // virtual   double Jac2D(const int ng, double x[], double InvJac[])=0;                ///< Jacobian (2D)
  // virtual  double Jac2D(const int ng, double x[], double InvJac[], double Jac[])=0;  ///< Jacobian (2D)
    // Jacobian at gaussian points for arbitrary fem element
// virtual    double JacG(
//     const int ng, double x[], double InvJac[], int FamilyType, int dimension)=0;  ///< Jacobian (3D-2D)
    // Jacobian at nodal points
 virtual   double Jac_nodes(const int ng,const double x[], double InvJac[])=0;    ///< Jacobian (3D-2D)
 // virtual   double Jac1D_nodes(const int ng, double x[], double InvJac[])=0;  ///< Jacobian (1D)
// virtual    double Jac2D_nodes(const int ng, double x[], double InvJac[])=0;  ///< Jacobian (2D)
 // virtual   double Jac3D_nodes(const int ng, double x[], double InvJac[])=0;  ///< Jacobian (3D)
 virtual   double JacOnGivenCanCoords(const int dim,const double ElemCoords[], double CanCoords[], double InvJac[], int FamilyType, int nShape)=0;

    // Jacobian dim -1
 virtual   double JacSur(const int ng,const double x[], double InvJac[]) =0;  ///< Boundary Jacobian (2D-1D)
 // virtual   double JacSur2D(const int ng, double x[], double InvJac[]) const=0;
 // virtual   double JacSur3D(const int n_gauss, double x[], double InvJac[]) const=0;  ///< Surface Jacobian (2D)
    //   double JacSur2D(const int ng,double x[]) const;     ///<  Line Jacobian    (1D)
  // virtual  double JacSur1D(const int ng, double x[], double InvJac[]) const=0;  ///<  P Jacobian    (1D)
  
  //=================================================================================
//--------------------------------------------------------------------------------------
    // functions -----------------------------------------------------------------------
    /// Compute normal normal_g[] at xx[] point
 virtual  void normal_g(const double* xx, double* normal_g) const=0;  ///< unit normal to the surface
 virtual  void normal_g(const double* xx, const double x_c[], double* normal_g) const=0;  ///<  and check x_c[] being interior point
 virtual void normal_g(const double* xx, const double x_c[], double* normal_g, int& sign_normal)  const=0;  ///<  and check x_c[] being interior point
 virtual   void Oxy_face(const double* xx, const double x_c[],  double *normal_g,double* tg_g, double* tg03) const=0;  ///<  and check x_c[] being interior point
    
  //===============================================================================
    // Reading - writing ---------------------------------------------------------------
    // const std::string& file1 << file name to write over
    // const std::string& namefem, << name fem geometry element
    /// Write
virtual    void write(const std::string& filename)=0;
virtual void print_rec_xdmf(const std::string& filename,const std::string& namefem,const std::string& namelem)=0;
//     void print_nodes_xdmf(const std::string& filename,const std::string& namefem)=0;
//     void print_gauss_xdmf(const std::string& filename,const std::string& namefem);
//     void print_function_rec_xdmf(const std::string& filename,const std::string& namefem);
//     void print_function_tetra_xdmf(const std::string& filename,const std::string& namefem);
//     void print_Dfunction_rec_xdmf(const std::string& filename,const std::string& namefem);
//     void print_Dfunction_tetra_xdmf(const std::string& filename,const std::string& namefem);

//======

// ================================================================
// Special functions
// =========================================================
  virtual  double Tj_interp_lin(double *u,int j, double *xi,int ndof,int idim)=0;
  virtual  double dTjdxi_lin(double *u,int j, double *xi,int i,int ndof,int idim)=0;
  virtual  double ComputeInverseMatrix_vol(double Matrix[], double InvMatrix[])=0;
  virtual  double ComputeInverseMatrix_sur(double Matrix[], double InvMatrix[])=0;
  virtual  bool phys_to_ref(double xref[],double xphys[],const double Xel[], int max_it = 30, double tol = 1e-12) = 0;

  // =================================================================================
protected: 
 
  // shape values at points [max shapes=27]
  double    phi[27]; /* Shape function values.  */
  double    dphi[81]; /* Shape function derivative values.  */
  double    dphidxi[27]; /* Shape function derivatives in the xi and eta direction.   */
  double    dphideta[27]; /* Shape function derivatives in the eta direction. (dphi_i/dxi) at point P */
  double    dphidzeta[27]; /* Shape function derivatives in the zeta direction.*/
  
  double    dphidx[27]; /* Shape function derivatives in the x direction.*/
  double    dphidy[27]; /* Shape function derivatives in the y direction.*/
  double    dphidz[27];  /* Shape function derivatives in the z direction.*/

  double    d2phi[135];        /* Shape function second derivative values.  */
  double    d2phidxi2[27];     /* Shape function second derivatives in the xi direction. */
  double    d2phidxideta[27];  /* Shape function second derivatives in the xi-eta direction.  */
  double    d2phidxidzeta[27]; /* Shape function second derivatives in the xi-zeta direction. */
  double    d2phideta2[27];    /* Shape function second derivatives in the eta direction.  */
  double    d2phidetadzeta[27];/* Shape function second derivatives in the eta-zeta direction. */
  double    d2phidzeta2[27];   /* Shape function second derivatives in the zeta direction. */
  
  double    d2phidx2[27];      /* Shape function second derivatives in the x direction.  */
  double    d2phidxdy[27];     /* Shape function second derivatives in the x-y direction.  */
  double    d2phidxdz[27];     /* Shape function second derivatives in the x-z direction.  */
  double    d2phidy2[27];      /* Shape function second derivatives in the y direction.  */
  double    d2phidydz[27];     /* Shape function second derivatives in the y-z direction.  */
  double    d2phidz2[27];      /* Shape function second derivatives in the z direction.  */
 
  double  curl_phi[81];        /* Shape function curl values. Only defined for vector types.   */
  double  div_phi[27];         /* Shape function divergence values. Only defined for vector types.  */
 
public:
 
   //  ----------------------------------------------------------------------------
 virtual void compute_shape_functions(const int  elem, const double  qp[]);
  // virtual computation at gaussian points
 virtual double compute_fe_gauss_points(const int  qp,const double xyz_vol[], double InvJac[])=0; //  compute  vol (xyz_vol[]=nodes with shapes)
 virtual double compute_fe_face_gauss_points(const int  qp, const double  xyz_bd[])=0; //  compute surface (xyz_bd[]=nodes with shapes)
 
  // ----------------------------------------------------------------------------------
  // \returns The shape function values at the quadrature points on the element.  ------------------------------->
  const double & get_phi() const  { /*assert(!_calculations_started || _calculate_phi);*/ /*_calculate_phi = true;*/ return *phi; }
  virtual void request_phi() const override  { get_phi(); }

  //  ----------------------------------------------------------------------------
  // First derivative ------------------------------------------------------------>  
  // \returns The shape function derivatives at the quadrature points.  ---------------------------------------------->
  const double & get_dphi() const  {/* assert(!_calculations_started || _calculate_dphi); *//*_calculate_dphi = _calculate_dphiref = true;*/ return *dphi; }
  virtual void request_dphi() const override  { get_dphi(); }
 // const std::vector<std::vector<OutputGradient>> & get_dual_dphi() const
  // { assert(!_calculations_started || _calculate_dphi);
  //   _calculate_dphi = calculate_dual = _calculate_dphiref = true; return dual_dphi; }
  // \returns The shape function x-derivative at the quadrature points.   */
  const  double & get_dphidx(int di) const {
    /* assert(!_calculations_started || _calculate_dphi); _calculate_dphi = _calculate_dphiref = true; */
      switch(di){ 
            case 0: assert(_dim>0); return *dphidx;
            case 1: assert(_dim>1); return *dphidy;
            case 2: assert(_dim>2); return *dphidz;
            default: abort(); return *dphidx;
          }
    }
  // /* \returns The shape function y-derivative at the quadrature points.   */
  // const double & get_dphidy() const{/* assert(!_calculations_started || _calculate_dphi); _calculate_dphi = _calculate_dphiref = true;*/ return *dphidy; }
  // /* \returns The shape function z-derivative at the quadrature points.  */
  // const double & get_dphidz() const{/*assert(!_calculations_started || _calculate_dphi); _calculate_dphi = _calculate_dphiref = true;*/ return *dphidz; }
    
  const  double & get_dphidxi(int di) const {
    /* assert(!_calculations_started || _calculate_dphi); _calculate_dphi = _calculate_dphiref = true; */
      switch(di){ 
            case 0: assert(_dim>0); return *dphidxi;
            case 1: assert(_dim>1); return *dphideta;
            case 2: assert(_dim>2); return *dphidzeta;
            default: abort(); return *dphidxi;
          }
    }  
    
  // /* \returns The shape function xi-derivative at the quadrature points.  */
  // const double & get_dphidxi() const
  // { /*assert(!_calculations_started || _calculate_dphiref);*//*_calculate_dphiref = true;*/ return *dphidxi; }
  // /* \returns The shape function eta-derivative at the quadrature points.  */
  // const double & get_dphideta() const
  // { /*assert(!_calculations_started || _calculate_dphiref);*//* _calculate_dphiref = true;*/ return *dphideta; }
  // /* \returns The shape function zeta-derivative at the quadrature points.  */
  // const double & get_dphidzeta() const
  // { /*assert(!_calculations_started || _calculate_dphiref);_calculate_dphiref = true; */return *dphidzeta; }

  // Special calculations 
  //  \returns The curl of the shape function at the quadrature points  when  _calculate_dphiref = true; 
  virtual    const  double & get_curl_phi() const  {assert(_calculate_curl_phi);return *curl_phi; };// _calculate_curl_phi = _calculate_dphiref = true;*/ 
  /* \returns The divergence of the shape function at the quadrature points when _calculate_div_phi= true   */
  const double & get_div_phi() const{assert(_calculate_div_phi); return *div_phi; }// = _calculate_dphiref ;

   // ---------------------------------------------------------------------------------------------------
  // Second derivative ------------------------------------------------------------> 
  //  ==================================================================================================== 
const double & get_d2phi() const  { /*assert(!_calculations_started || _calculate_phi);*/ /*_calculate_phi = true;*/ return *d2phi; }
  /* \returns The shape function second derivatives at the quadrature points, at reference coordinates   */
  //  di=0-> dxi di=1-> deta di=2-> dzeta
  const double & get_d2phidxidxj(
    int di,
    int dj
  ) const {  //  ====================================================================================================
      /*assert(!_calculations_started || _calculate_d2phi);
      _calculate_d2phi = _calculate_dphiref = true;*/
      switch(di){
        case 0: {  // di=0-> dx
          switch(dj){ assert(_dim>0);
            case 0: return *d2phidxi2;
            case 1: return *d2phidxideta;
            case 2:  return *d2phidxidzeta;
            default: abort(); return *d2phidxi2;
          }
        }
        case 1: {  // di=1-> dy
          switch(dj){ assert(_dim>1);
            case 0: return  *d2phidxideta;
            case 1: return  *d2phideta2;
            case 2:  return *d2phidetadzeta;
            default: abort(); return *d2phidxi2;
          }
        }
        case 2: {  // di=2-> dz
          switch(dj){ assert(_dim>2);
            case 0: return  *d2phidxidzeta;
            case 1: return  *d2phidetadzeta;
            case 2:  return *d2phidzeta2;
            default: abort(); return *d2phidxi2;
          }
        }
      default: abort(); return *d2phidxi2;
    }
  }
  //  ====================================================================================================
  /* \returns The shape function second derivatives at the quadrature points, at real coordinates   */
  //  di=0-> dx di=1-> dy di=2-> dz 
 const double & get_d2phididj(
   int di,int dj
 ) const {  //  ====================================================================================================
      /*assert(!_calculations_started || _calculate_d2phi);
      _calculate_d2phi = _calculate_dphiref = true;*/
      switch(di){
        case 0: {  // di=0-> dx
          switch(dj){ assert(_dim>0);
            case 0: return *d2phidx2;
            case 1: return *d2phidxdy;
            case 2:  return *d2phidxdz;
            default: abort(); return *d2phidx2;
          }
        }
        case 1: {  // di=1-> dy
          switch(dj){ assert(_dim>1);
            case 0: return  *d2phidxdy;
            case 1: return  *d2phidy2;
            case 2:  return *d2phidydz;
            default: abort(); return *d2phidx2;
          }
        }
        case 2: {  // di=2-> dz
          switch(dj){ assert(_dim>2);
            case 0: return  *d2phidxdz;
            case 1: return  *d2phidydz;
            case 2:  return *d2phidz2;
            default: abort(); return *d2phidx2;
          }
        }
      default: abort(); return *d2phidx2;
    }
  }  
  
//  ====================================================================================================
virtual void attach_quadrature_rule (QBase *q){};

  /* Prints the value of each shape function at each quadrature point.  */
  // void print_phi(std::ostream & os) const override;
  /* Prints the value of each shape function's derivative at each quadrature point.  */
   // void print_dphi(std::ostream & os) const override;
  /* Prints the value of each shape function's second derivatives at each quadrature point.  */
   // void print_d2phi(std::ostream & os) const override;

protected:

  /* Determine which values are to be calculated, for both the FE itself and for the FEMap.   */
  virtual  void determine_calculations();

  /* \returns true iff no calculations have been requested of this FE object or of its associated FEMap  */
  bool calculating_nothing() const  {
    return false;
    // _calculate_nothing &&
    //   !this->_calculate_phi && !this->_calculate_dphi &&
    //   !this->_calculate_d2phi &&
    //   !this->_calculate_curl_phi && !this->_calculate_div_phi &&
    //   !this->_calculate_map;
  }

  /**
   * After having updated the jacobian and the transformation
   * from local to global coordinates in \p FEMap::compute_map(),
   * the first derivatives of the shape functions are
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx, \p dphidy, and \p dphidz. This method
   * should rarely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected.
   */
 
  /* Object that handles computing shape function values, gradients, etc in the physical domain.   */
  // std::unique_ptr<FETransformationBase<OutputType>> _fe_trans;



};



// // Typedefs for convenience and backwards compatibility
// typedef MGFEbase<Real> FEBase;
// typedef MGFEbase<RealGradient> FEVectorBase;




// ------------------------------------------------------------
// MGFEbase class inline members
inline
MGFEbase::MGFEbase(const unsigned int d, const FEType & fet) :
  MGFE0(d,fet),
  // _fe_trans( FETransformationBase<OutputType>::build(fet) ),
  phi(),  dphi(),
  curl_phi(),  div_phi(),
  dphidxi(), dphideta(),  dphidzeta(),  
  dphidx(),  dphidy(),  dphidz(),
  d2phi(),
  d2phidxi2(),  d2phidxideta(),  d2phidxidzeta(),  d2phideta2(),  d2phidetadzeta(),  d2phidzeta2(),
  d2phidx2(),  d2phidxdy(),  d2phidxdz(),  d2phidy2(), d2phidydz(),  d2phidz2()
{}



inline
MGFEbase::~MGFEbase(){}


// #endif // LIBMESH_FE_BASE_H


