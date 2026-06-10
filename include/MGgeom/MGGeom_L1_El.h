#ifndef __mggeomel_h__
#define __mggeomel_h__

// std lib include ----------->
#include <string>
#include <math.h>
#include <memory>

// Femus lib include -------------------->
#include "MGGeom_L0_El.h"   // Geometric  multi level class basic Level L0


// ====================================================
// Only 3 geometric elements labelled ELTYPE (in #define -> MGFE_conf)
// ELTYPE == 27 (Hexahedron) ELTYPE == 10 (Tetrahedron) ELTYPE == 18 (Prism)
// MGFE_MDIM (in #define -> Domain_conf.h)=  FEM_MaxDIMension 
// 
// =====================================================



/// Class containing all the geometrical informations about the element.
/**Number of nodes and sides, both quadratic and linear, prolongation and embedding matrices
 and surface topology. */

// template<unsigned int DIM>
// class MGGeomElLag27: public MGGeomElBase {
//  public:
// 
//   ///@{ \name CONSTRUCTOR-DESTRUCTOR
//   MGGeomElLag27(const unsigned int dim, const GEOFamily fefam){};
//   ~MGGeomElLag27() = default;
//   ///@}
//  
//   // Multilevel indeces ---------------------------------------------
//    virtual const  int get_EL_NFACE()const=0;//{return EL_NFACE; }
//    virtual const  int get_EL_SUB() const=0;// { return EL_SUB; }
//    // Volume indeces
//     virtual const  int get_EL_NPTQ() const=0;// { return EL_NPTQ;}
//     virtual const  int get_EL_NPTL() const=0;// { return EL_NPTL;}
//     virtual const  int get_EL_NPTC() const=0;// { return EL_NPTC;}
//    // Face indeces 
//     virtual const  int get_ELB_NPTQ() const=0;// { return ELB_NPTQ;}
//     virtual const  int get_ELB_NPTL() const=0;// { return ELB_NPTL;}
//     virtual const  int get_ELB_NPTC() const=0;// { return ELB_NPTC;}  
// 
//   // ----------- embedding matrices// Hex27 family ------------------------------------
//   // HEX 27 family
//   virtual double Prol(int indx)=0;//{return Prol_static[indx];}; 
//   // static const double Prol_static[EL_NPTQ * EL_SUB];
//   virtual float _embedding_matrix_q(int indx1,int indx2,int indx3)=0;//{return embedding_matrix_q[indx1][indx2][indx3];}; 
//   // static const float embedding_matrix_q[EL_SUB][EL_NPTQ][EL_NPTQ];
//   virtual float _embedding_matrix_l(int indx1,int indx2,int indx3)=0;//{return embedding_matrix_l[indx1][indx2][indx3];}; 
//   // static const float embedding_matrix_l[EL_SUB][EL_NPTL][EL_NPTL];
//   virtual float _embedding_matrix_f(int indx1,int indx2,int indx3)=0;//{return embedding_matrix_f[indx1][indx2][indx3];}; 
//   // static const float embedding_matrix_f[EL_SUB][EL_NFACE][EL_NFACE];
//   
//   virtual int _surf_top(int indx)=0;//{return surf_top[indx];}; 
//   // static const int surf_top[EL_NFACE*ELB_NPTQ];  ///< surface topology (side-> nodes)
//   virtual int _em_ord(int indx)=0;//{return em_ord[indx];}; 
//  
// 
// 
// // #elif ELTYPE == 18  // Prism family ====================================>
// // #if FE_DIM == 1
// //   static const double Prol[3 * 2];
// //   static const float _embedding_matrix_q[2][3][3];
// //   static const float _embedding_matrix_l[2][2][2];
// //   static const int _surf_top[2];  ///< surface topology (side-> nodes)
// // #endif
// // #if FE_DIM == 2
// //   static const double Prol[6 * 3];
// //   static const float _embedding_matrix_q[4][6][6];  //
// //   static const float _embedding_matrix_l[4][3][3];  //
// //   static const int _surf_top[9];  ///< surface topology (side-> nodes)
// // #endif
// // #if FE_DIM == 3
// //   static const double Prol[18 * 6];
// //   static const float _embedding_matrix_q[8][18][18];  //
// //   static const float _embedding_matrix_l[8][6][6];    //
// //   static const int _surf_top[39];                     ///< surface topology (side-> nodes)
// // #endif
// // 
// // #endif
// };


//==========================  1D =====================================
class MGGeomElLag27_1D: public MGGeomElBase  {
 public:

  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  MGGeomElLag27_1D(const unsigned int dim, const GEOFamily fefam);
  ~MGGeomElLag27_1D(){};
 
  
 // static constinit std::string name override { "circle" };
  // Multilevel indeces
    const  int get_EL_NFACE()const{return 2; }
    const  int get_EL_SUB() const { return 2; }
   // Volume indeces
    const  int get_EL_NPTQ() const { return 3;}
    const  int get_EL_NPTL() const { return 2;}
    const  int get_EL_NPTC() const { return 2;}
   // Face indeces 
    const  int get_ELB_NPTQ() const { return 2;}
    const  int get_ELB_NPTL() const { return 1;}
    const  int get_ELB_NPTC() const { return 2;}  


  // HEX 27 family
  double Prol(int indx){return Prol_static[indx];} 
  static const double Prol_static[3 *2];
  float _embedding_matrix_q(int indx1,int indx2,int indx3){return embedding_matrix_q[indx1][indx2][indx3];}; 
  static const float embedding_matrix_q[2][3][3];
  float _embedding_matrix_l(int indx1,int indx2,int indx3){return embedding_matrix_l[indx1][indx2][indx3];}; 
  static const float embedding_matrix_l[2][2][2];
  float _embedding_matrix_f(int indx1,int indx2,int indx3){return embedding_matrix_f[indx1][indx2][indx3];}; 
  static const float embedding_matrix_f[2][2][2];
  
  int _surf_top(int indx){return surf_top[indx];}; 
  static const int surf_top[2*1];  ///< surface topology (side-> nodes)
  int _em_ord(int indx){return em_ord[indx];}; 
  static const int em_ord[2];
  
  int n_el_shapes(int order, int ib_flag=0);
  int n_el_subdomains(int ib_flag=0  );
   int n_el_sides(int ib_flag=0);
   std::string el_names(int type,int ibflag=0);
};

//==========================  2D =====================================
/// Class containing all the geometrical informations about the element.
/**Number of nodes and sides, both quadratic and linear, prolongation and embedding matrices
 and surface topology. */

class MGGeomElLag27_2D: public  MGGeomElBase {
 public:

  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  MGGeomElLag27_2D(const unsigned int dim, const GEOFamily fefam);
  ~MGGeomElLag27_2D() = default;
  
  
 // static constinit std::string name override { "circle" };
  // Multilevel indeces
   const  int get_EL_NFACE()const{return 4; }
    const  int get_EL_SUB() const { return 4; }
   // Volume indeces
    const  int get_EL_NPTQ() const { return 9;}
    const  int get_EL_NPTL() const { return 4;}
    const  int get_EL_NPTC() const { return 4;}
   // Face indeces 
   const  int get_ELB_NPTQ() const { return 3;}
    const  int get_ELB_NPTL() const { return 2;}
    const  int get_ELB_NPTC() const { return 4;}  


  // HEX 27 family  2-DIMEN
   double Prol(int indx){return Prol_static[indx];}; 
  static const double Prol_static[9 * 4];
   float _embedding_matrix_q(int indx1,int indx2,int indx3){return embedding_matrix_q[indx1][indx2][indx3];}; 
  static const float embedding_matrix_q[4][9][9];
   float _embedding_matrix_l(int indx1,int indx2,int indx3){return embedding_matrix_l[indx1][indx2][indx3];}; 
  static const float embedding_matrix_l[4][4][4];
   float _embedding_matrix_f(int indx1,int indx2,int indx3){return embedding_matrix_f[indx1][indx2][indx3];}; 
  static const float embedding_matrix_f[4][4][4];
  
  int _surf_top(int indx){return surf_top[indx];}; 
  static const int surf_top[4*3];  ///< surface topology (side-> nodes)
   int _em_ord(int indx){return em_ord[indx];}; 
  static const int em_ord[4];
  
    int n_el_shapes(int order, int ib_flag=0);
int n_el_subdomains(int ib_flag=0  );
 ; int n_el_sides(int ib_flag=0);
 std::string el_names(int type,int ibflag=0);
};
// 
// 
// // ===========================  3 DIM =================================
// /// Class containing all the geometrical informations about the element.
// /**Number of nodes and sides, both quadratic and linear, prolongation and embedding matrices
//  and surface topology. */
// 
class MGGeomElLag27_3D: public  MGGeomElBase {
 public:

  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  MGGeomElLag27_3D(const unsigned int dim, const GEOFamily fefam);
  ~MGGeomElLag27_3D() = default;
  ///@}
 // static constinit std::string name override { "circle" };
  // Multilevel indeces
    const  int get_EL_NFACE()const{return 6; }
    const  int get_EL_SUB() const { return 8; }
   // Volume indeces
    const  int get_EL_NPTQ() const { return 27;}
    const  int get_EL_NPTL() const { return 8;}
    const  int get_EL_NPTC() const { return 6;}
   // Face indeces 
    const  int get_ELB_NPTQ() const { return 9;}
    const  int get_ELB_NPTL() const { return 4;}
    const  int get_ELB_NPTC() const { return 1;}  


  // HEX 27 family
   double Prol(int indx){return Prol_static[indx];}; 
  static const double Prol_static[27 * 8];
   float _embedding_matrix_q(int indx1,int indx2,int indx3){return embedding_matrix_q[indx1][indx2][indx3];}; 
  static const float embedding_matrix_q[8][27][27];
   float _embedding_matrix_l(int indx1,int indx2,int indx3){return embedding_matrix_l[indx1][indx2][indx3];}; 
  static const float embedding_matrix_l[8][8][8];
   float _embedding_matrix_f(int indx1,int indx2,int indx3){return embedding_matrix_f[indx1][indx2][indx3];}; 
  static const float embedding_matrix_f[8][6][6];
  
   int _surf_top(int indx){return surf_top[indx];}; 
  static const int surf_top[6*9];  ///< surface topology (side-> nodes)
  int _em_ord(int indx){return em_ord[indx];}; 
  static const int em_ord[8];
  
  
 int n_el_shapes(int order, int ib_flag=0);
 int n_el_subdomains(int ib_flag=0  );
  int n_el_sides(int ib_flag=0);
  std::string el_names(int type,int ibflag=0);
};


//  ======================================================

// lagrangian 10
//=======================================================

// template<unsigned int DIM>
// class MGGeomElLag10: public MGGeomElBase {
//  public:
// 
//  
//   ///@{ \name CONSTRUCTOR-DESTRUCTOR
//   MGGeomElLag10(const unsigned int dim, const GEOFamily fefam):MGGeomElBase(dim,fefam){};
//    MGGeomElLag10():MGGeomElBase(){};
//   ~MGGeomElLag10() = default;
//   ///@}
// 
//  
// // #if ELTYPE == 10  // Tetra family =====================================>
//    // Multilevel indeces ---------------------------------------------
//    virtual const  int get_EL_NFACE()const=0;//{return EL_NFACE; }
//    virtual const  int get_EL_SUB() const=0;// { return EL_SUB; }
//    // Volume indeces
//     virtual const  int get_EL_NPTQ() const=0;// { return EL_NPTQ;}
//     virtual const  int get_EL_NPTL() const=0;// { return EL_NPTL;}
//     virtual const  int get_EL_NPTC() const=0;// { return EL_NPTC;}
//    // Face indeces 
//     virtual const  int get_ELB_NPTQ() const=0;// { return ELB_NPTQ;}
//     virtual const  int get_ELB_NPTL() const=0;// { return ELB_NPTL;}
//     virtual const  int get_ELB_NPTC() const=0;// { return ELB_NPTC;}  
// 
//   // ----------- embedding matrices// Hex27 family ------------------------------------
//   // HEX 27 family
//   virtual double Prol(int indx)=0;//{return Prol_static[indx];}; 
//   // static const double Prol_static[EL_NPTQ * EL_SUB];
//   virtual float _embedding_matrix_q(int indx1,int indx2,int indx3)=0;//{return embedding_matrix_q[indx1][indx2][indx3];}; 
//   // static const float embedding_matrix_q[EL_SUB][EL_NPTQ][EL_NPTQ];
//   virtual float _embedding_matrix_l(int indx1,int indx2,int indx3)=0;//{return embedding_matrix_l[indx1][indx2][indx3];}; 
//   // static const float embedding_matrix_l[EL_SUB][EL_NPTL][EL_NPTL];
//   virtual float _embedding_matrix_f(int indx1,int indx2,int indx3)=0;//{return embedding_matrix_f[indx1][indx2][indx3];}; 
//   // static const float embedding_matrix_f[EL_SUB][EL_NFACE][EL_NFACE];
//   
//   virtual int _surf_top(int indx)=0;//{return surf_top[indx];}; 
//   // static const int surf_top[EL_NFACE*ELB_NPTQ];  ///< surface topology (side-> nodes)
//   virtual int _em_ord(int indx)=0;//{return em_ord[indx];}; 
//  
//   
//  
// // #if FE_DIM == 1
// //   static const double Prol[3 * 2];
// //   static const float _embedding_matrix_q[2][3][3];
// //   static const float _embedding_matrix_l[2][2][2];
// //   static const int _surf_top[2];  ///< surface topology (side-> nodes)
// // #endif
// // #if FE_DIM == 2
// //   static const double Prol[6 * 3];
// //   static const float _embedding_matrix_q[4][6][6];  // TRI6
// //   static const float _embedding_matrix_l[4][3][3];  // TRI3
// //   static const int _surf_top[9];                    ///< surface topology (side-> nodes)
// // #endif
// // #if FE_DIM == 3
// //   static const double Prol[10 * 4];
// //   static const float _embedding_matrix_q[8][10][10];  // TeT10
// //   static const float _embedding_matrix_l[8][4][4];    // TET4
// //   static const int _surf_top[24];                     ///< surface topology (side-> nodes)
// // #endif
// 
// // #elif ELTYPE == 18  // Prism family ====================================>
// // #if FE_DIM == 1
// //   static const double Prol[3 * 2];
// //   static const float _embedding_matrix_q[2][3][3];
// //   static const float _embedding_matrix_l[2][2][2];
// //   static const int _surf_top[2];  ///< surface topology (side-> nodes)
// // #endif
// // #if FE_DIM == 2
// //   static const double Prol[6 * 3];
// //   static const float _embedding_matrix_q[4][6][6];  //
// //   static const float _embedding_matrix_l[4][3][3];  //
// //   static const int _surf_top[9];  ///< surface topology (side-> nodes)
// // #endif
// // #if FE_DIM == 3
// //   static const double Prol[18 * 6];
// //   static const float _embedding_matrix_q[8][18][18];  //
// //   static const float _embedding_matrix_l[8][6][6];    //
// //   static const int _surf_top[39];                     ///< surface topology (side-> nodes)
// // #endif
// 
// // #endif
// };


// ==========================1 DIM =================================
class MGGeomElLag10_1D: public MGGeomElBase {
 public:

 
  //  CONSTRUCTOR-DESTRUCTOR
  MGGeomElLag10_1D(const unsigned int dim, const GEOFamily fefam);
  ~MGGeomElLag10_1D(){};

   // Basic indeces -------------------------------
   // Multilevel indeces
   const  int get_EL_NFACE()const{return 2; }
   const  int get_EL_SUB() const {return 2; }
   // Volume indeces
   const  int get_EL_NPTQ() const { return 3;}
   const  int get_EL_NPTL() const { return 2;}
   const  int get_EL_NPTC() const { return 2;}
   // Face indeces 
   const  int get_ELB_NPTQ() const { return 1;}
   const  int get_ELB_NPTL() const { return 1;}
   const  int get_ELB_NPTC() const { return 1;} 
    
    
   // embedded matrices ----------------------------- 
   double Prol(int indx){return Prol_static[indx];}
   static const double Prol_static[3 * 2];
   float _embedding_matrix_q(int indx1,int indx2,int indx3){return embedding_matrix_q[indx1][indx2][indx3];} 
   static const float embedding_matrix_q[2][3][3];
   float _embedding_matrix_l(int indx1,int indx2,int indx3){return embedding_matrix_l[indx1][indx2][indx3];} 
   static const float embedding_matrix_l[2][2][2];
   float _embedding_matrix_f(int indx1,int indx2,int indx3){return embedding_matrix_f[indx1][indx2][indx3];}
   static const float embedding_matrix_f[2][2][2];
  
  int _surf_top(int indx){return surf_top[indx];}
  static const int surf_top[2*1];  ///< surface topology (side-> nodes)
   int _em_ord(int indx){return em_ord[indx];}
   static const int em_ord[2]; 
   
 int n_el_shapes(int order,int ib_flag=0);
  int n_el_subdomains(int ib_flag=0  );
   int n_el_sides(int ib_flag=0);
   std::string el_names(int type,int ibflag=0);
};

// // =================== 2 DIM  ========================================
class MGGeomElLag10_2D: public MGGeomElBase {
 public:

  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  MGGeomElLag10_2D(const unsigned int dim, const GEOFamily fefam);
  ~MGGeomElLag10_2D() = default;
  ///@}
 // level --------------------------------------------------
   // Multilevel indeces
     const  int get_EL_NFACE()const{return 3; }
     const  int get_EL_SUB() const { return 4; }
   // Volume indeces
      const  int get_EL_NPTQ() const { return 6;}
      const  int get_EL_NPTL() const { return 3;}
      const  int get_EL_NPTC() const { return 1;}
   // Face indeces 
      const  int get_ELB_NPTQ() const { return 3;}
      const  int get_ELB_NPTL() const { return 2;}
      const  int get_ELB_NPTC() const { return 2;}  
 // #if ELTYPE == 10  // Tetra family --------------------------->
      double Prol(int indx){return Prol_static[indx];}; 
   static const double Prol_static[6 * 3];
    float _embedding_matrix_q(int indx1,int indx2,int indx3){return embedding_matrix_q[indx1][indx2][indx3];}; 
   static const float embedding_matrix_q[4][6][6];
    float _embedding_matrix_l(int indx1,int indx2,int indx3){return embedding_matrix_l[indx1][indx2][indx3];}; 
   static const float embedding_matrix_l[4][3][3];
    float _embedding_matrix_f(int indx1,int indx2,int indx3){return embedding_matrix_f[indx1][indx2][indx3];}; 
   static const float embedding_matrix_f[4][3][3];
  
    int _surf_top(int indx){return surf_top[indx];}; 
   static const int surf_top[3*3];  ///< surface topology (side-> nodes)
    int _em_ord(int indx){return em_ord[indx];}; 
   static const int em_ord[3]; 
   
 int n_el_shapes(int order,int ib_flag=0);
 int n_el_subdomains(int ib_flag=0  );
 int n_el_sides(int ib_flag=0);
 std::string el_names(int type,int ibflag=0);
};

//=======================  3DIM =============================
class MGGeomElLag10_3D: public MGGeomElBase {
 public:

 
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  MGGeomElLag10_3D(const unsigned int dim, const GEOFamily fefam);
  ~MGGeomElLag10_3D() = default;
  ///@}


  
//  static constexpr int EL_NFACE=4; // // ELTYPE 10 Lagrange  faces
//  static constexpr int EL_SUB=8;
//  
//  static constexpr int EL_NPTQ =10;    // ELTYPE 10 Lagrange  quad points (fine)
//  static constexpr int EL_NPTL =4;              // ELTYPE 10 Lagrange lin points (coarse)
//  static constexpr int EL_NPTC =4;           // ELTYPE 10 Lagrange cells (childs)
// 
//  static constexpr int ELB_NPTQ =6; // ELTYPE 10 Lagrange quad boundary points (fine)
//  static constexpr int ELB_NPTL = 3;    // ELTYPE 10 Lagrange linear boundary points (coarse)
//  static constexpr int ELB_NPTC= 1;     // ELTYPE 10 Lagrange  boundary cells  
 
   // Multilevel indeces
   const  int get_EL_NFACE()const{return 4; }
   const  int get_EL_SUB() const { return 8; }
   // Volume indeces
   const  int get_EL_NPTQ() const { return 10;}
   const  int get_EL_NPTL() const { return 4;}
   const  int get_EL_NPTC() const { return 4;}
   // Face indeces 
   const  int get_ELB_NPTQ() const { return 6;}
   const  int get_ELB_NPTL() const { return 3;}
   const  int get_ELB_NPTC() const { return 1;}  
 
  // #if ELTYPE == 10  // Tetra family =====================================>
   double Prol(int indx){return Prol_static[indx];}; 
   static const double Prol_static[10 * 4];
   float _embedding_matrix_q(int indx1,int indx2,int indx3){return embedding_matrix_q[indx1][indx2][indx3];}; 
   static const float embedding_matrix_q[8][10][10];
   float _embedding_matrix_l(int indx1,int indx2,int indx3){return embedding_matrix_l[indx1][indx2][indx3];}; 
   static const float embedding_matrix_l[8][4][4];
   float _embedding_matrix_f(int indx1,int indx2,int indx3){return embedding_matrix_f[indx1][indx2][indx3];}; 
   static const float embedding_matrix_f[8][4][4];
  
   int _surf_top(int indx){return surf_top[indx];}; 
   static const int surf_top[4*6];  ///< surface topology (side-> nodes)
   int _em_ord(int indx){return em_ord[indx];}; 
   static const int em_ord[4];
   
   int n_el_shapes(int order,int ib_flag=0);
   int n_el_subdomains(int ib_flag=0  );
   int n_el_sides(int ib_flag=0);
   std::string el_names(int type,int ibflag=0);
   
};






// lagrangian 27


// template<unsigned int DIM>
// class MGGeomElLag27Lag27:
// public MGGeomElLag27<DIM,LAGRANGE27> {
//  public:
// 
//    
//  
// 
//   ///@{ \name CONSTRUCTOR-DESTRUCTOR
//   MGGeomElLag27Lag27<DIM>(const unsigned int dim, const GEOFamily fefam);
//   ~MGGeomElLag27Lag27() = default;
//   ///@}
// 
//   
//   
//   
// 
//   static constexpr int EL_NFACE=(DIM*2); // // ELTYPE 27 Lagrange  faces
//   static constexpr int EL_SUB=pow(2,DIM);
//  
//  static constexpr int EL_NPTQ =pow(3,DIM);    // ELTYPE 27 Lagrange  quad points (fine)
//  static constexpr int EL_NPTL =pow(2,DIM);    // ELTYPE 27 Lagrange lin points (coarse)
//  static constexpr int EL_NPTC =EL_NFACE;      // ELTYPE 27 Lagrange cells (childs)
// 
//  static constexpr int ELB_NPTQ =(EL_NPTQ+1)/3; // ELTYPE 27 Lagrange quad boundary points (fine)
//  static constexpr int ELB_NPTL = EL_NPTL/2;    // ELTYPE 27 Lagrange linear boundary points (coarse)
//  static constexpr int ELB_NPTC=  EL_NFACE;     // ELTYPE 27 Lagrange  boundary cells  
//   
//  // static constinit std::string name override { "circle" };
//   // Multilevel indeces
//    virtual const  int get_EL_NFACE() const { static const int EL_NFACE {6 }; return EL_NFACE; }
//    virtual const  int get_EL_SUB() const { static const int EL_SUB {9 }; return EL_SUB; }
//    // Volume indeces
//     virtual const  int get_EL_NPTQ() const { static const int EL_NPTQ {27}; return EL_NPTQ;}
//     virtual const  int get_EL_NPTL() const { static const int EL_NPTL { 9}; return EL_NPTL;}
//     virtual const  int get_EL_NPTC() const { static const int EL_NPTC { 6}; return EL_NPTC;}
//    // Face indeces 
//     virtual const  int get_ELB_NPTQ() const { static const int ELB_NPTQ { 9}; return ELB_NPTQ;}
//     virtual const  int get_ELB_NPTL() const { static const int ELB_NPTL { 4}; return ELB_NPTL;}
//     virtual const  int get_ELB_NPTC() const { static const int ELB_NPTC { 6}; return ELB_NPTC;}  
// 
// 
//   virtual double Prol(int indx){return Prol_static[indx];}; 
//   static const double Prol_static[EL_NPTQ * EL_NPTL];
//   virtual float _embedding_matrix_q(int indx1,int indx2,int indx3){return embedding_matrix_q[indx1][indx2][indx3];}; 
//   static const float embedding_matrix_q[EL_NPTL][EL_NPTQ][EL_NPTQ];
//   virtual float _embedding_matrix_l(int indx1,int indx2,int indx3){return embedding_matrix_l[indx1][indx2][indx3];}; 
//   static const float embedding_matrix_l[EL_NPTL][EL_NPTL][EL_NPTL];
//   virtual float _embedding_matrix_f(int indx1,int indx2,int indx3){return embedding_matrix_f[indx1][indx2][indx3];}; 
//   static const float embedding_matrix_f[EL_NPTL][EL_NFACE][EL_NFACE];
//   
//   virtual int _surf_top(int indx){return surf_top[indx];}; 
//   static const int surf_top[EL_NFACE*ELB_NPTQ];  ///< surface topology (side-> nodes)
//   virtual int _em_ord(int indx){return em_ord[indx];}; 
//   static const int em_ord[EL_NPTL];
// 
// 
// };













#endif
