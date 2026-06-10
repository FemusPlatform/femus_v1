#ifndef __mggeomel0_h__
#define __mggeomel0_h__

#include <string>
// #include "MGFE_conf.h"
#include <math.h>
// #include "MGFE.h"
#include <memory>
// ====================================================
// Only 3 geometric elements labelled ELTYPE (in #define -> MGFE_conf)
// ELTYPE == 27 (Hexahedron) ELTYPE == 10 (Tetrahedron) ELTYPE == 18 (Prism)
// MGFE_MDIM (in #define -> Domain_conf.h)=  FEM_MaxDIMension 
// 
// =====================================================
#include "MGGeom_L0_enum.h"
#include"MGFE_L0_II_enum_order.h"



// Father class to   template<unsigned int DIM,FEFamily T>MGGeomEl

class MGGeomElBase{
  
  //  // Constructor-Destructor  
  // MGGeomElBase (const unsigned int dim, const GEOFamily fefam);
public:

   /// GEOMETRIC CONFIGURATION OF THE ELEMENT
  int n_q[2];   ///< Quadratic number of nodes (volume+surface)
  int n_l[2];   ///< Linear number of nodes (volume+surface)
  int n_se[2];  ///< Number of linear subelements (volume+surface)
  int _n_sides[2]; /// Number of element sides (volume+surface)
  ElemType  fem_name[3];
  std::string name[3];  ///< Element name (volume+surface)
  std::string pname[3]; /// Print element name (volume+surface) for XDMF print (linear)
  
  
  // Constructor-Destructor  
  MGGeomElBase (const unsigned int dim, const GEOFamily fefam);
  MGGeomElBase () = default;// Constructor
  virtual ~MGGeomElBase() = default; // Destructor.
  // * Builds a specific GeomEl type with
  // dim=DIMENSION -> 1,2,3
  // fefam=FEFamily -> Lagrange27(hehaedral) Largrange14(tetrahedal) Lagrange18(Piramid)   
  static std::unique_ptr<MGGeomElBase> build (unsigned int dim, GEOFamily fefam);
  
    // Multilevel indeces
   virtual const  int get_EL_NFACE() const=0;
   virtual const  int get_EL_SUB() const=0;
   // Volume indeces
    virtual const  int get_EL_NPTQ() const=0;// { return 0;}
    virtual const  int get_EL_NPTL() const=0;// {  return 0;}
    virtual const  int get_EL_NPTC() const=0;// {  return 0;}
   // Face indeces 
    virtual const  int get_ELB_NPTQ() const=0;// {  return 0;}
    virtual const  int get_ELB_NPTL() const=0;// { return 0;}
    virtual const  int get_ELB_NPTC() const=0;// { return 0;}   
  
// -------------------------------------------------------------------------------   
  virtual double Prol(int indx)=0; 
  virtual float _embedding_matrix_q(int indx1,int indx2,int indx3)=0; 
  virtual float _embedding_matrix_l(int indx1,int indx2,int indx3)=0; 
  virtual float _embedding_matrix_f(int indx1,int indx2,int indx3)=0; 
  
  virtual int _surf_top(int indx)=0; 
  virtual int _em_ord(int indx)=0;
  
  
  virtual int n_el_shapes(int order, int ib_flag=0)=0;
   virtual int n_el_subdomains(int ib_flag=0)=0;
  virtual  int n_el_sides(int ib_flag=0)=0;
  virtual std::string el_names(int type=0, int ibflag=0)=0;
 
};



#endif
