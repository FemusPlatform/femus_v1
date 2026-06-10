// std include ------------> 
#include <iostream>
// Femus include --------------->
 #include "MGGeom_L1_El.h"
 
 

 
//Constructor no  use !!!!!!!!!!!!!!!
 // virtual int n_el_shapes(int order, int ib_flag=0)=0;
 //  virtual int n_el_subdomains(int ib_flag=0)=0;
 // virtual  int n_el_sides(int ib_flag=0)=0;
 // virtual std::string el_names(int type=0, int ibflag=0)=0;
MGGeomElBase::MGGeomElBase (const unsigned int dim, const GEOFamily fefam) {
  // n_q[0]=n_el_shapes(2,0);n_q[1]=n_el_shapes(2,1);   ///< Quadratic number of nodes (volume+surface)
  // n_l[0]=n_el_shapes(1,0);n_l[1]=n_el_shapes(1,1);   ///< Linear number of nodes (volume+surface)
  // n_se[0]= n_el_subdomains(0);n_se[1]= n_el_subdomains(1);  ///< Number of linear subelements (volume+surface)
  //  _n_sides[0]=n_el_sides(0); _n_sides[1]=n_el_sides(1);  /// Number of element sides (volume+surface)
  // name[0]=el_names(0,0);name[1]=el_names(0,1);
  // pname[0]=el_names(1,0); pname[1]=el_names(1,1);
  // // std::string name[2];  ///< Element name (volume+surface)
  // // std::string pname[2]; //
  
};
 
 
//  ==========================================================
// This build function allows the use of template 
// class MGGeomEl<unsigned const int DIM, GEOFamily geofam > calling 
// thos father MGGeomElBase(const unsigned int dim,     const GEOFamily geofam)
std::unique_ptr<MGGeomElBase> MGGeomElBase::build (
    const unsigned int dim, 
    const GEOFamily geofam
){//  ==========================================================
  switch (dim){
    case 1:   {  // 1D -------------------------------------
        switch (geofam){
          case LAGRANGE27:
            return std::unique_ptr<MGGeomElBase>(new MGGeomElLag27_1D(1, LAGRANGE27));
            case LAGRANGE10:
             return std::unique_ptr<MGGeomElBase>(new MGGeomElLag10_1D(1, LAGRANGE10));
              //  case LAGRANGE14:
              // return std::unique_ptr<MGGeomElBase>(new MGGeomElLag14_1D(1, LAGRANGE14));
                default:
              std::cout <<" MGGeomEl  family are LAGRANGE27-10-14 "; abort();
          }
      }
      case 2:   { // 2D -------------------------------------
        switch (geofam) {
          case LAGRANGE27:
            return std::unique_ptr<MGGeomElBase>(new MGGeomElLag27_2D(2, LAGRANGE27));
            case LAGRANGE10:
             return std::unique_ptr<MGGeomElBase>(new MGGeomElLag10_2D(2, LAGRANGE10));
           // case LAGRANGE14:
           //  return std::unique_ptr<MGGeomElBase>(new MGGeomElLag14_2D(2, LAGRANGE14)); 
              default:
              std::cout <<" MGGeomEl  family are LAGRANGE27-10-14 "; abort();
          }
      }
      case 3:  {     // 3D -------------------------------------
    
        switch (geofam) {
          case LAGRANGE27:
            return std::unique_ptr<MGGeomElBase>(new MGGeomElLag27_3D(3, LAGRANGE27));
            case LAGRANGE10: 
             return std::unique_ptr<MGGeomElBase>(new MGGeomElLag10_3D(3, LAGRANGE10));
            // case LAGRANGE14: 
            // return std::unique_ptr<MGGeomElBase>(new MGGeomElLag14_3D(3, LAGRANGE14));
            default:
              std::cout <<" MGGeomEl  family are LAGRANGE27-10-14 "; abort();
          }
      }
      
    }// switch
    return nullptr;
}

