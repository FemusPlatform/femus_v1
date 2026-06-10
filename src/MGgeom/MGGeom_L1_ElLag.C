
#include  <iostream>
#include "MGGeom_L0_enum.h"
#include "MGGeom_L0_enum.h"

 #include "MGGeom_L1_El.h"


// int el_lag_nshape(const unsigned int dim,int order, const GEOFamily fefam  ,int ib_flag=0);
// int el_lag_nsub(const unsigned int dim, const GEOFamily fefam  ,int ib_flag=0){
// switch(fefam);
// 
// int el_lag_nside(const unsigned int dim, const GEOFamily fefam  ,int ib_flag=0);
// std::string el_lag_name(const unsigned int dim, const GEOFamily fefam  ,int type=0, int ibflag=0);


// =========================================================
int el_lag_nshape(const unsigned int dim,int order, const GEOFamily fefam,int ib_flag){
switch(fefam){    
    case LAGRANGE27: {
         switch(dim){ // dim ------------------------------------
         case 0:   return 1;   
         case 1: {  
             switch(order){
             case 0:  return 1*(1-ib_flag)+1*ib_flag;       
             case 1:  return 2*(1-ib_flag)+1*ib_flag;                                                                                                                                                                                      
             case 2:  return 3*(1-ib_flag)+1*ib_flag;                        
             default: std::cout << "vol order <3"; abort();
          }
         }
          case 2: {  
              switch(order){ // dim=2
             case 0:     return 1*(1-ib_flag)+2*ib_flag;                    
             case 1:     return 4*(1-ib_flag)+2*ib_flag;                                                                                                                                                                                 
             case 2:     return 9*(1-ib_flag)+3*ib_flag;      
             default: std::cout << "bd order <3"; abort();
          }
          }
          case 3:{  
              switch(order){ // dim=3
             case 0:            return 1*(1-ib_flag)+2*ib_flag;                        
             case 1:            return 8*(1-ib_flag)+4*ib_flag;                                                                                                                                                                                     
             case 2:             return 27*(1-ib_flag)+9*ib_flag;                     
             default: std::cout << "bd order <3"; abort();return -1;
          }
          }
          default:  std::cout << " 0< dim <4"; abort();return -1;
          }
    }
    case LAGRANGE10: {
           switch(dim){ // dim ------------------------------------
         case 0:   return 1;   
         case 1:   switch(order){
             case 0:  return 1*(1-ib_flag)+1*ib_flag;       
             case 1:  return 2*(1-ib_flag)+1*ib_flag;                                                                                                                                                                                      
             case 2:  return 3*(1-ib_flag)+1*ib_flag;                        
             default: std::cout << "vol order <3"; abort();
          }
          case 2:  {
              switch(order){ // dim=2
             case 0:     return 1*(1-ib_flag)+2*ib_flag;                    
             case 1:     return 3*(1-ib_flag)+2*ib_flag;                                                                                                                                                                                 
             case 2:     return 6*(1-ib_flag)+3*ib_flag;      
             default: std::cout << "bd order <3"; abort();
          }
          }
          case 3: {
              switch(order){ // dim=3
             case 0:            return 1*(1-ib_flag)+3*ib_flag;                        
             case 1:            return 4*(1-ib_flag)+3*ib_flag;                                                                                                                                                                                     
             case 2:             return 10*(1-ib_flag)+6*ib_flag;                     
             default: std::cout << "bd order <3"; abort();
              }
          }
          default: std::cout << " 0< dim <4"; abort();
          }
    }
     default: std::cout << "This family element was not implemented"; abort();
//eltype
}
    return 0;
}
// =========================================================

int el_lag_nsub(const unsigned int dim, const GEOFamily fefam,int ib_flag){
switch(fefam){    
    case LAGRANGE27: {
         switch(dim){ // dim ------------------------------------
         case 0:   return 1;   
         case 1:   return 2*(1-ib_flag)+1*ib_flag;                   
         case 2:   return 4*(1-ib_flag)+2*ib_flag;                                                                                                                                                                                  
         case 3:   return 8*(1-ib_flag)+4*ib_flag;                     
         default:  std::cout << " 0< dim <4"; abort();
         }
    }
    case LAGRANGE10: {
         switch(dim){ // dim ------------------------------------
         case 0:   return 1;   
         case 1:   return 2*(1-ib_flag)+1*ib_flag;                   
         case 2:   return 4*(1-ib_flag)+2*ib_flag;                                                                                                                                                                                  
         case 3:   return 8*(1-ib_flag)+4*ib_flag;                     
         default: std::cout << " 0< dim <4"; abort();
         }
    }       
     default: std::cout << "This family element was not implemented"; abort();
//eltype
}
    return 0;
}

// =========================================================
int el_lag_nside(const unsigned int dim, const GEOFamily fefam  ,int ib_flag=0){
switch(fefam){    
    case LAGRANGE27: {
         switch(dim){ // dim ------------------------------------
         case 0:   return 1;   
         case 1:   return 2*(1-ib_flag)+0*ib_flag;                   
         case 2:   return 4*(1-ib_flag)+2*ib_flag;                                                                                                                                                                                  
         case 3:   return 6*(1-ib_flag)+4*ib_flag;                     
         default:  std::cout << " 0< dim <4"; abort();
         }
    }
    case LAGRANGE10: {
         switch(dim){ // dim ------------------------------------
         case 0:   return 1;   
         case 1:   return 2*(1-ib_flag)+0*ib_flag;                   
         case 2:   return 3*(1-ib_flag)+2*ib_flag;                                                                                                                                                                                  
         case 3:   return 4*(1-ib_flag)+3*ib_flag;                     
         default:  std::cout << " 0< dim <4"; abort();
         }
    }       
         default: std::cout << "This family element was not implemented"; abort();
//eltype
}
    return 0;
}



    
// =========================================================
std::string el_lag_name(const unsigned int dim, const GEOFamily fefam,int type, int ibflag){
switch(fefam){    
    case LAGRANGE27: {
        switch(type){ // dim --
          case 0: {  //hdf5   
            switch(dim){ // dim ------------------------------------
            case 0:   return "";   
            case 1:  if(ibflag==0) return "Edge_3";  
                     if(ibflag==1) return "Edge_1";
            case 2:  if(ibflag==0) return "Quad_9";  
                     if(ibflag==1) return "Edge_3"; 
            case 3:  if(ibflag==0) return "Hexahedron_27";  
                     if(ibflag==1) return "Quad_9";                     
            default: std::cout << " 0< dim <4"; abort();
            }
          }
          case 1: {  //vtk    
            switch(dim){ // dim ------------------------------------
            case 0:   return ""; 
            case 1:  if(ibflag==0) return "Polyline";  
                     if(ibflag==1) return "Polyvertex"; 
            case 2:  if(ibflag==0) return "Quadrilateral";  
                     if(ibflag==1) return "Polyline";          
            case 3:  if(ibflag==0) return "Hexahedron";  
                     if(ibflag==1) return "Quadrilateral";                                                                                                                                                                                                 
            default:  std::cout << " 0< dim <4"; abort();
            }
          }  
        }
    }
    case LAGRANGE10: {
               switch(type){ // dim --
          case 0: {  //hdf5   
            switch(dim){ // dim ------------------------------------
            case 0:   return "";   
            case 1:  if(ibflag==0) return "Edge_3";  
                     if(ibflag==1) return "Edge_1";
            case 2:  if(ibflag==0) return "Tri_6";  
                     if(ibflag==1) return "Edge_3"; 
            case 3:  if(ibflag==0) return "Tet_10";  
                     if(ibflag==1) return "Tri_6";                     
            default:  std::cout << " 0< dim <4"; abort();
            }
          }
          case 1: {  //vtk    
            switch(dim){ // dim ------------------------------------
            case 0:   return ""; 
            case 1:  if(ibflag==0) return "Polyline";  
                     if(ibflag==1) return "Polyvertex"; 
            case 2:  if(ibflag==0) return "Triangle";  
                     if(ibflag==1) return "Polyline";          
            case 3:  if(ibflag==0) return "Tetrahedron";  
                     if(ibflag==1) return "Triangle";                                                                                                                                                                                                 
            default: std::cout << " 0< dim <4"; abort();
            }
          }  
        }
    }

 
     default: std::cout << "This family element was not implemented"; abort();
//eltype
}
    return 0;
} 
//=======================================
//  n_el_shapes
int MGGeomElLag27_3D::n_el_shapes(int order,int ib_flag  ){ return el_lag_nshape(3,order,LAGRANGE27,ib_flag);}
int MGGeomElLag27_2D::n_el_shapes(int order,int ib_flag  ){ return el_lag_nshape(2,order,LAGRANGE27,ib_flag);}
int MGGeomElLag27_1D::n_el_shapes(int order, int ib_flag  ){ return el_lag_nshape(1,order,LAGRANGE27,ib_flag);}

int MGGeomElLag10_3D::n_el_shapes(int order,int ib_flag){ return el_lag_nshape(3,order,LAGRANGE10,ib_flag);}
int MGGeomElLag10_2D::n_el_shapes(int order,int ib_flag  ){ return el_lag_nshape(2,order,LAGRANGE10,ib_flag);}
int MGGeomElLag10_1D::n_el_shapes(int order, int ib_flag){ return el_lag_nshape(1,order,LAGRANGE10,ib_flag);}
   
 //=======================================
//  subdomains
     int MGGeomElLag27_3D::n_el_subdomains(int ib_flag  ){return el_lag_nsub(3,LAGRANGE27,ib_flag); }
    int MGGeomElLag27_2D::n_el_subdomains(int ib_flag  ) {return el_lag_nsub(2,LAGRANGE27,ib_flag); }
    int MGGeomElLag27_1D::n_el_subdomains(int ib_flag  ) {return el_lag_nsub(1,LAGRANGE27,ib_flag); }  
   
  int MGGeomElLag10_3D::n_el_subdomains(int ib_flag) {return el_lag_nsub(3,LAGRANGE10,ib_flag); }
    int MGGeomElLag10_2D::n_el_subdomains(int ib_flag  ) {return el_lag_nsub(2,LAGRANGE10,ib_flag); }
    int MGGeomElLag10_1D::n_el_subdomains(int ib_flag  ) {return el_lag_nsub(1,LAGRANGE10,ib_flag); }
    
 //=======================================
//  n_el_sides   
    int MGGeomElLag27_3D::n_el_sides(int ib_flag  ){return el_lag_nside(3,LAGRANGE27,ib_flag);}
    int MGGeomElLag27_2D::n_el_sides(int ib_flag  ){return el_lag_nside(2,LAGRANGE27,ib_flag);}
    int MGGeomElLag27_1D::n_el_sides(int ib_flag  ){return el_lag_nside(1,LAGRANGE27,ib_flag);}
   
    int MGGeomElLag10_3D::n_el_sides(int ib_flag  ){return el_lag_nside(3,LAGRANGE10,ib_flag);}
    int MGGeomElLag10_2D::n_el_sides(int ib_flag  ){return el_lag_nside(2,LAGRANGE10,ib_flag);}
    int MGGeomElLag10_1D::n_el_sides(int ib_flag  ){return el_lag_nside(1,LAGRANGE10, ib_flag);}
    
  //=======================================
//  names 

std::string MGGeomElLag27_3D::el_names(int type,int ibflag){return el_lag_name(3, LAGRANGE27,type,ibflag);}
std::string MGGeomElLag27_2D::el_names(int type,int ibflag){return el_lag_name(2, LAGRANGE27,type,ibflag);}
std::string MGGeomElLag27_1D::el_names(int type,int ibflag){return el_lag_name(1, LAGRANGE27,type,ibflag);}
 
std::string MGGeomElLag10_3D::el_names(int type,int ibflag){return el_lag_name(3, LAGRANGE10,type,ibflag);}
std::string MGGeomElLag10_2D::el_names(int type,int ibflag){return el_lag_name(2, LAGRANGE10,type,ibflag);}
std::string MGGeomElLag10_1D::el_names(int type,int ibflag){return el_lag_name(1, LAGRANGE10,type,ibflag);}
 

   
   
