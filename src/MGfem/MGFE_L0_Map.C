#include "MGFE_L0_Map.h"
#include "MGquadrature_L1_gauss.h"




// =====================================================================
 MGFEMap::MGFEMap(
   int dim_in,
   GEOFamily geo_fam_in, 
   Order  qorder_in
 ): dim_(dim_in),
    geo_fam_(geo_fam_in),  // FEFamily::LAGRANGE
    Q_order_(qorder_in) 
 {// ====================================================================
 
 
  char* buff;//automatically cleaned when it exits scope
  // filaname_femap_= std::string(getcwd(buff,255) )+"./../MESH/fem.h5";
  filename_femap_= "./../MESH/fem.h5";
  

  int dim_space =  dim_; // MGFE<2,LAGRANGE>::MDIM_FEM=fem  dim_space
  int shapes_quad;  int shapes_lin;  int shapes_k;  int shapes_linrt;// Fem order input
  int face_shapes_quad;  int face_shapes_lin;  int face_shapes_k;
    // file hdf5
   // hid_t file_id = H5Fcreate(filename_femap_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  /*int face_shapes_linrt*/ // Fem order input
  int geom_type_in;   // Fem type: 1 27(Quad,Hex) or 0 10(Tri,Tet)
  int fam_type_in;   // Fem type: 0(Lagrangian) or 1 (Raviart Thomas linear)
  ElemType elem_type_q; ElemType face_type_q;
  ElemType elem_type_l; ElemType face_type_l;
  ElemType elem_type_k; ElemType face_type_k;
  if(geo_fam_==LAGRANGE27){
    switch(dim_){
      case 1: 
        elem_type_q=EDGE3;    face_type_q=NODEELEM; 
        elem_type_l=EDGE2;    face_type_l=NODEELEM;
        elem_type_k=NODEELEM; face_type_k=NODEELEM;
        shapes_quad=3;   shapes_lin=2;    shapes_k=1;    shapes_linrt=1;
    face_shapes_quad=1;    face_shapes_lin=1;    face_shapes_k=1;    /*face_shapes_linrt=1; */
    geom_type_in=1;    fam_type_in=0;
      break;
      case 2: 
        elem_type_q=QUAD9;    face_type_q=EDGE3;
        elem_type_l=QUAD4;    face_type_l=EDGE2;
        elem_type_k=NODEELEM; face_type_k=NODEELEM;
       shapes_quad=9;    shapes_lin=4;    shapes_k=1;    shapes_linrt=4;    
    face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
    geom_type_in=1;    fam_type_in=0;
      
      break;
      case 3: elem_type_q=HEX27; face_type_q=QUAD9; 
      elem_type_l=HEX8; face_type_l=QUAD4;
      elem_type_k=NODEELEM; face_type_k=NODEELEM;
       shapes_quad=27;    shapes_lin=8;    shapes_k=1;    shapes_linrt=6;
    face_shapes_quad=9;    face_shapes_lin=4;    face_shapes_k=1;   /*face_shapes_linrt=1; */
    geom_type_in=1;    fam_type_in=0;
      break;
      default: std::cout << "Error dim"; exit(2);
     }
    }
      if(geo_fam_==LAGRANGE10){
    switch(dim_){
      case 1: 
        elem_type_q=EDGE3;    face_type_q=NODEELEM;
        elem_type_l=EDGE2;    face_type_l=NODEELEM;
        elem_type_k=NODEELEM; face_type_k=NODEELEM;
       shapes_quad=3;   shapes_lin=2;    shapes_k=1;    shapes_linrt=0;
    face_shapes_quad=1;    face_shapes_lin=1;    face_shapes_k=1;    /*face_shapes_linrt=1; */
    geom_type_in=0;    fam_type_in=0;
      
      break;
      case 2: 
        elem_type_q=TRI6;     face_type_q=EDGE3; 
        elem_type_l=TRI3;     face_type_l=EDGE2;
        elem_type_k=NODEELEM; face_type_k=NODEELEM;
        shapes_quad=6;    shapes_lin=3;    shapes_k=1;    shapes_linrt=0;    
    face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
    geom_type_in=0;    fam_type_in=0;
      break;
      
      case 3: 
        elem_type_q=TET10;    face_type_q=TRI6;
        elem_type_l=TET4;     face_type_l=TRI3;
        elem_type_k=NODEELEM; face_type_k=NODEELEM;
        shapes_quad=10;    shapes_lin=4;    shapes_k=1;    shapes_linrt=0;
    face_shapes_quad=6;    face_shapes_lin=3;    face_shapes_k=1;   /*face_shapes_linrt=1; */
    geom_type_in=0;    fam_type_in=0;
      break;
      default: std::cout << "Error dim"; exit(2);
     }
    }
  // group1_id = H5Gcreate(file_id, "/lagrange", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  

  
   fetq_.reset(new  FEType(Order::SECOND,FEFamily::LAGRANGE));
   fetl_.reset(new  FEType(Order::FIRST,FEFamily::LAGRANGE));
   fetk_.reset(new  FEType(Order::CONSTANT,FEFamily::LAGRANGE));
   
    qrule_.reset(new  QGauss(dim_,  Q_order_));
    qrule_->init(elem_type_q,0); 
  face_qrule_.reset(new  QGauss(dim_-1, Q_order_));
    face_qrule_->init(face_type_q,0); 
//   if(dim_space==3){ // volume MGFE<2,LAGRANGE>::MDIM_FEM=3
//     shapes_quad=27;    shapes_lin=8;    shapes_k=1;    shapes_linrt=6;
//     face_shapes_quad=9;    face_shapes_lin=4;    face_shapes_k=1;   /*face_shapes_linrt=1; */
//     geom_type_in=1;    fam_type_in=0;
//     fetq_->set(ElemType::HEX27); fetl_->set(ElemType::HEX8);fetk_->set(ElemType::NODEELEM);
//   }
//   if(dim_space==2){// surface MGFE<2,LAGRANGE>::MDIM_FEM=2
//     shapes_quad=9;    shapes_lin=4;    shapes_k=1;    shapes_linrt=4;    
//     face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
//     geom_type_in=1;    fam_type_in=0;
//     fetq_->set(ElemType::QUAD9); fetl_->set(ElemType::QUAD4);fetk_->set(ElemType::NODEELEM);
//     //  shapes_quad=6;    shapes_lin=3;    shapes_k=1;    shapes_linrt=0;    
//     // face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
//     // geom_type_in=0;    fam_type_in=0;
//     // fetq_->set(ElemType::TRI6); fetl_->set(ElemType::TRI3);fetk_->set(ElemType::NODEELEM);
//   }
//   if(dim_space==1){ // line MGFE<2,LAGRANGE>::MDIM_FEM=1
//     shapes_quad=3;   shapes_lin=2;    shapes_k=1;    shapes_linrt=1;
//     face_shapes_quad=1;    face_shapes_lin=1;    face_shapes_k=1;    /*face_shapes_linrt=1; */
//     geom_type_in=1;    fam_type_in=0;
//      fetq_->set(ElemType::EDGE3); fetl_->set(ElemType::EDGE2);fetk_->set(ElemType::NODEELEM);
//   }
// }
   
// if(geo_fam_==LAGRANGE10){
//     if(dim_==2){  qrule_->init(TRI6,0);  face_qrule_->init(EDGE3,0);}
//     if(dim_==3){  qrule_->init(TET10,0);  face_qrule_->init(TRI6,0);}
//   
//    fetq_.reset(new  FEType(Order::SECOND,FEFamily::LAGRANGE));
//    fetl_.reset(new  FEType(Order::FIRST,FEFamily::LAGRANGE));
//    fetk_.reset(new  FEType(Order::CONSTANT,FEFamily::LAGRANGE));
//    
//   if(dim_space==3){ // volume MGFE<2,LAGRANGE>::MDIM_FEM=3
//     shapes_quad=10;    shapes_lin=4;    shapes_k=1;    shapes_linrt=0;
//     face_shapes_quad=6;    face_shapes_lin=3;    face_shapes_k=1;   /*face_shapes_linrt=1; */
//     geom_type_in=0;    fam_type_in=0;
//     fetq_->set(ElemType::TET10); fetl_->set(ElemType::TET4);fetk_->set(ElemType::NODEELEM);
//   }
//   if(dim_space==2){// surface MGFE<2,LAGRANGE>::MDIM_FEM=2
//     shapes_quad=6;    shapes_lin=3;    shapes_k=1;    shapes_linrt=0;    
//     face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
//     geom_type_in=0;    fam_type_in=0;
//     fetq_->set(ElemType::TRI6); fetl_->set(ElemType::TRI3);fetk_->set(ElemType::NODEELEM);
//     //  shapes_quad=6;    shapes_lin=3;    shapes_k=1;    shapes_linrt=0;    
//     // face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
//     // geom_type_in=0;    fam_type_in=0;
//     // fetq_->set(ElemType::TRI6); fetl_->set(ElemType::TRI3);fetk_->set(ElemType::NODEELEM);
//   }
//   if(dim_space==1){ // line MGFE<2,LAGRANGE>::MDIM_FEM=1
//     shapes_quad=3;   shapes_lin=2;    shapes_k=1;    shapes_linrt=0;
//     face_shapes_quad=1;    face_shapes_lin=1;    face_shapes_k=1;    /*face_shapes_linrt=1; */
//     geom_type_in=0;    fam_type_in=0;
//      fetq_->set(ElemType::EDGE3); fetl_->set(ElemType::EDGE2);fetk_->set(ElemType::NODEELEM);
//   }
// }
  // --------------------------------------------------------------------------------------------> 
  //volume fem classes -----------------  dim_space FEM ------------------>
  std::cout<< "DIMENSION= "<< dim_space  <<  " QGauss on FEM ELEM TYPE= "<< elem_type_q<< " qrule= "<<Q_order_ <<" \n";
   // dfe_q_.reset(new MGFE<2,LAGRANGE>(dim_space,2,shapes_quad,geom_type_in,fam_type_in)); //!!!new
  fetq_->set(elem_type_q);
  std::unique_ptr<MGFEbase> dfe_q_=MGFEbase::build(dim_space,*fetq_);
  // MGFEbase *   dfe_q_=(MGFEbase::build(dim_space,*fetq_)).get();
    dfe_q_->attach_quadrature_rule (qrule_.get()); // Tell the finite element object to use our quadrature rule.
    dfe_q_->set_par(filename_femap_,dim_space,2,shapes_quad,geom_type_in,fam_type_in);
 
   // dfe_q_->init_qua(); //quad --------------------------
   // dfe_l_.reset(new MGFE<2,LAGRANGE>(fetl_,dim_space,1,shapes_lin,geom_type_in,fam_type_in));
   // dfe_l_.reset(new MGFE<2,LAGRANGE>(dim_space,1,shapes_lin,geom_type_in,fam_type_in)); //!!!new
    fetl_->set(elem_type_l);
    std::unique_ptr<MGFEbase>  dfe_l_= MGFEbase::build(dim_space,*fetl_);
    dfe_l_->attach_quadrature_rule (qrule_.get());
    dfe_l_->set_par(filename_femap_,dim_space,1,shapes_lin,geom_type_in,fam_type_in);
   // dfe_l_->init_lin(); // linear -------------------
   // dfe_k_.reset(new MGFE<2,LAGRANGE>(dim_space,0,shapes_k,geom_type_in,fam_type_in));//!!!new
  //  dfe_k_.reset(new MGFE<2,LAGRANGE>((fetk_,dim_space,0,shapes_k,geom_type_in,fam_type_in));
   // dfe_k_.reset(new MGFE<2,LAGRANGE>(*fetk_));
     fetk_->set(elem_type_k);
     std::unique_ptr<MGFEbase>   dfe_k_= MGFEbase::build(dim_space,*fetk_);
      dfe_k_->attach_quadrature_rule (qrule_.get());
    dfe_k_->set_par(filename_femap_,dim_space,0,shapes_k,geom_type_in,fam_type_in);
   // dfe_k_->init_pie(); // cosnt -----------------------------------
  this->set_FE(dfe_q_,2);  // quadratic fem
  this->set_FE(dfe_l_,1);  // linear fem
  this->set_FE(dfe_k_,0);  // piecewise fem
  // <------------------------------------------------------------------

  // --------------------------  dim_space-1 FEM --------------------------->
   std::cout<<  "FACE DIMENSION= "<< dim_space-1  << "QGauss on  on FEM ELEM TYPE" << elem_type_q<< "qrule"<<Q_order_ <<" \n";
  std::unique_ptr<MGFEbase>   face_dfe_q_, face_dfe_l_, face_dfe_k_;
  // if(dim_space-1 > 0) { //face fem classes
     fetq_->set(face_type_q);
    face_dfe_q_ = MGFEbase::build(dim_space-1,*fetq_);
    face_dfe_q_->attach_quadrature_rule (face_qrule_.get());
    face_dfe_q_->set_par(filename_femap_,dim_space-1,2,face_shapes_quad,geom_type_in,fam_type_in);
    // face_dfe_q_.reset(new MGFE<2,LAGRANGE>(*fetq_,dim_space-1,2,face_shapes_quad,geom_type_in,fam_type_in));
    // face_dfe_q_->init_qua(); //quad
      fetl_->set(face_type_l);
     face_dfe_l_= MGFEbase::build(dim_space-1,*fetl_);
      face_dfe_l_->attach_quadrature_rule (face_qrule_.get());
     face_dfe_l_->set_par(filename_femap_,dim_space-1,1,face_shapes_lin,geom_type_in,fam_type_in);
    // face_dfe_l_.reset(new MGFE<2,LAGRANGE>(*fetl_,dim_space-1,1,face_shapes_lin,geom_type_in,fam_type_in));
    // face_dfe_l_->init_lin(); // linear
      fetk_->set(face_type_k);
     face_dfe_k_=MGFEbase::build(dim_space-1,*fetk_);
    face_dfe_k_->attach_quadrature_rule (face_qrule_.get());
    face_dfe_k_->set_par(filename_femap_,dim_space-1,0,shapes_k,geom_type_in,fam_type_in);
    // face_dfe_k_.reset(new MGFE<2,LAGRANGE>(*fetk_,dim_space-1,0,face_shapes_k,geom_type_in,fam_type_in));
    // face_dfe_k_->init_pie(); // cosnt
    this->set_FE(face_dfe_q_,5);  // quadratic fem
    this->set_FE(face_dfe_l_,4);  // linear fem
    this->set_FE(face_dfe_k_,3);  // piecewise fem
  // }
 
  return;      
}


//  MGFEMap::MGFEMap(){
//         char* buff;//automatically cleaned when it exits scope
//         filaname_femap_= std::string(getcwd(buff,255) )+"/../MESH/fem.h5";
//    
//        int  dim_space=2;
//           QGauss qrule (dim_space, FIFTH); qrule.init(QUAD9,0);
//   QGauss face_qrule (dim_space-1, FIFTH); face_qrule.init(EDGE3,0);
//  
//   int dim_space_in =  dim_space; // MGFE<2,LAGRANGE>::MDIM_FEM=fem  dim_space
//   int shapes_quad;
//   int shapes_lin;
//   int shapes_k;
//   int shapes_linrt;// Fem order input
//   int face_shapes_quad;
//   int face_shapes_lin;
//   int face_shapes_k;
//     // file hdf5
//   // std::string file_name_femap=mg_femap->get_filename_femap();
//    // hid_t file_id = H5Fcreate(file_name_femap.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//   /*int face_shapes_linrt*/ // Fem order input
//   int geom_type_in;   // Fem type: 1 27(Quad,Hex) or 0 10(Tri,Tet)
//   int fam_type_in;   // Fem type: 0(Lagrangian) or 1 (Raviart Thomas linear)
//   // group1_id = H5Gcreate(file_id, "/lagrange", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//    FEType *fetq_=new  FEType(Order::SECOND,FEFamily::LAGRANGE);
//   
//      
//    FEType *fetl_=new  FEType(Order::FIRST,FEFamily::LAGRANGE);
//    FEType *fetk_=new  FEType(Order::CONSTANT,FEFamily::LAGRANGE);
//   if(dim_space_in==3){ // volume MGFE<2,LAGRANGE>::MDIM_FEM=3
//     shapes_quad=27;    shapes_lin=8;    shapes_k=1;    shapes_linrt=6;
//     face_shapes_quad=9;    face_shapes_lin=4;    face_shapes_k=1;   /*face_shapes_linrt=1; */
//     geom_type_in=1;    fam_type_in=0;
//     fetq_->set(ElemType::HEX27); fetl_->set(ElemType::HEX8);fetk_->set(ElemType::NODEELEM);
//   }
//   if(dim_space_in==2){// surface MGFE<2,LAGRANGE>::MDIM_FEM=2
//     shapes_quad=9;    shapes_lin=4;    shapes_k=1;    shapes_linrt=4;    
//     face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
//     geom_type_in=1;    fam_type_in=0;
//     fetq_->set(ElemType::QUAD9); fetl_->set(ElemType::QUAD4);fetk_->set(ElemType::NODEELEM);
//     //  shapes_quad=6;    shapes_lin=3;    shapes_k=1;    shapes_linrt=0;    
//     // face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
//     // geom_type_in=0;    fam_type_in=0;
//     // fetq_->set(ElemType::TRI6); fetl_->set(ElemType::TRI3);fetk_->set(ElemType::NODEELEM);
//   }
//   if(dim_space_in==1){ // line MGFE<2,LAGRANGE>::MDIM_FEM=1
//     shapes_quad=3;   shapes_lin=2;    shapes_k=1;    shapes_linrt=1;
//     face_shapes_quad=1;    face_shapes_lin=1;    face_shapes_k=1;    /*face_shapes_linrt=1; */
//     geom_type_in=1;    fam_type_in=0;
//      fetq_->set(ElemType::EDGE3); fetl_->set(ElemType::EDGE2);fetk_->set(ElemType::NODEELEM);
//   }
//   
//   //volume fem classes -----------------  dim_space_in FEM ------------------>
//   std::cout<< "DIMENSION= "<< dim_space_in  <<  " QGauss on FEM ELEM TYPE= QUAD9 FIFTH pol \n";
//    // dfe_q_.reset(new MGFE<2,LAGRANGE>(dim_space_in,2,shapes_quad,geom_type_in,fam_type_in)); //!!!new
//   std::unique_ptr<MGFEbase> dfe_q_=MGFEbase::build(dim_space,*fetq_);
//   // MGFEbase *   dfe_q_=(MGFEbase::build(dim_space_in,*fetq_)).get();
//     dfe_q_->attach_quadrature_rule (&qrule); // Tell the finite element object to use our quadrature rule.
//     dfe_q_->set_par(filaname_femap_,dim_space_in,2,shapes_quad,geom_type_in,fam_type_in);
//  
//    // dfe_q_->init_qua(); //quad --------------------------
//    // dfe_l_.reset(new MGFE<2,LAGRANGE>(fetl_,dim_space_in,1,shapes_lin,geom_type_in,fam_type_in));
//    // dfe_l_.reset(new MGFE<2,LAGRANGE>(dim_space_in,1,shapes_lin,geom_type_in,fam_type_in)); //!!!new
//      std::unique_ptr<MGFEbase>  dfe_l_= MGFEbase::build(dim_space_in,*fetl_);
//       dfe_l_->attach_quadrature_rule (&qrule);
//     dfe_l_->set_par(filaname_femap_,dim_space_in,1,shapes_lin,geom_type_in,fam_type_in);
//    // dfe_l_->init_lin(); // linear -------------------
//    // dfe_k_.reset(new MGFE<2,LAGRANGE>(dim_space_in,0,shapes_k,geom_type_in,fam_type_in));//!!!new
//   //  dfe_k_.reset(new MGFE<2,LAGRANGE>((fetk_,dim_space_in,0,shapes_k,geom_type_in,fam_type_in));
//    // dfe_k_.reset(new MGFE<2,LAGRANGE>(*fetk_));
//      std::unique_ptr<MGFEbase>   dfe_k_= MGFEbase::build(dim_space_in,*fetk_);
//       dfe_k_->attach_quadrature_rule (&qrule);
//     dfe_k_->set_par(filaname_femap_,dim_space_in,0,shapes_k,geom_type_in,fam_type_in);
//    // dfe_k_->init_pie(); // cosnt -----------------------------------
//   set_FE(dfe_q_.get(),2);  // quadratic fem
//   set_FE(dfe_l_.get(),1);  // linear fem
//   set_FE(dfe_k_.get(),0);  // piecewise fem
//   // <------------------------------------------------------------------
// 
//   // --------------------------  dim_space_in-1 FEM --------------------------->
//    std::cout<<  "FACE DIMENSION= "<< dim_space_in-1  << "QGauss on  on FEM ELEM TYPE= EDGE3 FIFTH pol \n";
//   // if(dim_space_in-1 > 0) { //face fem classes
//      std::unique_ptr<MGFEbase>   face_dfe_q_=MGFEbase::build(dim_space_in-1,*fetq_);
//       face_dfe_q_->attach_quadrature_rule (&face_qrule);
//     face_dfe_q_->set_par(filaname_femap_,dim_space_in-1,2,face_shapes_quad,geom_type_in,fam_type_in);
//     // face_dfe_q_.reset(new MGFE<2,LAGRANGE>(*fetq_,dim_space_in-1,2,face_shapes_quad,geom_type_in,fam_type_in));
//     // face_dfe_q_->init_qua(); //quad
//      std::unique_ptr<MGFEbase>   face_dfe_l_=MGFEbase::build(dim_space_in-1,*fetl_);
//       face_dfe_l_->attach_quadrature_rule (&face_qrule);
//      face_dfe_l_->set_par(filaname_femap_,dim_space_in-1,1,face_shapes_lin,geom_type_in,fam_type_in);
//     // face_dfe_l_.reset(new MGFE<2,LAGRANGE>(*fetl_,dim_space_in-1,1,face_shapes_lin,geom_type_in,fam_type_in));
//     // face_dfe_l_->init_lin(); // linear
//      std::unique_ptr<MGFEbase>  face_dfe_k_=MGFEbase::build(dim_space_in-1,*fetk_);
//     face_dfe_k_->attach_quadrature_rule (&face_qrule);
//     face_dfe_k_->set_par(filaname_femap_,dim_space_in-1,0,shapes_k,geom_type_in,fam_type_in);
//     // face_dfe_k_.reset(new MGFE<2,LAGRANGE>(*fetk_,dim_space_in-1,0,face_shapes_k,geom_type_in,fam_type_in));
//     // face_dfe_k_->init_pie(); // cosnt
//   
//   // }
//     set_FE(face_dfe_q_.get(),5);  // quadratic fem
//     set_FE(face_dfe_l_.get(),4);  // linear fem
//     set_FE(face_dfe_k_.get(),3);  // piecewise fem
//   
// //   int dim_space_in =  dim_space; 
// //   int shapes_quad;  int shapes_lin;  int shapes_k;  int shapes_linrt;// Fem order input
// //   int face_shapes_quad;  int face_shapes_lin;  int face_shapes_k;  /*int face_shapes_linrt*/ // Fem order input
// //   int geom_type_in;   // Fem type: 1 27(Quad,Hex) or 0 10(Tri,Tet)
// //   int fam_type_in;   // Fem type: 0(Lagrangian) or 1 (Raviart Thomas linear)
// //  
// //    FEType *fetq_=new  FEType(Order::SECOND,FEFamily::LAGRANGE);
// //    FEType *fetl_=new  FEType(Order::FIRST,FEFamily::LAGRANGE);
// //    FEType *fetk_=new  FEType(Order::CONSTANT,FEFamily::LAGRANGE);
// //   if(dim_in==3){ // volume MGFE<2,LAGRANGE>::MDIM_FEM=3
// //     shapes_quad=27;    shapes_lin=8;    shapes_k=1;    shapes_linrt=6;
// //     face_shapes_quad=9;    face_shapes_lin=4;    face_shapes_k=1;   /*face_shapes_linrt=1; */
// //     geom_type_in=1;    fam_type_in=0;
// //     fetq_->set(ElemType::HEX27); fetl_->set(ElemType::HEX8);fetk_->set(ElemType::NODEELEM);
// //   }
// //   if(dim_in==2){// surface MGFE<2,LAGRANGE>::MDIM_FEM=2
// //     shapes_quad=9;    shapes_lin=4;    shapes_k=1;    shapes_linrt=4;    
// //     face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
// //     geom_type_in=1;    fam_type_in=0;
// //     fetq_->set(ElemType::QUAD9); fetl_->set(ElemType::QUAD4);fetk_->set(ElemType::NODEELEM);
// //   }
// //   if(dim_in==1){ // line MGFE<2,LAGRANGE>::MDIM_FEM=1
// //     shapes_quad=3;   shapes_lin=2;    shapes_k=1;    shapes_linrt=1;
// //     face_shapes_quad=1;    face_shapes_lin=1;    face_shapes_k=1;    /*face_shapes_linrt=1; */
// //     geom_type_in=1;    fam_type_in=0;
// //      fetq_->set(ElemType::EDGE3); fetl_->set(ElemType::EDGE2);fetk_->set(ElemType::NODEELEM);
// //   }
// //   
// //   //volume fem classes -----------------  dim_in FEM ------------------>
// //    // dfe_q_.reset(new MGFE<2,LAGRANGE>(dim_in,2,shapes_quad,geom_type_in,fam_type_in)); //!!!new
// //    std::unique_ptr<MGFEbase> dfe_q_=MGFEbase::build(dim_in,*fetq_);
// //     dfe_q_->set_par(dim_in,2,shapes_quad,geom_type_in,fam_type_in);
// //    // dfe_q_->init_qua(); //quad --------------------------
// //    // dfe_l_.reset(new MGFE<2,LAGRANGE>(fetl_,dim_in,1,shapes_lin,geom_type_in,fam_type_in));
// //    // dfe_l_.reset(new MGFE<2,LAGRANGE>(dim_in,1,shapes_lin,geom_type_in,fam_type_in)); //!!!new
// //     std::unique_ptr<MGFEbase> dfe_l_= MGFEbase::build(dim_in,*fetl_);
// //     dfe_l_->set_par(dim_in,1,shapes_lin,geom_type_in,fam_type_in);
// //    // dfe_l_->init_lin(); // linear -------------------
// //    // dfe_k_.reset(new MGFE<2,LAGRANGE>(dim_in,0,shapes_k,geom_type_in,fam_type_in));//!!!new
// //   //  dfe_k_.reset(new MGFE<2,LAGRANGE>((fetk_,dim_in,0,shapes_k,geom_type_in,fam_type_in));
// //    // dfe_k_.reset(new MGFE<2,LAGRANGE>(*fetk_));
// //     std::unique_ptr<MGFEbase> dfe_k_= MGFEbase::build(dim_in,*fetk_);
// //     dfe_k_->set_par(dim_in,0,shapes_k,geom_type_in,fam_type_in);
// //    // dfe_k_->init_pie(); // cosnt -----------------------------------
// //   mg_femap->set_FE(dfe_q_.get(),2);  // quadratic fem
// //   mg_femap->set_FE(dfe_l_.get(),1);  // linear fem
// //   mg_femap->set_FE(dfe_k_.get(),0);  // piecewise fem
// //   // <------------------------------------------------------------------
// // 
// //   // --------------------------  dim_in-1 FEM --------------------------->
// //   if(dim_in-1 > 0) { //face fem classes
// //    std::unique_ptr<MGFEbase> face_dfe_q_=MGFEbase::build(dim_in-1,*fetq_);
// //     face_dfe_q_->set_par(dim_in-1,2,face_shapes_quad,geom_type_in,fam_type_in);
// //     // face_dfe_q_.reset(new MGFE<2,LAGRANGE>(*fetq_,dim_in-1,2,face_shapes_quad,geom_type_in,fam_type_in));
// //     // face_dfe_q_->init_qua(); //quad
// //    std::unique_ptr<MGFEbase>  face_dfe_l_=MGFEbase::build(dim_in-1,*fetl_);
// //      face_dfe_l_->set_par(dim_in-1,1,face_shapes_lin,geom_type_in,fam_type_in);
// //     // face_dfe_l_.reset(new MGFE<2,LAGRANGE>(*fetl_,dim_in-1,1,face_shapes_lin,geom_type_in,fam_type_in));
// //     // face_dfe_l_->init_lin(); // linear
// //    std::unique_ptr<MGFEbase>   face_dfe_k_=MGFEbase::build(dim_in-1,*fetk_);
// //     face_dfe_k_->set_par(dim_in-1,0,shapes_k,geom_type_in,fam_type_in);
// //     // face_dfe_k_.reset(new MGFE<2,LAGRANGE>(*fetk_,dim_in-1,0,face_shapes_k,geom_type_in,fam_type_in));
// //     // face_dfe_k_->init_pie(); // cosnt
// //     mg_femap->set_FE(face_dfe_q_.get(),5);  // quadratic fem
// //     mg_femap->set_FE(face_dfe_l_.get(),4);  // linear fem
// //     mg_femap->set_FE(face_dfe_k_.get(),3);  // piecewise fem
// //   }
// //   // <------------------------------------------------------------------
// //   // Raviat-Thomas
// //   // std::unique_ptr<MGFEbase> dfe_rt_=MGFEbase::build(dim_in,*fetl_);
// //   // std::unique_ptr<MGFEbase>  dfe_rt_->set_par(dim_in,1,shapes_linrt,geom_type_in,1);
// //   // // dfe_rt_.reset(new MGFE<2,LAGRANGE>(*fetq_,dim_in,1,shapes_linrt,geom_type_in,1));
// //   // // dfe_rt_->init_lin();
// //   //  mg_femap->set_FE(dfe_rt_.get(),6); // Raviart-Thomas linear fem
//         
//         
//         
//         
//         
//         
//         
//         
//         
//         
//         
//         return;
//   }
