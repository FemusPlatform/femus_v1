// c++ libraries ------------------>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

// config includes ------------------>
#include "Printinfo_conf.h"
// #include "Domain_conf.h"


// FEMUS class header
#include "MGFemus_L1.h"
#include "MGFemus_L0_Init.h"
#include "MGSystem_L2.h"
#include "MGSystem_L3_TimeLoop.h"
// #include "MGenum_quadrature_type.h"
#include "MGquadrature_L1_gauss.h"
// Petsc ---------------------->
#ifdef HAVE_PETSCM
#include "petsc.h"  // for Petsc solver
#endif
// Mpi
#ifdef HAVE_MPI
#include <mpi.h>  //For MPI_COMM_WORLD
#endif


// ****************************************************************************
// ****************  Constructor Destructor ***********************************

// // ============================================================================
// // //  This function is the Basic constructor
// FEMUS::FEMUS() : FEMUS0(), _comm(MPI_COMM_WORLD) { // MPI_COMM_NULL or MPI_COMM_WORLD communicator
//   std::abort();
//   // init_femus();
//   return;
// }

// // ============================================================================
// // This function is a constructor with  communicator
// FEMUS::FEMUS(MPI_Comm comm): FEMUS0(), _comm(comm) { // use communicator
//   std::abort();
//   // init_femus();  // Call init class variables
//   return;
// }

// // ================================================================================================
// // This function is a constructor with  fem and mesh classes
// FEMUS::FEMUS(int &argc, char** & argv, MGUtils& mgutils)
//   : FEMUS0(),_comm(MPI_COMM_WORLD) { // MPI_COMM_WORLD communicator _comm(MPI_COMM_WORLD)
// TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::FEMUS(constructor) \n");)
// 
// //  // MPI
//   init_femus(argc, argv);            // Init class variables
// //  // parameters
//   init_param(mgutils);     // Init parameters
// //  // geom + fem
//   init_fem();              // Init finite element
// //  // mesh (extended= med+mpi+hdf5)
//   set_mesh();              // Set mesh
// //  // equation_system
//   init_equation_system();  // Init equation system
// //  // Equation map
//   eqnMap_.reset(new EquationsMap(mgutils._app_dir));
//   mgutils.AddFieldClass(eqnMap_.get());
//   eqnMap_->FillEquationMap(get_MGExtSystem()); // Init the equation class objects and set systems for each class
// //  // construct MG class and operator
//   init_systems();       // Init system data
// 
//    TRACKING_FUN(printf(" <--- FEMUS.C:  FEMUS::FEMUS(constructor) \n");)
//   return;
// }

 // FEMUS(int &argc, char** & argv,  MGUtils& mgutils, std::unique_ptr<MGGeomElBase> &mggeomel,  MGFEMap& mgfemap);
// ================================================================================================
// This function is a constructor with  fem and mesh classes
FEMUS::FEMUS(
  int &argc, char** & argv, 
  MGUtils& mgutils,
  MGGeomMap  & mggeomap,  ///< geom class
  MGFEMap& mgfemap     ///< fem class)
       ):  FEMUS0(),_mg_utils(mgutils), _mg_geomap(mggeomap),_mg_geomel(*mggeomap.get_Geom(0)),_mg_femap(mgfemap), _comm(MPI_COMM_WORLD) { // MPI_COMM_WORLD communicator _comm(MPI_COMM_WORLD)
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::FEMUS(constructor) \n");)

//  // MPI
  init_femus(argc, argv);            // Init class variables
//  // parameters
  init_param(mgutils);     // Init parameters
//  // geom + fem
   // // A) setting MGGeomEl
  
  
  int dim=  _mg_utils._DIMENSION;
  // _mg_geomel= MGGeomElBase::build(dim,LAGRANGE27);
   // if(_mg_geomel == NULL) {std::cout << "FEMUS::init_fem: no _mg_geomel";   abort(); }  //MGeomEl
  // /// B) setting MGFEMap (fem)
 
    // _mg_femap = &mgfemap;
   // if(_mg_femap == NULL) {std::cout << "FEMUS::init_fem: no _mg_femap";  abort(); } //MGFem
  init_fem();              // Init finite element
  // std::unique_ptr<MGGeomElBase> mg_geomel= MGGeomElBase::build(dim,LAGRANGE27);
  // _mg_geomap->set_Geom(mg_geomel.get_Geom(),0);  // quadratic fem
  // init_fem();
//  // mesh (extended= med+mpi+hdf5)
  set_mesh();              // Set mesh
//  // equation_system
  init_equation_system();  // Init equation system
//  // Equation map
  eqnMap_.reset(new EquationsMap(mgutils._app_dir));
  mgutils.AddFieldClass(eqnMap_.get());
  eqnMap_->FillEquationMap(get_MGExtSystem()); // Init the equation class objects and set systems for each class
//  // construct MG class and operator
  init_systems();       // Init system data

   TRACKING_FUN(printf(" <--- FEMUS.C:  FEMUS::FEMUS(constructor) \n");)
  return;
}


// ============================================================================
// This function initializes the class
void FEMUS::init_femus(int &argc, char** & argv) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::init_femus \n");)

  int flag = 0;
  MPI_Initialized(&flag);
  if(flag) _local_MPI_Init = false;       // set _local_MPI_Init
  else _local_MPI_Init = true;            //_local_MPI_Init

  // femus init
  _start.reset(new MGFemusInit(argc, argv));

  _MgEquationMapInitialized = false;
  _MgMeshInitialized = false;

TRACKING_FUN(printf(" <- FEMUS.C:  FEMUS::init_femus \n");)
  return;
}
// =======================================================================
// This function
void FEMUS::init_param(
  MGUtils& mgutils,  ///< utils class
  int name           ///< utils class name
) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::init_param \n");)

  // _mg_utils = &mgutils;
  _mg_utils.set_name(name);  // set ultils class

TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::init_param \n");)
  return;
}

// ================================================================================================
// This init-function sets _mg_geomel and _mg_femap private variable
void FEMUS::init_fem(
  // std::unique_ptr<MGGeomElBase> mggeomel,  ///< geom class
  MGFEMap& mgfemap     ///< fem class
) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::init_fem \n");)
int dim=  _mg_utils._DIMENSION;
  // // A) setting MGGeomEl
 
  // _mg_geomap.reset(mggeomap.release());
   // if(_mg_geomel == NULL) {std::cout << "FEMUS::init_fem: no _mg_geomel";   abort(); }  //MGeomEl
  // /// B) setting MGFEMap (fem)
   // _mg_femap = &mgfemap;
   // if(_mg_femap == NULL) {std::cout << "FEMUS::init_fem: no _mg_femap";  abort(); } //MGFem

TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::init_fem \n");)
  return;
}

// ================================================================================================
void FEMUS::init_fem() {

TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::init_fem \n");)
//   int dim=  _mg_utils._DIMENSION;
//   /// A) setting MGGeomEl
//   // _mg_geomel = MGGeomElBase::build(dim,LAGRANGE27);
//   // if(_mg_geomel == NULL) { std::cout << "FEMUS::init_fem: no _mg_geomel";  abort(); }
//  // _mg_geomap->set_Geom(_mg_geomel.get(),0);  
//   /// B) setting MGFEMap (fem)
//   // _mg_femap = new MGFEMap();
//   
//   QGauss qrule (dim, FIFTH); qrule.init(_mg_geomel.fem_name[0],0);
//   QGauss face_qrule (dim-1, FIFTH); qrule.init(_mg_geomel.fem_name[0],0);
//  
//   int dim_in =  dim; // MGFE<2,LAGRANGE>::MDIM_FEM=fem  dim
//   int shapes_quad;
//   int shapes_lin;
//   int shapes_k;
//   int shapes_linrt;// Fem order input
//   int face_shapes_quad;
//   int face_shapes_lin;
//   int face_shapes_k;
//     // file hdf5
//   std::string file_name_femap=_mg_femap.get_filename_femap();
//    hid_t file_id = H5Fcreate(file_name_femap.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//   /*int face_shapes_linrt*/ // Fem order input
//   int geom_type_in;   // Fem type: 1 27(Quad,Hex) or 0 10(Tri,Tet)
//   int fam_type_in;   // Fem type: 0(Lagrangian) or 1 (Raviart Thomas linear)
//   // group1_id = H5Gcreate(file_id, "/lagrange", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//    FEType *fetq=new  FEType(Order::SECOND,FEFamily::LAGRANGE);
//   
//      
//    FEType *fetl=new  FEType(Order::FIRST,FEFamily::LAGRANGE);
//    FEType *fetk=new  FEType(Order::CONSTANT,FEFamily::LAGRANGE);
//   if(dim_in==3){ // volume MGFE<2,LAGRANGE>::MDIM_FEM=3
//     shapes_quad=27;    shapes_lin=8;    shapes_k=1;    shapes_linrt=6;
//     face_shapes_quad=9;    face_shapes_lin=4;    face_shapes_k=1;   /*face_shapes_linrt=1; */
//     geom_type_in=1;    fam_type_in=0;
//     fetq->set(ElemType::HEX27); fetl->set(ElemType::HEX8);fetk->set(ElemType::NODEELEM);
//   }
//   if(dim_in==2){// surface MGFE<2,LAGRANGE>::MDIM_FEM=2
//     shapes_quad=9;    shapes_lin=4;    shapes_k=1;    shapes_linrt=4;    
//     face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
//     geom_type_in=1;    fam_type_in=0;
//     fetq->set(ElemType::QUAD9); fetl->set(ElemType::QUAD4);fetk->set(ElemType::NODEELEM);
//     //  shapes_quad=6;    shapes_lin=3;    shapes_k=1;    shapes_linrt=0;    
//     // face_shapes_quad=3;    face_shapes_lin=2;    face_shapes_k=1;    /*face_shapes_linrt=1;   */
//     // geom_type_in=0;    fam_type_in=0;
//     // fetq->set(ElemType::TRI6); fetl->set(ElemType::TRI3);fetk->set(ElemType::NODEELEM);
//   }
//   if(dim_in==1){ // line MGFE<2,LAGRANGE>::MDIM_FEM=1
//     shapes_quad=3;   shapes_lin=2;    shapes_k=1;    shapes_linrt=1;
//     face_shapes_quad=1;    face_shapes_lin=1;    face_shapes_k=1;    /*face_shapes_linrt=1; */
//     geom_type_in=1;    fam_type_in=0;
//      fetq->set(ElemType::EDGE3); fetl->set(ElemType::EDGE2);fetk->set(ElemType::NODEELEM);
//   }
//   
//   //volume fem classes -----------------  dim_in FEM ------------------>
//   std::cout<< "DIMENSION= "<< dim_in  <<  " QGauss on FEM ELEM TYPE=" << _mg_geomel.fem_name[0]<< " FIFTH pol \n";
//    // dfe_q_.reset(new MGFE<2,LAGRANGE>(dim_in,2,shapes_quad,geom_type_in,fam_type_in)); //!!!new
//     dfe_q_=MGFEbase::build(dim_in,*fetq);
//     dfe_q_->attach_quadrature_rule (&qrule); // Tell the finite element object to use our quadrature rule.
//     dfe_q_->set_par(file_name_femap,dim_in,2,shapes_quad,geom_type_in,fam_type_in);
//  
//    // dfe_q_->init_qua(); //quad --------------------------
//    // dfe_l_.reset(new MGFE<2,LAGRANGE>(fetl,dim_in,1,shapes_lin,geom_type_in,fam_type_in));
//    // dfe_l_.reset(new MGFE<2,LAGRANGE>(dim_in,1,shapes_lin,geom_type_in,fam_type_in)); //!!!new
//      dfe_l_= MGFEbase::build(dim_in,*fetl);
//       dfe_l_->attach_quadrature_rule (&qrule);
//     dfe_l_->set_par(file_name_femap,dim_in,1,shapes_lin,geom_type_in,fam_type_in);
//    // dfe_l_->init_lin(); // linear -------------------
//    // dfe_k_.reset(new MGFE<2,LAGRANGE>(dim_in,0,shapes_k,geom_type_in,fam_type_in));//!!!new
//   //  dfe_k_.reset(new MGFE<2,LAGRANGE>((fetk,dim_in,0,shapes_k,geom_type_in,fam_type_in));
//    // dfe_k_.reset(new MGFE<2,LAGRANGE>(*fetk));
//      dfe_k_= MGFEbase::build(dim_in,*fetk);
//       dfe_k_->attach_quadrature_rule (&qrule);
//     dfe_k_->set_par(file_name_femap,dim_in,0,shapes_k,geom_type_in,fam_type_in);
//    // dfe_k_->init_pie(); // cosnt -----------------------------------
//   _mg_femap.set_FE(dfe_q_.get(),2);  // quadratic fem
//   _mg_femap.set_FE(dfe_l_.get(),1);  // linear fem
//   _mg_femap.set_FE(dfe_k_.get(),0);  // piecewise fem
//   // <------------------------------------------------------------------
// 
//   // --------------------------  dim_in-1 FEM --------------------------->
//    std::cout<<  "FACE DIMENSION= "<< dim_in-1  << "QGauss on  on FEM ELEM TYPE=" << _mg_geomel.fem_name[1]<< " FIFTH pol \n";
//   if(dim_in-1 > 0) { //face fem classes
//     face_dfe_q_=MGFEbase::build(dim_in-1,*fetq);
//       face_dfe_q_->attach_quadrature_rule (&face_qrule);
//     face_dfe_q_->set_par(file_name_femap,dim_in-1,2,face_shapes_quad,geom_type_in,fam_type_in);
//     // face_dfe_q_.reset(new MGFE<2,LAGRANGE>(*fetq,dim_in-1,2,face_shapes_quad,geom_type_in,fam_type_in));
//     // face_dfe_q_->init_qua(); //quad
//      face_dfe_l_=MGFEbase::build(dim_in-1,*fetl);
//       face_dfe_l_->attach_quadrature_rule (&face_qrule);
//      face_dfe_l_->set_par(file_name_femap,dim_in-1,1,face_shapes_lin,geom_type_in,fam_type_in);
//     // face_dfe_l_.reset(new MGFE<2,LAGRANGE>(*fetl,dim_in-1,1,face_shapes_lin,geom_type_in,fam_type_in));
//     // face_dfe_l_->init_lin(); // linear
//       face_dfe_k_=MGFEbase::build(dim_in-1,*fetk);
//       face_dfe_k_->attach_quadrature_rule (&face_qrule);
//     face_dfe_k_->set_par(file_name_femap,dim_in-1,0,shapes_k,geom_type_in,fam_type_in);
//     // face_dfe_k_.reset(new MGFE<2,LAGRANGE>(*fetk,dim_in-1,0,face_shapes_k,geom_type_in,fam_type_in));
//     // face_dfe_k_->init_pie(); // cosnt
//     _mg_femap.set_FE(face_dfe_q_.get(),5);  // quadratic fem
//     _mg_femap.set_FE(face_dfe_l_.get(),4);  // linear fem
//     _mg_femap.set_FE(face_dfe_k_.get(),3);  // piecewise fem
//   }
  // <------------------------------------------------------------------
  // Raviat-Thomas
//    dfe_rt_=MGFEbase::build(dim_in,*fetl);
//    
//     dfe_rt_->set_par(file_name_femap,dim_in,1,shapes_linrt,geom_type_in,1);
//   // dfe_rt_.reset(new MGFE<2,LAGRANGE>(*fetq,dim_in,1,shapes_linrt,geom_type_in,1));
//   // dfe_rt_->init_lin();
//    _mg_femap.set_FE(dfe_rt_.get(),6); // Raviart-Thomas linear fem

   // status = H5Gclose (group_id);
  // if(_mg_femap == NULL) {    std::cout << "FEMUS::init_fem: no _mg_femap";    abort();  }
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::init_fem \n");)
  return;
}
// ================================================================================================
// This function is the destructor
FEMUS::~FEMUS() {

TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::~FEMUS (destructor****) \n");)

  delete _mg_equations_map;
//  // DO NOT TOUCH ================
  delete _mg_time_loop;
//  //==============================

  delete _mg_mesh;
  // delete _mg_geomel;
  // delete _mg_femap;
//  delete _start;
#ifdef HAVE_MED
//   if(_med_mesh) _med_mesh->decrRef();        // med-mesh
#endif
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::~FEMUS (destructor****) \n");)
  return;
}
// ============================================================================
// This function is the problem destructor
void FEMUS::terminate() {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS:: terminate \n");)
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS:: terminate \n");)
}
// // // ****************************************************************************
// // // ****************    end Constructor Destructor *****************************
// //
// // // ****************************************************************************
// // // ****************    Set    *************************************************
#ifdef TWO_PHASE
void FEMUS::set_mgcc(MGSolCC& cc) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS:: set_mgcc \n");)
  _mg_equations_map->set_mgcc(cc);
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS:: set_mgcc \n");)
  return;
}
#endif
//=============================================================================
// This function sets the controlled domain
void FEMUS::setCtrlDomain(const std::vector<std::vector<double>> coord) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS:: setCtrlDomain \n");)
  _mg_equations_map->eqnmap_ctrl_domain(coord);
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS:: setCtrlDomain \n");)
  return;
}
//=============================================================================
// This function sets the controllable domain
void FEMUS::setControlledDomain(const std::vector<std::vector<double>> coord) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS:: setControlledDomain \n");)
  _mg_equations_map->eqnmap_controlled_domain(coord);
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS:: setControlledDomain \n");)
  return;
}
//=============================================================================
// This function sets the lift domain
// void FEMUS::setLiftDomain(
//   const double xMin, const double xMax, const double yMin, const double yMax, const double zMin,
//   const double zMax) { // ========================================================================
// TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS:: setLiftDomain \n");)
//   _mg_equations_map->eqnmap_lift_domain(xMin, xMax, yMin, yMax, zMin, zMax);
// TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS:: setLiftDomain \n");)
//   return;
// }
void FEMUS::setLiftDomain(const std::vector<std::vector<double>> coord) { 
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS:: setLiftDomain \n");)
  _mg_equations_map->eqnmap_lift_domain(coord);
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS:: setLiftDomain \n");)
  return;
}
//=============================================================================
/// This function returns a value from the systems
double FEMUS::GetValue(const int& ff, int flag) {return _mg_equations_map->GetValue(ff, flag);}
//=============================================================================
/// This function sets a value in the systems
void FEMUS::SetValue(const int& ff, double value) {_mg_equations_map->SetValue(ff, value);}
//=============================================================================
/// This function sets a set of values in the systems
void FEMUS::SetValueVector(const int& ff, std::vector<double> value) {
  _mg_equations_map->SetValueVector(ff, value);
}
// =============================================================================
/// This function sets the type of problem
void FEMUS::init_equation_system(int n_data_points, int n_data_cell) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS:: init_equation_system \n");)
  _mg_equations_map = new EquationSystemsExtendedM(_mg_utils, *_mg_mesh, _mg_femap, n_data_points, n_data_cell);  // MGEquationsMap class
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS:: init_equation_system \n");)  //     return *_mg_equations_map;
  return;
}
// =============================================================================
void FEMUS::init_systems() {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS:: init_systems \n");)

 // check for initialization
 if(_start    == NULL){std::cout<<"FEMUS::init_system: no MGFemusInit* _start ";abort();}
 // if(_mg_utils == NULL){std::cout<<"FEMUS::init_system: no MGUtils* _mg_utils ";abort();}
 // if(_mg_geomel== NULL){std::cout<<"FEMUS::init_system: no MGGeomEl* _mg_geomel ";abort();}
 // if(_mg_femap == NULL){std::cout<<"FEMUS::init_system: no MGFEMap* _mg_femap ";abort();}
 if(_mg_mesh  == NULL){std::cout<<"FEMUS::init_system: no MeshExtended* _mg_mesh ";abort();}
 if(_mg_equations_map == NULL){std::cout<< "FEMUS::init_system: no EquationSystemsExtendedM* _mg_equations_map ";abort();}

   // start function
  _mg_equations_map->init_data(0); // set point data structures
  _mg_equations_map->setDofBcOpIc();  // set operators
  _mg_equations_map->set_mesh_mg(*_mg_mesh);
#ifdef HAVE_MED
  _mg_equations_map->set_mesh_med(*_med_mesh);
#endif
  // time loop
  _mg_time_loop = new MGTimeLoop(_mg_utils, *_mg_equations_map);
  if(_mg_time_loop == NULL) { std::cout << "FEMUS::setSystem: no _mg_time_loop";abort();}
  _MgEquationMapInitialized = true;
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS:: init_systems \n");)

  return;
}
// =============================================================================
// This function sets the mesh from med-mesh (m) to libmesh
void FEMUS::set_mesh() {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS:: set_mesh \n");)

  const int NoLevels = (int)_mg_utils._dict_geom["nolevels"];
  _mg_mesh = new MeshExtended(_start->comm(), _mg_utils, _mg_geomel);
  // check insanity
  if(_mg_mesh == NULL) {std::cout << "FEMUS::setMesh: no _mg_mesh";    abort();}
  if(NoLevels != _mg_mesh->_NoLevels) {
    std::cout << "Inconsistent Number of Levels between Mesh and SolBase" << std::endl;  abort();}
  // print mesh at level NoLevels-1 (linear connectivity)
  _mg_mesh->print(NoLevels - 1, 0);

#ifdef HAVE_MED
  // prind mesh at level NoLevels-1 (med format)
  std::string mesh_name = _mg_utils.get_file("F_MESH_READ");  //://= _mg_utils.get_file("F_MESH_READ");
  unsigned pos = mesh_name.find(".");                          // position of "live" in str
  std::ostringstream name;
  if(_mg_mesh->_iproc == 0) {
    name << _mg_utils._mesh_dir << mesh_name.substr(0, pos) << "_fine.med";
  }
#endif
  _MgMeshInitialized = true;

TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS:: set_mesh \n");)

  return;
}
// *******************************************************************
int FEMUS::get_proc() const { return _mg_mesh->_iproc; }
// *******************************************************************
// *******************************************************************
// **************** Solve  *******************************************
// *******************************************************************
/// This function sets up the intial set
void FEMUS::solve_setup(int& t_in, double& time) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::solve_setup \n"); )

  const int restart = stoi(_mg_utils._dict_config["restart"]);  // restart or not
  // FEMUS:set_time() -> i_in <- initial time iteration; _time <- actual time
  _mg_time_loop->transient_setup(restart,_i_in, _time);          //  MGTimeLoop: setup
  t_in = _i_in;
  time = _time;
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::solve_setup \n"); )

  return;
}
//=============================================================================
// This function solves one step  for transient problems
int FEMUS::solve_onestep(
  const int& eq_min,      ///< eq min to solve --> enum  FIELDS (equations_conf.h)
  const int& eq_max       ///< eq max to solve --> enum  FIELDS (equations_conf.h)
) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::solve_onestep \n");)
  int err = _mg_time_loop->transient_onestep(_i_in,_i_step,_i_print_step,_time,_dt,eq_min,eq_max);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::solve_onestep \n");)

  return err;
}
// ========================================================================
void FEMUS::solve_and_update(
  const int& eq_min,      ///< eq min to solve --> enum  FIELDS (equations_conf.h)
  const int& eq_max       ///< eq max to solve --> enum  FIELDS (equations_conf.h)
) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::solve_and_update \n");)
  _mg_time_loop->transient_solve_and_update(_i_in,_i_step,_i_print_step,_time,_dt,eq_min,eq_max);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::solve_and_update \n");)
  return;
}
// ========================================================================
// This function solves one step  for transient problems
void FEMUS::solve_steady(
  const double& toll,     ///< tolerance
  const int& eq_min,      ///< eq min to solve --> enum  FIELDS (equations_conf.h)
  const int& eq_max       ///< eq max to solve --> enum  FIELDS (equations_conf.h)
) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::solve_steady \n");)
  _mg_time_loop->steady(_n_steps,toll,_i_step,_i_print_step,_dt,eq_min,eq_max);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::solve_steady \n");)
  return;
}
// ========================================================================
// This function solves one step  for transient problems
void FEMUS::set_uooold(
  const int& vec_from,           ///< source vector to be copied
  const int& vec_to,             ///< target vector
  const double& toll,            ///< tolerance
  const double delta_t_step_in,  //   (in)
  const int& eq_min,             ///< eq min to solve --> enum  FIELDS (equations_conf.h) (in)
  const int& eq_max              ///< eq max to solve --> enum  FIELDS (equations_conf.h) (in)
) {
TRACKING_FUN(printf(" --> FEMUS.C: FEMUS::set_uooold \n");)
  _mg_time_loop->set_uooold(vec_from,vec_to,toll,delta_t_step_in,eq_min,eq_max);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::set_uooold \n");)
  return;
}
//=============================================================================
// This function solves one step  for transient problems
void FEMUS::solve_control_onestep(
  const int& nmax_step,   ///< number max of steps         (in)
  const int& eq_min,      ///< eq min to solve --> enum  FIELDS (equations_conf.h) (in)
  const int& eq_max,      ///< eq max to solve --> enum  FIELDS (equations_conf.h) (in)
  bool& converged         ///< check if the solution converged (1->converged)     (out)
) {
TRACKING_FUN(printf(" --> FEMUS.C: FEMUS::solve_control_onestep1 \n");)
  _mg_time_loop->transient_control_onestep(nmax_step,_i_in,_i_step,_i_print_step,_time,_dt,eq_min,eq_max,converged);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C: FEMUS::solve_control_onestep1 \n");)
  return;
}
//=============================================================================
// This function solves one step  for transient problems
void FEMUS::solve_control_onestep(
  const int& max_iter,              ///< max non linear iterations
  const int& eq_min,               ///< eq min to solve --> enum  FIELDS (equations_conf.h) (in)
  const int& eq_max,               ///< eq max to solve --> enum  FIELDS (equations_conf.h) (in)
  std::vector<int> controlled_eq,  ///< equations to solve and to control convergence      (in)
  bool& converged,                 ///< check if the solution converged (1->converged)     (out)
  const double& toll               ///< tolerance
) {
TRACKING_FUN(printf(" --> FEMUS.C: FEMUS::solve_control_onestep2 \n");)
  _mg_time_loop->transient_control_onestep(max_iter,_i_in,_i_step,_i_print_step,_time,_dt,eq_min,eq_max,controlled_eq, converged,toll);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C: FEMUS::solve_control_onestep2 \n");)
  return;
}
// ========================================================================
double FEMUS::System_functional(
  const int& ff,     ///< eqn_number
  double parameter,  ///< functional parameter
  double& control    ///< step control
) {
TRACKING_FUN(printf(" --> FEMUS.C: FEMUS::System_functional \n");)
  return _mg_equations_map->System_functional(ff, parameter, control);
TRACKING_FUN(printf(" <-- FEMUS.C: FEMUS::System_functional \n");)
}
// ========================================================================
double FEMUS::System_functional(
  const int& ff     ///< eqn_number
) {
TRACKING_FUN(printf(" --> FEMUS.C: FEMUS::System_functional \n");)
  return _mg_equations_map->System_functional(ff);
TRACKING_FUN(printf(" <-- FEMUS.C: FEMUS::System_functional \n");)
}
//=============================================================================
// This function solves one step  for transient problems
void FEMUS::dummy_step() {
TRACKING_FUN(printf(" --> FEMUS.C: FEMUS::dummy_step \n");)
  _mg_time_loop->dummy_step(_i_in,_i_step,_i_print_step,_time,_dt);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C: FEMUS::dummy_step \n");)
  return;
}
//=============================================================================
// This function solves one step for under relaxations problems (no loop inside)
void FEMUS::solve_underrelaxed_onestep(
  const int& eq_min,                  ///< eq min to solve --> enum  FIELDS (equations_conf.h) (in)
  const int& eq_max,                  ///< eq max to solve --> enum  FIELDS (equations_conf.h) (in)
  std::vector<double> controlled_eq,  ///< vector containing numbers of controlled equations
  bool& converged,                    ///< check if the solution converged (1->converged)     (out)
  const double& toll                  ///< tolerance
) {
TRACKING_FUN(printf(" --> FEMUS.C: FEMUS::solve_underrelaxed_onestep \n");)
  _mg_time_loop->transient_underrelaxed_onestep(_i_in,_i_step,_i_print_step,_time,_dt,eq_min,eq_max,controlled_eq,converged,toll);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C: FEMUS::solve_underrelaxed_onestep \n");)
  return;
}
//=============================================================================
// Non linear algorithm
//! non linear= solve_non_linear_onestep + update_non_linear_onestep
double  FEMUS::solve_non_linear_onestep(
  double toll,            ///< tolerance
  int max_iter,              ///< max non linear iterations
  const int& eq_min,      ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
  const int& eq_max      ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::solve_and_update \n");)
  double norm= _mg_time_loop->transient_non_linear_solve(_i_in,_i_step,_i_print_step,_time,_dt,toll,max_iter,eq_min,eq_max);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::solve_and_update \n");)
  return norm;
}
// ================================================================================================
//! non linear= solve_non_linear_onestep + update_non_linear_onestep
void FEMUS::update_non_linear_onestep(
  double toll,            ///< tolerance
  int max_iter,              ///< max non linear iterations
  const int& eq_min,      ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
  const int& eq_max      ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
) {
TRACKING_FUN(printf(" --> FEMUS.C:  FEMUS::update_non_linear_onestep \n");)
 _mg_time_loop->transient_update(_i_in, _i_step, _i_print_step, _time, _dt,toll,max_iter,eq_min,eq_max);  ///< step time
TRACKING_FUN(printf(" <-- FEMUS.C:  FEMUS::update_non_linear_onestep \n");)
  return;
}

