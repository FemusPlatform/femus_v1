// std lib
#include <cmath>
#include <iomanip>
#include <sstream>

// conf includes ------------------
// #include "MGFE_conf.h"       // FEM conf
#include "Printinfo_conf.h"  // print config file
#include "Solverlib_conf.h"  // petsc conf
// #include "Equations_conf.h"
// class includes --------------------->
#include "MGGeom_L1_El.h"
#include "MGMesh_L0.h"
#include "MGUtils_L0.h"  // print filenames
// #include "MeshExtended.h"
// #include "MGSolver_L2_DA.h"
// // #include "parallelM.h"
#ifdef HAVE_PETSCM
#include <mpi.h>  // This is needed in the constructor
#endif



// ===============================================
/// This function manages the printing in Xdmf format
void MGMesh::print(
const int Level,     // level
const int t_step) {
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::print  \n  ");)
// ===============================================================================
  std::string connlin = _mgutils.get_file("CONNLIN");
  const int iproc = _iproc;
  if(iproc == 0) {
    if(connlin=="_conn_lin")     print_conn_lin_hf5(Level);
    else print_conn_quad_hf5(Level);  // print new connection mesh in hdf5 format
    print_vol(Level, t_step);   // print new connection mesh in xdmf format (vol+bd)
  }
// ===================================
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print  \n  ");)

  return;
}

// ==============================================
/// This funcTion prints the volume/boundary
/// Mesh (connectivity) in Xdmf format
void MGMesh::print_vol(
const int Level,       // level
const int /*t_step */  // time
) {                        // =========================================
#ifdef TRACKING_FUN
  printf(" --> MGMesh.C: MGMesh::print_vol  \n  "); // ===================================
#endif
  std::string inout_dir = _mgutils._inout_dir;
  std::string basemesh = _mgutils.get_file("BASEMESH");
  std::string connlin = _mgutils.get_file("CONNLIN");

  std::ostringstream namefile;
  namefile << inout_dir << basemesh << ".xmf";
  std::ostringstream conn_file;
  conn_file /* << inout_dir << "/"*/ << basemesh;
  std::ostringstream topol_file;
  topol_file << /*inout_dir << */ basemesh << connlin << ".h5";
  conn_file << ".h5";

  std::ofstream out(namefile.str().c_str());
  //   int nvrt[2];  std::string mtype[2];
  std::string grid_mesh[2];
  grid_mesh[0] = "Mesh";
  grid_mesh[1] = "Boundary";

  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" \n";
  out << " [ <!ENTITY HeavyData \"mesh.h5 \"> ] ";
  out << ">\n";
  out << " \n";
  out << "<Xdmf> \n"
      << "<Domain> \n";

  for(int ik = 0; ik < _NoFamFEM; ik++) {
    out << "<Grid Name=\"" << grid_mesh[ik].c_str() << "\"> \n";
    // Topology
    print_xmf_topology(out, topol_file.str(), Level, ik);
    // Geometry
    print_xmf_geometry(out, conn_file.str(), Level, ik);
    out << "</Grid> \n";
  }

  out << "</Domain> \n"
      << "</Xdmf> \n";
  out.close();

  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print_vol  \n  ");)
  // ==================================================================
  return;
}
// ==============================================
/// This funcTion prints the Xdmf format mesh topology
void MGMesh::print_xmf_topology(
std::ofstream& out,      //  file xdmf
std::string store_file,  // file where the topology is
int Level,               // Level
int ik                   // Fem family (boundary (1) /volume (0) mesh)
) {                          // =========================================
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::print_xmf_topology  \n  ");)
// ===================================
 int print_xdmf=_mgutils._dict_geom["LINEAR_XDFM"];
  std::cout <<print_xdmf <<"\n";
  if( print_xdmf==1){
  // Topology
  int tot_el =  _GeomEl.n_se[ik]*_NoElements[ik][Level];
  out << "<Topology Type=\"" << _GeomEl.pname[ik] << "\"  Dimensions=\"" << tot_el << "\"> \n";
  out << "<DataStructure DataType=\"Int\" Dimensions=\" " << tot_el << "  " << _GeomEl.n_l[ik]
      << "\" Format=\"HDF\">  \n";
  out << store_file << ":MSH" << ik << "CONN \n";
  out << "</DataStructure> \n"
      << "</Topology> \n";
  }
   else{
  
   // Topology
  int tot_el = /* _GeomEl.n_se[ik]**/_NoElements[ik][Level];
  out << "<Topology Type=\"" << _GeomEl.name[ik] << "\"  Dimensions=\"" << tot_el << "\"> \n";
  out << "<DataStructure DataType=\"Int\" Dimensions=\" " << tot_el << "  " << _GeomEl.n_q[ik]
      << "\" Format=\"HDF\">  \n";
      
  out << "mesh.msh1.h5:/ELEMS"<< _GeomEl.n_q[ik]<<"/CONN/CONN_L"<< _NoLevels-1<<" \n";
  out << "</DataStructure> \n"
      << "</Topology> \n";
}  
    
    
    
  
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print_xmf_topology  \n  ");)
// ===================================

  return;
}

// ==============================================
/// This function prints mesh Xdmf format geometry
void MGMesh::print_xmf_geometry(
std::ofstream& out,      //  file xdmf
std::string store_file,  // file where the geometry is
int Level,               // Level
int /*ik*/               // Fem family (boundary (1) /volume (0) mesh)
) {      // =========================================
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::print_xmf_geometry  \n  ");)
// ===================================

  // Geometry
  out << "<Geometry Type=\"X_Y_Z\"> \n";

  for(int ix = 1; ix < 4; ix++) {
    out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << _NoNodes[Level] << "  " << 1
        << "\" Format=\"HDF\">  \n";
    out << store_file << ":/NODES/COORD/X" << ix << "\n";
    out << "</DataStructure> \n";
  }
  out << " </Geometry>\n";
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print_xmf_geometry  \n  ");)
// ===================================

  return;
}

// ===============================================
//           HDF5 FUNCTION
// ===============================================

void MGMesh::readf_hdf5_NPOINTS(std::string file_name, std::string setdata_name, int /*dim_data*/,int data[]) {
// TRACK( printf(" --> MGMesh.C: MGMesh::readf_hdf5_NPOINTS  \n  "); )
  //      std::string    inout_dir = _mgutils._inout_dir;
  //       std::string    basemesh = _mgutils.get_file ( "BASEMESH" );
  //       const int ndigits = stoi ( _mgutils._dict_config["ndigits"] );
  //       std::ostringstream meshname;
  //       meshname << inout_dir << basemesh << ".h5";

  //       hid_t  file_id     = H5Fopen ( meshname.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t status = 0;
  //       int  topdata[2];
  status = H5Dread(H5Dopen(file_id, setdata_name.c_str(), H5P_DEFAULT),
                   H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  P_INFO(std::cout << "\n Readf_NPOINTS:  Reading mesh from= " << file_name << std::endl;  if(status != 0) abort();)

    H5Sclose(file_id);
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::readf_hdf5_NPOINTS  \n  ");)
  return;
}

/// Read mesh from hdf5 file (namefile)
///           as MGMesh class (MGMesh.h):
/// data parameters : DFL,_type_FEM
/// Verteces        : _NoNodes,_xyz
/// Connectivity    : _NoElements,_conn_map

void MGMesh::read_c() {  // ======================
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::read_c  \n  ");)
  // ===================================

 const int el_nnodes = _GeomEl.n_q[0];  // element nodes
 const int el_nnodes_bd = _GeomEl.n_q[1];  // element nodes
  // Reading parameters from  GeometrySetting.in ---------------------------------------------------
  int err=_mgutils._dict_geom.size();
  double theta = _mgutils._dict_geom["theta"];
  double Coxyz[3];
  Coxyz[0]=  _mgutils._dict_geom["Ox"];
  Coxyz[1]=  _mgutils._dict_geom["Oy"];
  Coxyz[2]=  _mgutils._dict_geom["Oz"];
  // label[0,1,2] defined by  VART0,1,2 
  int type_nonodes[3];
  type_nonodes[0]=(int)(_mgutils._dict_geom["VART0"]);
  type_nonodes[1]=(int)(_mgutils._dict_geom["VART1"]);
  type_nonodes[2]=(int)(_mgutils._dict_geom["VART2"]);
  std::string label[3]; std::string dict_label[4]={"Q","V","F","C",};
 for(int ktype=0;ktype<3;ktype++) label[ktype]= dict_label[type_nonodes[ktype]];
 // restart
 const int restart = stoi(_mgutils._dict_config["restart"]);
if(_mgutils._dict_geom.size()-err) {
  std::cout<< "MGMesh::MGMesh: Missing parameter in GeometrySetting.in"; std::exit(EXIT_FAILURE);}
  //  ---------------------------------------------------
   
//   int ns_type= stoi(_mgutils._dict_config["MG_NavierStokes"]);
//   set_mesh_attribute(1,"F");// 
//   std::string _label[3]={"Q","V","C"};  if(ns_type==3) _label[1]="F";
  // restart from GeometricSetting.in (Mesh:Mesh())
  if(restart == 0) {   //=================================
    std::string input_dir = _mgutils._inout_dir;
    std::string basemesh = _mgutils.get_file("BASEMESH");

    // Open an existing file. ---------------
    std::ostringstream meshname;
    meshname << input_dir << basemesh << ".h5"; 
    std::cout << " Reading mesh from= " << meshname.str() << std::endl; 
    hid_t file_id = H5Fopen(meshname.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t status = 0;
    // Reading DFL -------------------------
    int topdata[4];
    status = H5Dread(H5Dopen(file_id, "/DFLS", H5P_DEFAULT), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, topdata);
    _dim = topdata[0];       // dimension
    _NoFamFEM = topdata[1];  // Volume/Boundary
//     _NoFamFEM = 2;  // Volume/Boundary
    _NoLevels = topdata[2];  // number of levels
    _n_subdom = topdata[3];  // number of subdomains

    // Reading _type_FEM ------------
//     int n_meshes = _NoFamFEM * _NoLevels;
    _type_FEM = new int[_NoFamFEM];
    status = H5Dread(H5Dopen(file_id, "/GEOM_FAM", H5P_DEFAULT), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _type_FEM);
    P_INFO(std::cout << " NoFamFEM= " << _NoFamFEM << " NoLevels= " << _NoLevels << std::endl;)

    // Reading _NoNodes ------------------
    _NoNodes = new int[4*_NoLevels];
    status = H5Dread(H5Dopen(file_id, "/NODES/MAP/SIZE_QVFC", H5P_DEFAULT), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_NoNodes);
    
    // Reading  _xyz ----------------
    double ILref = 1.;
    int n_nodes = _NoNodes[_NoLevels - 1];
    _xyz = new double[_dim * n_nodes];   // double scale[3]
    _xyzo = new double[_dim * n_nodes];  // double scale[3]
    _dxdydz = new double[_dim * n_nodes];
    
   for(int kdim = 0; kdim < _dim ; kdim++) {
    for(int inode = 0; inode <  n_nodes; inode++) {
      _xyz[inode+kdim*n_nodes] = Coxyz[kdim];     // nondimensionalization
      _xyzo[inode+kdim*n_nodes] = Coxyz[kdim];    // nondimensionalization
      _dxdydz[inode+kdim*n_nodes] = Coxyz[kdim];  // nondimensionalization
    }
   }

    std::cout << "-----------------------------------------------\n";
    std::cout << "-----------------------------------------------\n";
    std::cout << "        DIMENSION   " << _dim << std::endl;
    std::cout << "-----------------------------------------------\n";
    std::cout << "-----------------------------------------------\n";

    double* coord;  coord = new double[n_nodes];
  
    for(int kc = 0; kc <2; kc++) {
      std::ostringstream Name;
      Name << "NODES/COORD/X" << kc + 1;
      status = H5Dread(H5Dopen(file_id, Name.str().c_str(), H5P_DEFAULT), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord);
//  double _theta[2]={0.,0.}; double Coxyz[3]={0.,0.,0.};
      
      for(int inode = 0; inode < n_nodes; inode++) {
        _xyz[inode+0*n_nodes]  +=((1-kc)*cos(theta)-(kc)*sin(theta))*coord[inode]*ILref;
        _xyzo[inode+0*n_nodes] +=((1-kc)*cos(theta)-(kc)*sin(theta))*coord[inode]*ILref;
        _xyz[inode+1*n_nodes]  +=((1-kc)*sin(theta)+(kc)*cos(theta))*coord[inode]*ILref;
        _xyzo[inode+1*n_nodes] +=((1-kc)*sin(theta)+(kc)*cos(theta))*coord[inode]*ILref;

        _dxdydz[inode+kc*n_nodes]=0.; //nondimensionalization

//         _xyz[inode + kc * n_nodes] = coord[inode] * ILref;   // nondimensionalization
//         _xyzo[inode + kc * n_nodes] = coord[inode] * ILref;  // nondimensionalization
//         _dxdydz[inode + kc * n_nodes] = 0.;                  // nondimensionalization
      }
    }
      if(_dim==3) {
      std::ostringstream Name;
      Name << "NODES/COORD/X3";
      status = H5Dread(H5Dopen(file_id, Name.str().c_str(), H5P_DEFAULT), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord);
      for(int inode = 0; inode < n_nodes; inode++) {
        _xyz[inode+2*n_nodes]  =coord[inode]*ILref+Coxyz[2];
        _xyzo[inode+2*n_nodes] = coord[inode]*ILref+Coxyz[2];
        _dxdydz[inode+2*n_nodes]=0.; //nondimensionalization
      }
    }
    

    delete[] coord;
    std::cout << " Reading Multimesh with  " << _NoFamFEM << " familes, " << n_nodes << " nodes and "
              << _dim * n_nodes << " coordinates" << std::endl;


    for(int type=0; type<3; type++) {
      _max_nd[type] = new int[_n_subdom * _NoLevels];
      _min_nd[type] = new int[_n_subdom * _NoLevels];

      status = H5Dread(H5Dopen(file_id, ("/NODES/MAP/MIN"+label[type]+"_SL").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_min_nd[type]);
      status = H5Dread(H5Dopen(file_id, ("/NODES/MAP/MAX"+label[type]+"_SL").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_max_nd[type]);
    }

    // Reading  _NoElements --------------
    _off_el = new int*[_NoFamFEM];      // _NoFamFEM (volume)
    _NoElements = new int*[_NoFamFEM];  //  Number of elements (for each level)

    for(int ifem = 0; ifem < _NoFamFEM; ifem++) {
      _off_el[ifem] = new int[_n_subdom * _NoLevels + 1];  // offset for VOLUME
      std::string family=std::to_string(_GeomEl.n_q[0]);
      status = H5Dread(H5Dopen(file_id,("/ELEMS"+family+"/OFFSET_SL").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_el[ifem]);
      _NoElements[ifem] = new int[_NoLevels];
      status = H5Dread(H5Dopen(file_id,("/ELEMS"+family+"/SIZE_L").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_NoElements[ifem]);
    }

    // node mapping ------------------------------------------------------------
//     int ns_type= stoi(_mgutils._dict_config["MG_NavierStokes"]);
    _node_map = new int*[3*_NoLevels];        // 3 node maps; _NoLevels+1 inverse node maps
    int n_nodes_top = _NoNodes[_NoLevels - 1];  //  same length in file (n_nodes_top)
    int* temp = new int[n_nodes_top];
    // temp to read
    for(int type=0; type<3; type++) {
      for(int ilev = 0; ilev < _NoLevels; ilev++) {   // maps from o to 2*_NoLevels-1
//         int type_nonodes = (type == 2) ? 3 : type;
//         if(ns_type==3) type_nonodes = (type == 1) ? 2 : type;
        _node_map[ilev + _NoLevels*type] = new int[_NoNodes[ilev + _NoLevels*type_nonodes[type]]];    // inverse node map
        for(int inode = 0; inode < _NoNodes[ilev + _NoLevels*type_nonodes[type]]; inode++)  _node_map[ilev + _NoLevels*type][inode] = -1;
        // hdf5 reading
        status = H5Dread(H5Dopen(file_id,("/NODES/MAP/MAP" + label[type] + "_L"+std::to_string(ilev)).c_str(), H5P_DEFAULT),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
        // inverse map
        for(int inode = 0; inode < n_nodes_top; inode++)  if(temp[inode] != -1) _node_map[ilev + _NoLevels*type][temp[inode]] = inode;
      }
    }
    delete[] temp;
    // 0 1 2 3 4 5 6 7 8 9 10 11 16 17 18 19 12 13 14 15 24 22 21 23 20 25 26
    // linear coarse map (storage in [2*_NoLevels])
//     _node_map[_NoLevels] = new int[_NoNodes[_NoLevels]];// first linear block
    // hdf5 reading
//     status = H5Dread(   H5Dopen(file_id, ("/NODES/MAP/MAPV_L"+std::to_string(0)).c_str(), H5P_DEFAULT),
//                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);  // read from file
//
//     // inverse map
//     for (int inode = 0; inode < n_nodes_top; inode++) if (temp[inode] != -1) _node_map[_NoLevels][temp[inode]] = inode;

    // connectivity map -----------------------------------------------------
    _el_map = new int*[_NoFamFEM];
    _el_neighbor = new int*[_NoFamFEM];// double tmppp[1125];
    for(int ifem = 0; ifem < _NoFamFEM; ifem++) {
      _el_map[ifem] = new int[_off_el[ifem][_NoLevels * _n_subdom] * el_nnodes];  // Volume

      std::string family=std::to_string(_GeomEl.n_q[0]);
      status = H5Dread(H5Dopen(file_id, ("/ELEMS"+family+"/CONN/CONN_ALL").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_el_map[ifem]);

      // neighbour element map -----------------------------------------------------
      const int n_sides=_GeomEl._n_sides[0];
      _el_neighbor[ifem] = new int[_off_el[ifem][_NoLevels * _n_subdom] * n_sides];  // Volume
      status = H5Dread(H5Dopen(file_id, ("/ELEMS"+family+"/EL_NEIG").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_neighbor[ifem]);
    }
    // Close the file. ---------------------

    if(status < 0) {      std::cout << " File Mesh input in data_in is missing";      abort();    }
    H5Fclose(file_id); 
  }

  else {  //======restart != 0 ======================
    std::string inout_dir = _mgutils._inout_dir;    std::string basemesh = _mgutils.get_file("BASEMESH");
    const int ndigits = stoi(_mgutils._dict_config["ndigits"]);
    std::ostringstream meshname_xyz;
    meshname_xyz << inout_dir << basemesh << "." << std::setw(ndigits) << std::setfill('0') << 0 << ".h5";
    std::ostringstream meshname;    meshname << inout_dir << basemesh << ".h5";
    std::cout << " Reading mesh from= " << meshname.str() << std::endl;
    std::cout << " Reading mesh coordinates from = " << meshname_xyz.str() << std::endl;
    hid_t file_id = H5Fopen(meshname.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t file_id_xyz = H5Fopen(meshname_xyz.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t status = 0;

    // Reading DFL -------------------------
    int topdata[4];
    status = H5Dread(H5Dopen(file_id, "/DFLS",H5P_DEFAULT),H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, topdata);
    _dim = topdata[0];       // dimension
    _NoFamFEM = topdata[1];  // Volume/Boundary
    _NoLevels = topdata[2];  // number of levels
    _n_subdom = topdata[3];  // number of subdomains
    
    // Reading _type_FEM ------------
    int n_meshes = _NoFamFEM * _NoLevels;
    _type_FEM = new int[n_meshes];
//  status = H5Dread(H5Dopen(file_id, "/FEM", H5P_DEFAULT),H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _type_FEM);
    status = H5Dread(H5Dopen(file_id, "/GEOM_FAM", H5P_DEFAULT), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _type_FEM);
    P_INFO(std::cout << " NoFamFEM= " << _NoFamFEM << " NoLevels= " << _NoLevels << std::endl;)
    // Reading _NoNodes ------------------
//  _NoNodes = new int[n_meshes + 1];
    _NoNodes = new int[4*_NoLevels];
//  status = H5Dread(H5Dopen(file_id, "/NODES/MAP/NDxLEV",H5P_DEFAULT),H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _NoNodes);
    status = H5Dread(H5Dopen(file_id, "/NODES/MAP/SIZE_QVFC", H5P_DEFAULT), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_NoNodes);

    // Reading  _xyz ----------------
    double ILref = 1.;
    int n_nodes = _NoNodes[_NoLevels - 1];
    _xyz = new double[_dim * n_nodes];   // double scale[3]
    _xyzo = new double[_dim * n_nodes];  // double scale[3]
    _dxdydz = new double[_dim * n_nodes];
    
    double* coord;   coord = new double[n_nodes];

    for(int kc = 0; kc < _dim; kc++) {
      std::ostringstream Name;     
      Name << "NODES/COORD/X" << kc + 1;
      status = H5Dread(H5Dopen(file_id_xyz, Name.str().c_str(),H5P_DEFAULT),H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,coord);

      for(int inode = 0; inode < n_nodes; inode++) {
        _xyz[inode + kc * n_nodes] = coord[inode] * ILref;   // nondimensionalization
        _xyzo[inode + kc * n_nodes] = coord[inode] * ILref;  // nondimensionalization
        _dxdydz[inode + kc * n_nodes] = 0.;                  // nondimensionalization
      }
    }

    delete[] coord;
#ifdef PRINT_INFO
    std::cout << " Reading Multimesh with  " << n_meshes << " meshes, " << n_nodes << " nodes and "
              << _dim * n_nodes << " coordinates" << std::endl;
#endif
    for(int type=0; type<3; type++) {
      _max_nd[type] = new int[_n_subdom * _NoLevels];
      _min_nd[type] = new int[_n_subdom * _NoLevels];

      status = H5Dread(H5Dopen(file_id, ("/NODES/MAP/MIN"+label[type]+"_SL").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_min_nd[type]);
      status = H5Dread(H5Dopen(file_id, ("/NODES/MAP/MAX"+label[type]+"_SL").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_max_nd[type]);
    }

    // Reading  _NoElements --------------
//     _off_el = new int*[_NoFamFEM];                    // volume and boundary
//     _off_el[0] = new int[_n_subdom * _NoLevels + 1];  // offset for VOLUME
//     status = H5Dread(H5Dopen(file_id, "/ELEMS/FEM0/OFF_EL", H5P_DEFAULT),
//                      H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_el[0]);
//     _off_el[1] = new int[_n_subdom * _NoLevels + 1];  // offset for boundary
//     status = H5Dread(H5Dopen(file_id, "/ELEMS/FEM1/OFF_EL", H5P_DEFAULT),
//                      H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_el[1]);
// 
//     // Number of elements (for each level)
//     _NoElements = new int*[_NoFamFEM];  // quadratic and linear
// 
//     for(int ifem = 0; ifem < _NoFamFEM; ifem++) {
//       _NoElements[ifem] = new int[_NoLevels];
//       std::ostringstream Name;    Name << "/ELEMS/FEM" << ifem << "/NExLEV";
//       status = H5Dread(H5Dopen(file_id, Name.str().c_str(), H5P_DEFAULT),
//                        H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _NoElements[ifem]);
//     }

    _off_el = new int*[_NoFamFEM];      // _NoFamFEM (volume)
    _NoElements = new int*[_NoFamFEM];  //  Number of elements (for each level)

    for(int ifem = 0; ifem < _NoFamFEM; ifem++) {
      _off_el[ifem] = new int[_n_subdom * _NoLevels + 1];  // offset for VOLUME
      std::string family=std::to_string(_GeomEl.n_q[0]);
      status = H5Dread(H5Dopen(file_id,("/ELEMS"+family+"/OFFSET_SL").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_el[ifem]);
      _NoElements[ifem] = new int[_NoLevels];
      status = H5Dread(H5Dopen(file_id,("/ELEMS"+family+"/SIZE_L").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_NoElements[ifem]);
    }

    // node mapping ------------------------------------------------------------
//     _node_map = new int*[_NoLevels + 1];        // _NoLevels+1 inverse node maps
//     int n_nodes_top = _NoNodes[_NoLevels - 1];  //  same length in file (n_nodes_top)
//     int* temp = new int[n_nodes_top];           // temp to read
// 
//     for(int ilev = 0; ilev < _NoLevels; ilev++) {   // maps from o to 2*_NoLevels-1
//       _node_map[ilev] = new int[_NoNodes[ilev]];    // inverse node map
//       // hdf5 reading
//       std::ostringstream Name;   
//       Name << "/NODES/MAP/MAP" << ilev;  // filename
//       status = H5Dread(H5Dopen(file_id, Name.str().c_str(),H5P_DEFAULT),H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
//       // inverse map
//       for(int inode = 0; inode < n_nodes_top; inode++)
//         if(temp[inode] != -1) _node_map[ilev][temp[inode]] = inode;
//     }
// 
//     // linear coarse map (storage in [2*_NoLevels])
//     _node_map[_NoLevels] = new int[_NoNodes[2 * _NoLevels]];
//     // hdf5 reading
//     std::ostringstream Name;  Name << "/NODES/MAP/MAP" << _NoLevels;  // filename
//     status = H5Dread(H5Dopen(file_id, Name.str().c_str(),H5P_DEFAULT),
//                      H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);  // read from file
// 
//     // inverse map
//     for(int inode = 0; inode < n_nodes_top; inode++)
//       if(temp[inode] != -1) _node_map[_NoLevels][temp[inode]] = inode;
    _node_map = new int*[3*_NoLevels];        // 3 node maps; _NoLevels+1 inverse node maps
    
    int n_nodes_top = _NoNodes[_NoLevels - 1];  //  same length in file (n_nodes_top)
    int* temp = new int[n_nodes_top];
    // temp to read
    for(int type=0; type<3; type++) {
      for(int ilev = 0; ilev < _NoLevels; ilev++) {   // maps from o to 2*_NoLevels-1
//         int type_nonodes = (type == 2) ? 3 : type;
//         if(ns_type==3) type_nonodes = (type == 1) ? 2 : type;
        _node_map[ilev + _NoLevels*type] = new int[_NoNodes[ilev + _NoLevels*type_nonodes[type]]];    // inverse node map
        for(int inode = 0; inode < _NoNodes[ilev + _NoLevels*type_nonodes[type]]; inode++)  _node_map[ilev + _NoLevels*type][inode] = -1;
        // hdf5 reading
        status = H5Dread(H5Dopen(file_id,("/NODES/MAP/MAP" + label[type] + "_L"+std::to_string(ilev)).c_str(), H5P_DEFAULT),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
        // inverse map
        for(int inode = 0; inode < n_nodes_top; inode++)  if(temp[inode] != -1) _node_map[ilev + _NoLevels*type][temp[inode]] = inode;
      }
    }
      delete[] temp;
//     // connectivity map -----------------------------------------------------
//     _el_map = new int*[_NoFamFEM];
//     _el_map[0] = new int[_off_el[0][_NoLevels * _n_subdom] * el_nnodes];  // Volume
//     status = H5Dread(H5Dopen(file_id, "/ELEMS/FEM0/MSH",H5P_DEFAULT),
//                      H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_map[0]); // Read from file
// 
//     _el_map[1] = new int[_off_el[1][_NoLevels * _n_subdom] * el_nnodes_bd];  // Boundary
//     status = H5Dread(H5Dopen(file_id, "/ELEMS/FEM1/MSH", H5P_DEFAULT),
//                      H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_map[1]);// Read from file
// 
//     // neighbour element map -----------------------------------------------------
//     //     const int n_sides=_GeomEl._n_sides[0];
//     _el_neighbor = new int*[_NoFamFEM - 1];  // now only for the volume family
//     _el_neighbor[0] = new int[_off_el[0][_NoLevels * _n_subdom] * el_nnodes];  // Volume
//     status = H5Dread(H5Dopen(file_id, "/ELEMS/FEM0/EL_NEIG",H5P_DEFAULT),
//                      H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_neighbor[0]); // Read from file
      
         // connectivity map -----------------------------------------------------
    _el_map = new int*[_NoFamFEM];
    _el_neighbor = new int*[_NoFamFEM];// double tmppp[1125];
    for(int ifem = 0; ifem < _NoFamFEM; ifem++) {
      _el_map[ifem] = new int[_off_el[ifem][_NoLevels * _n_subdom] * el_nnodes];  // Volume

      std::string family=std::to_string(_GeomEl.n_q[0]);
      status = H5Dread(H5Dopen(file_id, ("/ELEMS"+family+"/CONN/CONN_ALL").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,_el_map[ifem]);

      // neighbour element map -----------------------------------------------------
      const int n_sides=_GeomEl._n_sides[0];
      _el_neighbor[ifem] = new int[_off_el[ifem][_NoLevels * _n_subdom] * n_sides];  // Volume
      status = H5Dread(H5Dopen(file_id, ("/ELEMS"+family+"/EL_NEIG").c_str(), H5P_DEFAULT),
                       H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_neighbor[ifem]);
    } 
      

    // Close the files. ---------------------
    H5Fclose(file_id);    H5Fclose(file_id_xyz);    
    if(status < 0) { std::cout << " File Mesh input in data_in is missing"; abort(); }
  }  
  //=================== end restart != 0 ===================

  // READ AND STORE BOUNDARY GROUP IDS
  const int offset = _NoNodes[_NoLevels - 1];
  std::ostringstream file_bc;
  file_bc << _mgutils._inout_dir << _mgutils.get_file("INMESH");  //"/mesh.h5";
  P_INFO(std::cout << " Reading bc_id from= " << file_bc.str() << std::endl;)
  hid_t file_id = H5Fopen(file_bc.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dims[2];

  // face id vector ---------------------------------------------------------
  _NodeBDgroup = new int[offset];
  for(int i = 0; i < offset; i++) { _NodeBDgroup[i] = 0; }

  // Getting dataset
  hid_t dtset = H5Dopen(file_id, "NODES/BC", H5P_DEFAULT);
  hid_t filespace = H5Dget_space(dtset); /* Get filespace handle first. */
  hid_t status = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status < 0) {    std::cerr << "GenIc::read dims not found";  }
  else {  // reading
    assert((int)dims[0] == offset);
    status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_NodeBDgroup[0]);
  }
  H5Dclose(dtset);  H5Sclose(filespace);

  int n_nodes = _NoNodes[_NoLevels - 1];
  _normalb = new double[n_nodes * _dim];
  set_node_normal(); 

  // clean --------------------------------------------------------------------------------------
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::read_c  \n  ");)
  return;
}

/// Write mesh to hdf5 file (namefile)
///              as MGMesh class (MGMesh.h):
/// data parameters : DFL,_type_FEM
/// Vertices        : _NoNodes,_xyz
/// Connectivity    : _NoElements,_conn_map

void MGMesh::write_c(const int t_step) {
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::write_c  \n  ");)
// ===================================

  // this is my version commented. His version commented is there.
  // why is that commented?
  std::string inout_dir = _mgutils._inout_dir;
  std::string basemesh = _mgutils.get_file("BASEMESH");
  const int ndigits = stoi(_mgutils._dict_config["ndigits"]);

  // Open an existing file. ---------------
  std::ostringstream namefile;
  namefile << inout_dir << basemesh << "." << std::setw(ndigits) << std::setfill('0') << t_step << ".h5";

  // Open file to write multilevel mesh
  std::cout << " Mesh written to=     " << namefile.str() << std::endl;
  hid_t file = H5Fcreate(namefile.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  int n_meshes = _NoFamFEM * _NoLevels;

  // Writing DFL ---------------------
  int topdata[3];  topdata[0] = _dim;  topdata[1] = _NoFamFEM;  topdata[2] = _NoLevels;
  hsize_t dimsf[2];  dimsf[0] = 3;  dimsf[1] = 1;
  hid_t dtsp = H5Screate_simple(2, dimsf, NULL);
  hid_t dtset = H5Dcreate(file, "DFL", H5T_NATIVE_INT, dtsp, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, topdata);
  H5Sclose(dtsp);  H5Dclose(dtset);

  // Writing _type_FEM --------------
  dimsf[0] = n_meshes;
  dimsf[1] = 1;
  dtsp = H5Screate_simple(2, dimsf, NULL);
  dtset = H5Dcreate(file, "FEM", H5T_NATIVE_INT, dtsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _type_FEM);
  H5Sclose(dtsp);
  H5Dclose(dtset);

  // Create a group named "/MyGroup" in the file.
  hid_t group_id = H5Gcreate(file, "NODES", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t group_id2 = H5Gcreate(file, "NODES/COORD", H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);

  // writing _NoNodes ---------------------
  dimsf[0] = n_meshes;  dimsf[1] = 1;
  dtsp = H5Screate_simple(2, dimsf, NULL);
  dtset = H5Dcreate(file, "/NNODES", H5T_NATIVE_INT, dtsp, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _NoNodes);
  H5Sclose(dtsp);  H5Dclose(dtset);

  //  Writing  _xyz -------------------
  int n_nodes = _NoNodes[_NoLevels - 1];
  dimsf[0] = n_nodes;  dimsf[1] = 1;
  double* coord;  coord = new double[n_nodes];

  for(int kc = 0; kc < _dim; kc++) {
    std::ostringstream Name;    Name << "NODES/COORD/X" << kc + 1;

    dtsp = H5Screate_simple(2, dimsf, NULL);
    for(int v = 0; v < n_nodes; v++) coord[v] = _xyz[v + kc * n_nodes] /*+ _dxdydz[v+kc*n_nodes]*/;
    dtset = H5Dcreate(file, Name.str().c_str(), H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT,     H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord);
    H5Sclose(dtsp);    H5Dclose(dtset);
  }

if(_dim == 1){
  std::ostringstream Name;  Name << "NODES/COORD/X" << 2;

  dtsp = H5Screate_simple(2, dimsf, NULL);
  for(int v = 0; v < n_nodes; v++) coord[v] = 0.;
  dtset = H5Dcreate(file, Name.str().c_str(), H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord);
  H5Sclose(dtsp);  H5Dclose(dtset);

  std::ostringstream Nameb;  Nameb << "NODES/COORD/X" << 3;
  dtsp = H5Screate_simple(2, dimsf, NULL);
  for(int v = 0; v < n_nodes; v++) coord[v] = 0.;
  dtset = H5Dcreate(file, Nameb.str().c_str(), H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord);
  H5Sclose(dtsp);  H5Dclose(dtset);
}

if(_dim == 2){
  std::ostringstream Name;  Name << "NODES/COORD/X" << 3;
  dtsp = H5Screate_simple(2, dimsf, NULL);
  for(int v = 0; v < n_nodes; v++) coord[v] = 0.;
  dtset = H5Dcreate(file, Name.str().c_str(), H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord);
  H5Sclose(dtsp);  H5Dclose(dtset);
}

  // Close the group
  H5Gclose(group_id2);  H5Gclose(group_id);  delete[] coord;

  // // //
  // // //
  // // //   // writing _NoElements --------------
  // // //   dimsf[0] = n_meshes;  dimsf[1] = 1;
  // // //   dtsp =H5Screate_simple (2, dimsf, NULL);
  // // //   dtset=H5Dcreate(file,"/CONN/NELEMENTS",H5T_NATIVE_INT,dtsp,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  // // //   H5Dwrite (dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT, _NoElements);
  // // //   H5Sclose (dtsp);  H5Dclose (dtset);
  // // //
  // // //   // Writing  Connectivity ---------------
  // // //   for (int  ifem=0; ifem<_NoFamFEM;ifem++){
  // // //     for (int  ilev=0; ilev<_NoLevels;ilev++) {
  // // //       // int   *dconn;dconn=new int [_NoElements[ilev+ifem*_NoLevels]*8];
  // // //       // int  n_len=sub_conn(dconn,ifem,ilev,8);
  // // //       //std::cout <<  n_len << " n_len " << _NoElements[ilev+ifem*_NoLevels]*8<< std::endl;
  // // //       std::ostringstream Name; Name << "/CONN/FEM"<< ifem+1 << "/MSH"<< ifem*_NoLevels+ilev;
  // // //       dimsf[0] = _NoElements[ifem][ilev]*_type_FEM[ifem];
  // // //       dimsf[1] = 1;
  // // //       //dimsf[0] = n_len;  dimsf[1] = 1;
  // // //       dtsp = H5Screate_simple (2, dimsf, NULL);
  // // //       dtset = H5Dcreate
  // (file,Name.str().c_str(),H5T_NATIVE_INT,dtsp,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  // // //       H5Dwrite(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,_el_map[ifem]);
  // // //       // H5Dwrite(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,dconn);
  // // //       H5Sclose (dtsp);  H5Dclose (dtset);
  // // //     }
  // // //   }

  // Close the file -------------
  H5Fclose(file);
  // _normalb= new double [n_nodes*DIMENSION];
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::write_c  \n  ");)
  // ===================================

  return;
}

// ================================================================================================
/// This function prints the connectivity in hdf5 format
/// The changes are only for visualization of quadratic FEM
void MGMesh::print_conn_quad_hf5(
const int Level  // Level <-
) {  // ===========================================================================================
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::print_conn_quad_hf5  \n  ");)
// ================================================================================================
  std::string basemesh = _mgutils.get_file("BASEMESH");
  std::string input_dir = _mgutils._inout_dir;
  std::string connlin = _mgutils.get_file("CONNLIN");

//   int conn[8][8];
  int* gl_conn;
  // storage in hf5 (Xdmf)
  std::ostringstream namefile;  namefile << input_dir << basemesh << connlin << ".h5";
  P_INFO(std::cout << "Print mesh connectivity to =" << namefile.str() << std::endl;)
  hid_t file = H5Fcreate(namefile.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  for(int ik = 0; ik < _NoFamFEM; ik++) {
    int icount = 0;
    int indx_mesh = Level + _NoLevels * ik;
    int mode = _type_FEM[indx_mesh];
    int n_elements = _NoElements[ik][Level];  //  int  n_elements = _NoElements[indx_mesh];
//     int nsubel, nnodes;
    gl_conn = new int[n_elements * mode];
//     switch (mode) {
//         // #if DIMENSION==2
//         // -----------------------------------
//       case 9:  // Quad 9  0-4-1-5-2-6-3-7-8
//         nsubel = 4;
//         nnodes = 4;  // Quad9 -> 4 Quad4
//         gl_conn = new int[n_elements * 4 * 4];
//         conn[0][0] = 0;
//         conn[0][1] = 4;
//         conn[0][2] = 8;
//         conn[0][3] = 7;  // quad4  0-4-8-7
//         conn[1][0] = 4;
//         conn[1][1] = 1;
//         conn[1][2] = 5;
//         conn[1][3] = 8;  // quad4  4-1-5-8
//         conn[2][0] = 8;
//         conn[2][1] = 5;
//         conn[2][2] = 2;
//         conn[2][3] = 6;  // quad4  8-5-2-6
//         conn[3][0] = 7;
//         conn[3][1] = 8;
//         conn[3][2] = 6;
//         conn[3][3] = 3;  // quad4  7-8-6-3
//         break;
//
//         //=====================
//       case 6:  //  Tri6  0-3-1-4-2-5
//         nsubel = 4;
//         nnodes = 3;  //  Tri6 -> 4 Tri3
//         gl_conn = new int[n_elements * 4 * 3];
//         conn[0][0] = 0;
//         conn[0][1] = 3;
//         conn[0][2] = 5;  // quad4  0-4-8-7
//         conn[1][0] = 3;
//         conn[1][1] = 4;
//         conn[1][2] = 5;  // quad4  4-1-5-8
//         conn[2][0] = 3;
//         conn[2][1] = 1;
//         conn[2][2] = 4;  // quad4  8-5-2-6
//         conn[3][0] = 4;
//         conn[3][1] = 2;
//         conn[3][2] = 5;  // quad4  7-8-6-3
//         break;
//
//         //===================
//       case 3:  // Edge 3  0-2-1
//         nsubel = 2;
//         nnodes = 2;  // Edge 3 -> 2 Edge 2
//         gl_conn = new int[n_elements * 2 * 2];
//         conn[0][0] = 0;
//         conn[0][1] = 2;  // element 0-2
//         conn[1][0] = 2;
//         conn[1][1] = 1;  // element 1-2
//         break;
//
//         // #else
//         // ----------------------
//       case 27:  //  Hex 27 (8 Hex8)
//         nsubel = 8;
//         nnodes = 8;  //  Hex 27 -> 8 Hex8
//         gl_conn = new int[n_elements * 8 * 8];
//         conn[0][0] = 0;  conn[0][1] = 8; conn[0][2] = 20; conn[0][3] = 11;
//         conn[0][4] = 12; conn[0][5] = 21;conn[0][6] = 26; conn[0][7] = 24;
//         conn[1][0] = 8;
//         conn[1][1] = 1;
//         conn[1][2] = 9;
//         conn[1][3] = 20;
//         conn[1][4] = 21;
//         conn[1][5] = 13;
//         conn[1][6] = 22;
//         conn[1][7] = 26;
//         conn[2][0] = 11;
//         conn[2][1] = 20;
//         conn[2][2] = 10;
//         conn[2][3] = 3;
//         conn[2][4] = 24;
//         conn[2][5] = 26;
//         conn[2][6] = 23;
//         conn[2][7] = 15;
//         conn[3][0] = 20;
//         conn[3][1] = 9;
//         conn[3][2] = 2;
//         conn[3][3] = 10;
//         conn[3][4] = 26;
//         conn[3][5] = 22;
//         conn[3][6] = 14;
//         conn[3][7] = 23;
//         conn[4][0] = 12;
//         conn[4][1] = 21;
//         conn[4][2] = 26;
//         conn[4][3] = 24;
//         conn[4][4] = 4;
//         conn[4][5] = 16;
//         conn[4][6] = 25;
//         conn[4][7] = 19;
//         conn[5][0] = 21;
//         conn[5][1] = 13;
//         conn[5][2] = 22;
//         conn[5][3] = 26;
//         conn[5][4] = 16;
//         conn[5][5] = 5;
//         conn[5][6] = 17;
//         conn[5][7] = 25;
//         conn[6][0] = 24;
//         conn[6][1] = 26;
//         conn[6][2] = 23;
//         conn[6][3] = 15;
//         conn[6][4] = 19;
//         conn[6][5] = 25;
//         conn[6][6] = 18;
//         conn[6][7] = 7;
//         conn[7][0] = 26;
//         conn[7][1] = 22;
//         conn[7][2] = 14;
//         conn[7][3] = 23;
//         conn[7][4] = 25;
//         conn[7][5] = 17;
//         conn[7][6] = 6;
//         conn[7][7] = 18;
//         break;
//
//       // ---------------------------------------
//       case 10:  // Tet10 -> 8 Tet4
//         nsubel = 8;
//         nnodes = 4;  // Tet10 -> 8 Tet4
//         gl_conn = new int[n_elements * 8 * 4];
//         conn[0][0] = 0;
//         conn[0][1] = 4;
//         conn[0][2] = 6;
//         conn[0][3] = 7;
//         conn[1][0] = 4;
//         conn[1][1] = 1;
//         conn[1][2] = 5;
//         conn[1][3] = 8;
//         conn[2][0] = 5;
//         conn[2][1] = 2;
//         conn[2][2] = 6;
//         conn[2][3] = 9;
//         conn[3][0] = 7;
//         conn[3][1] = 8;
//         conn[3][2] = 9;
//         conn[3][3] = 3;
//         conn[4][0] = 4;
//         conn[4][1] = 8;
//         conn[4][2] = 6;
//         conn[4][3] = 7;
//         conn[5][0] = 4;
//         conn[5][1] = 5;
//         conn[5][2] = 6;
//         conn[5][3] = 8;
//         conn[6][0] = 5;
//         conn[6][1] = 9;
//         conn[6][2] = 6;
//         conn[6][3] = 8;
//         conn[7][0] = 7;
//         conn[7][1] = 6;
//         conn[7][2] = 9;
//         conn[7][3] = 8;
//         break;
//
//       // ---------------------------------------
//       // -----------------------------------------
//       default:  // interior 3D
//         nsubel = 1;
//         nnodes = mode;
//         gl_conn = new int[n_elements * nsubel * nnodes];
//
//         for (int n = 0; n < mode; n++) conn[0][n] = n;
//
//         break;
//     }

    // mapping
    for(int iproc = 0; iproc < _n_subdom; iproc++) {
      for(int el=_off_el[ik][iproc*_NoLevels+Level]; el<_off_el[ik][iproc*_NoLevels+Level+1]; el++) {
//         for (int se = 0; se < nsubel; se++) {
        for(int i = 0; i < mode; i++) {
          gl_conn[icount] = _el_map[ik][el * mode +_ord_FEM[i]/* conn[se][i]*/];
          icount++;
        }
      }
    }
//     }

    // Print mesh in hdf files
    std::ostringstream Name;    Name << "MSH" << ik << "CONN";
    hsize_t dimsf[2];    dimsf[0] = n_elements *mode;    dimsf[1] = 1;
    hid_t dtsp = H5Screate_simple(2, dimsf, NULL);    hid_t dtset = H5Dcreate(
        file, Name.str().c_str(),H5T_NATIVE_INT,dtsp,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gl_conn);
    H5Sclose(dtsp);    H5Dclose(dtset);

    delete[] gl_conn;
  }

  H5Fclose(file);
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print_conn_quad_hf5  \n  ");)
  return;
}


// ================================================================================================
/// This function prints the connectivity in hdf5 format
/// The changes are only for visualization of quadratic FEM
void MGMesh::print_conn_lin_hf5(const int Level  // Level <-
                               ) {     // ========================================================================================
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::print_conn_lin_hf5  \n  ");)
// ================================================================================================
  std::string basemesh = _mgutils.get_file("BASEMESH");
  std::string input_dir = _mgutils._inout_dir;
  std::string connlin = _mgutils.get_file("CONNLIN");

  int conn[8][8];
  int* gl_conn;
  // storage in hf5 (Xdmf)
  std::ostringstream namefile;
  namefile << input_dir << basemesh << connlin << ".h5";
  P_INFO(std::cout << "Print mesh connectivity to =" << namefile.str() << std::endl;)
  hid_t file = H5Fcreate(namefile.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  for(int ik = 0; ik < _NoFamFEM; ik++) {
    int icount = 0;
//     int indx_mesh = Level + _NoLevels * ik;
    int mode = _type_FEM[ik];
    int n_elements = _NoElements[ik][Level];  //  int  n_elements = _NoElements[indx_mesh];
    int nsubel, nnodes;

    switch(mode) {
    // #if DIMENSION==2
    // -----------------------------------
    case 9:  // Quad 9  0-4-1-5-2-6-3-7-8
      nsubel = 4;
      nnodes = 4;  // Quad9 -> 4 Quad4
      gl_conn = new int[n_elements * 4 * 4];
      conn[0][0] = 0;        conn[0][1] = 4;        conn[0][2] = 8;        conn[0][3] = 7;  // quad4  0-4-8-7
      conn[1][0] = 4;        conn[1][1] = 1;        conn[1][2] = 5;        conn[1][3] = 8;  // quad4  4-1-5-8
      conn[2][0] = 8;        conn[2][1] = 5;        conn[2][2] = 2;        conn[2][3] = 6;  // quad4  8-5-2-6
      conn[3][0] = 7;        conn[3][1] = 8;        conn[3][2] = 6;        conn[3][3] = 3;  // quad4  7-8-6-3
      break;

    //=====================
    case 6:  //  Tri6  0-3-1-4-2-5
      nsubel = 4;
      nnodes = 3;  //  Tri6 -> 4 Tri3
      gl_conn = new int[n_elements * 4 * 3];
      conn[0][0] = 0;        conn[0][1] = 3;        conn[0][2] = 5;  // quad4  0-4-8-7
      conn[1][0] = 3;        conn[1][1] = 4;        conn[1][2] = 5;  // quad4  4-1-5-8
      conn[2][0] = 3;        conn[2][1] = 1;        conn[2][2] = 4;  // quad4  8-5-2-6
      conn[3][0] = 4;        conn[3][1] = 2;        conn[3][2] = 5;  // quad4  7-8-6-3
      break;

    //===================
    case 3:  // Edge 3  0-2-1
      nsubel = 2;
      nnodes = 2;  // Edge 3 -> 2 Edge 2
      gl_conn = new int[n_elements * 2 * 2];
      conn[0][0] = 0;        conn[0][1] = 2;  // element 0-2
      conn[1][0] = 2;        conn[1][1] = 1;  // element 1-2
      break;

    // #else
    // ----------------------
    case 27:  //  Hex 27 (8 Hex8)
      nsubel = 8;
      nnodes = 8;  //  Hex 27 -> 8 Hex8
      gl_conn = new int[n_elements * 8 * 8];
      conn[0][0] = 0;  conn[0][1] = 8; conn[0][2] = 20; conn[0][3] = 11;
      conn[0][4] = 12; conn[0][5] = 21; conn[0][6] = 26; conn[0][7] = 24;
      conn[1][0] = 8;        conn[1][1] = 1;        conn[1][2] = 9;        conn[1][3] = 20;
      conn[1][4] = 21;        conn[1][5] = 13;        conn[1][6] = 22;        conn[1][7] = 26;
      conn[2][0] = 11;        conn[2][1] = 20;        conn[2][2] = 10;        conn[2][3] = 3;
      conn[2][4] = 24;        conn[2][5] = 26;        conn[2][6] = 23;        conn[2][7] = 15;
      conn[3][0] = 20;        conn[3][1] = 9;        conn[3][2] = 2;        conn[3][3] = 10;
      conn[3][4] = 26;        conn[3][5] = 22;        conn[3][6] = 14;        conn[3][7] = 23;
      conn[4][0] = 12;        conn[4][1] = 21;        conn[4][2] = 26;        conn[4][3] = 24;
      conn[4][4] = 4;        conn[4][5] = 16;        conn[4][6] = 25;        conn[4][7] = 19;
      conn[5][0] = 21;        conn[5][1] = 13;        conn[5][2] = 22;        conn[5][3] = 26;
      conn[5][4] = 16;        conn[5][5] = 5;        conn[5][6] = 17;        conn[5][7] = 25;
      conn[6][0] = 24;        conn[6][1] = 26;        conn[6][2] = 23;        conn[6][3] = 15;
      conn[6][4] = 19;        conn[6][5] = 25;       conn[6][6] = 18;        conn[6][7] = 7;
      conn[7][0] = 26;        conn[7][1] = 22;        conn[7][2] = 14;        conn[7][3] = 23;
      conn[7][4] = 25;        conn[7][5] = 17;        conn[7][6] = 6;        conn[7][7] = 18;
      break;

    // ---------------------------------------
    case 10:  // Tet10 -> 8 Tet4
      nsubel = 8;
      nnodes = 4;  // Tet10 -> 8 Tet4
      gl_conn = new int[n_elements * 8 * 4];
      conn[0][0] = 0;        conn[0][1] = 4;        conn[0][2] = 6;        conn[0][3] = 7;
      conn[1][0] = 4;        conn[1][1] = 1;        conn[1][2] = 5;        conn[1][3] = 8;
      conn[2][0] = 5;        conn[2][1] = 2;        conn[2][2] = 6;        conn[2][3] = 9;
      conn[3][0] = 7;        conn[3][1] = 8;        conn[3][2] = 9;        conn[3][3] = 3;
      conn[4][0] = 4;        conn[4][1] = 8;        conn[4][2] = 6;        conn[4][3] = 7;
      conn[5][0] = 4;        conn[5][1] = 5;        conn[5][2] = 6;        conn[5][3] = 8;
      conn[6][0] = 5;        conn[6][1] = 9;        conn[6][2] = 6;        conn[6][3] = 8;
      conn[7][0] = 7;        conn[7][1] = 6;        conn[7][2] = 9;        conn[7][3] = 8;
      break;

    // ---------------------------------------
    // -----------------------------------------
    default:  // interior 3D
      nsubel = 1;
      nnodes = mode;
      gl_conn = new int[n_elements * nsubel * nnodes];
      for(int n = 0; n < mode; n++) conn[0][n] = n;
      break;
    }

    // mapping
    for(int iproc = 0; iproc < _n_subdom; iproc++) {
      for(int el = _off_el[ik][iproc * _NoLevels + Level]; el < _off_el[ik][iproc * _NoLevels + Level + 1];
          el++) {
        for(int se = 0; se < nsubel; se++) {
          for(int i = 0; i < nnodes; i++) {
            gl_conn[icount] = _el_map[ik][el * mode + conn[se][i]];
            icount++;
          }
        }
      }
    }

    // Print mesh in hdf files
    std::ostringstream Name;
    Name << "MSH" << ik << "CONN";
    hsize_t dimsf[2];    dimsf[0] = n_elements * nsubel * nnodes;    dimsf[1] = 1;
    hid_t dtsp = H5Screate_simple(2, dimsf, NULL);
    hid_t dtset = H5Dcreate(file, Name.str().c_str(),H5T_NATIVE_INT,dtsp,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gl_conn);
    H5Sclose(dtsp);    H5Dclose(dtset);    delete[] gl_conn;
  }

  H5Fclose(file);
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print_conn_lin_hf5  \n  ");)
  return;
}


// ================================================================================================
/// This function prints the subdomain subdivision (proc)
void MGMesh::print_subdom_hf5(std::string filename  // filename
                             ) const { // ==========================================================
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::print_subdom_hf5  \n  ");)
// ================================================================================================
  // setup ucoord
  int n_elements = _NoElements[0][_NoLevels - 1];
  int n_subdivision = (_dim > 1) ? 4 * (_dim - 1) : 2;
  double* ucoord;  ucoord = new double[n_subdivision * n_elements];
  // storage PID -> ucoord -----------------------------
  int cel = 0;
  for(int iproc = 0; iproc < _n_subdom; iproc++) {
    for(int iel = _off_el[0][_NoLevels - 1 + iproc * _NoLevels];
        iel < _off_el[0][_NoLevels - 1 + iproc * _NoLevels + 1]; iel++) {
      for(int is = 0; is < n_subdivision; is++) ucoord[cel * n_subdivision + is] = iproc;
      cel++;
    }
  }

  // print to hdf5 format ------------------------------
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dimsf[2];  dimsf[0] = n_subdivision * n_elements;  dimsf[1] = 1;
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(
                  file_id, "PID", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ucoord);
  // clean ----------------------------------------------
  assert(status == 0);  H5Sclose(dataspace);  H5Dclose(dataset);  H5Fclose(file_id);
  delete[] ucoord;
  //   }
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print_subdom_hf5  \n  ");)
  return;
}
// ================================================================================================
/// This function prints the distance from the wall
void MGMesh::print_dist_hf5(
std::string filename,  // filename
double dist[],
std::string dir_name) const {  // =============================================================
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::print_dist_hf5  \n  ");)
// ================================================================================================
  // setup ucoord ---------------------------------------
  int n_elements = _NoElements[0][_NoLevels - 1];
  int n_subdivision = (_dim > 1) ? 4 * (_dim - 1) : 2;
  double* ucoord;  ucoord = new double[n_elements * n_subdivision];
  // storage DIST -> ucoord -----------------------------
  int cel = 0;
  for(int iproc = 0; iproc < _n_subdom; iproc++) {
    const int nel_b = _off_el[0][_NoLevels - 1 + iproc * _NoLevels];
    for(int iel = 0; iel < _off_el[0][_NoLevels - 1 + iproc * _NoLevels + 1] - nel_b; iel++) {
      for(int is = 0; is < n_subdivision; is++) ucoord[cel * n_subdivision + is] = dist[iel + nel_b];
      cel++;
    }
  }

  // print to hdf5 format -----------------------------------------
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dimsf[2];  dimsf[0] = n_elements * n_subdivision;  dimsf[1] = 1;
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(
                  file_id, dir_name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ucoord);


  // clean --------------------------------------------------------
  assert(status == 0);  H5Sclose(dataspace);  H5Dclose(dataset);  H5Fclose(file_id);
  delete[] ucoord;
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print_dist_hf5  \n  ");)
  return;
}


// ================================================================================================
void MGMesh::print_VolFrac_hf5(
std::string filename,          // filename
std::string dir_name) const {  // =============================================================
  TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::print_VolFrac_hf5  \n  ");)
// ================================================================================================
  // setup ucoord ---------------------------------------
  int n_elements = _NoElements[0][_NoLevels - 1];
  int n_subdivision = (_dim > 1) ? 4 * (_dim - 1) : 2;
  double* ucoord;  ucoord = new double[n_elements * n_subdivision];
  // storage DIST -> ucoord -----------------------------
  int cel = 0;
  for(int iproc = 0; iproc < _n_subdom; iproc++) {
    const int nel_b = _off_el[0][_NoLevels - 1 + iproc * _NoLevels];
    for(int iel = 0; iel < _off_el[0][_NoLevels - 1 + iproc * _NoLevels + 1] - nel_b; iel++) {
      for(int is = 0; is < n_subdivision; is++) ucoord[cel * n_subdivision + is] = _VolFrac[iel + nel_b];
      cel++;
    }
  }

  // print to hdf5 format -----------------------------------------
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dimsf[2];  dimsf[0] = n_elements * n_subdivision;  dimsf[1] = 1;
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(
                  file_id, dir_name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ucoord);

  // clean --------------------------------------------------------
  assert(status == 0);  H5Sclose(dataspace);  H5Dclose(dataset);  H5Fclose(file_id);
  delete[] ucoord;
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print_VolFrac_hf5  \n  ");)
  return;
}

// ================================================================================================
/// This function prints the control domain
void MGMesh::print_ctrl_dom_hf5(
std::string filename,  // filename
double ctrl_dom[],
std::string dir_name
)  // ============================================================================================
const {TRACKING_FUN(printf(" --> MGMesh.C: MGMesh::print_ctrl_dom_hf5  \n  ");)
// ================================================================================================
  // setup ucoord ---------------------------------------
  int n_elements = _NoElements[0][_NoLevels - 1];
  int n_subdivision = (_dim > 1) ? 4 * (_dim - 1) : 2;
  double* ucoord;  ucoord = new double[n_elements * n_subdivision];
  // storage ctrl_dom -> ucoord -----------------------------
  int cel = 0;
  for(int iproc = 0; iproc < _n_subdom; iproc++) {
    const int nel_b = _off_el[0][_NoLevels - 1 + iproc * _NoLevels];
    for(int iel = 0; iel < _off_el[0][_NoLevels - 1 + iproc * _NoLevels + 1] - nel_b; iel++) {
      for(int is = 0; is < n_subdivision; is++) ucoord[cel * n_subdivision + is] = ctrl_dom[iel + nel_b];
      cel++;
    }
  }

  // print to hdf5 format -----------------------------------------
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dimsf[2];  dimsf[0] = n_elements * n_subdivision;  dimsf[1] = 1;
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(
                  file_id, dir_name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ucoord);

  // clean --------------------------------------------------------
  assert(status == 0);  H5Sclose(dataspace);  H5Dclose(dataset);  H5Fclose(file_id);
  delete[] ucoord;
  TRACKING_FUN(printf(" <-- MGMesh.C: MGMesh::print_ctrl_dom_hf5  \n  ");)
  return;
}
