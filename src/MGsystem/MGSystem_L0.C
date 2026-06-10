// std library 
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// this class config 
#include "Printinfo_conf.h"




// conf files
// #include "MGFE_conf.h"  // Fem approximation conf file

// local class 
#include "MGSystem_L0.h"
#include "MGMesh_L0.h"
#include "MGUtils_L0.h"

// =================================================
/// This function builds the class
MGSystem::MGSystem(
    MGUtils& mgutils_in,  // file name class in
    MGMesh& mgmesh_in,    // mesh class in
    int n_pt_in,          // # system point data field
    int n_cell_in         // # system cell data  field
    )
    : _mgutils(mgutils_in),  // file name class in
      _mgmesh(mgmesh_in),    // mesh class in
      parameters() {         // parameter

    _n_data[1] = n_cell_in;  //< n cell data  (<=3)
    _n_data[0] = n_pt_in;    //< n vertex data (<=3)
                           
  return;
}

// =========================================
/// This function clears substructures
void MGSystem::clear_data() { 
  // point and cell data field
  for (int i = 0; i < _n_data[0] + _n_data[1]; i++) { delete[] _sys_data[i]; }
  return;
}

// ========================================
/// This function allocates meshdata
void MGSystem::init_data(const int mode) {        // # of data files

  const int top_lev = _mgmesh._NoLevels;
  const int n_nodes = _mgmesh._NoNodes[top_lev - 1];
  //   const int n_elements=_mgmesh._NoElements[0][top_lev-1];
  const int n_elements_tot = _mgmesh._off_el[0][_mgmesh._NoLevels * _mgmesh._n_subdom];
  //< point and cell data field
  for (int i = 0; i < _n_data[0]; i++) {
    _sys_data[i] = new double[n_nodes];
    for (int j = 0; j < n_nodes; j++) { _sys_data[i][j] = 0.; }
  }
  for (int i = _n_data[0]; i < _n_data[0] + _n_data[1]; i++) {
    _sys_data[i] = new double[n_elements_tot];
    for (int j = 0; j < n_elements_tot; j++) { _sys_data[i][j] = 0.; }
  }
  for (int i = _n_data[0] + _n_data[1]; i < 6; i++) { _sys_data[i] = NULL; }

#ifdef PRINT_INFO
  std::cout << " MGSystem::init_data with " << _n_data[0] << " point data field and " << _n_data[1]
            << " cell data field \n";
#endif
  // generating or reading data
  if (mode == 1) { read_data("input/data_sys.h5"); }
  // generating, reading and printing
  if (mode == 2) {
    read_data("input/data_sys.h5");
    print_data("../../data_in/data_sys_case.h5");
  }
  return;
}
// ========================================
/// This function reads all the  parameters
void MGSystem::read_par() {
  std::ostringstream filename;
  filename << _mgutils._data_dir << _mgutils.get_file("BASEPARAM");

  // reading from file
  std::ifstream fin;
  fin.open(filename.str().c_str());
  std::string buf = "";
  double value;

  if (fin.is_open()) {
    while (buf != "//") { fin >> buf; }
    fin >> buf;
    while (buf != "//") {
      if (buf == "#") {
        getline(fin, buf);
      }  // comment line
      else {
        fin >> value;
        set_par(buf, value);
      }  // set new parameter
      fin >> buf;
    }
  }

  else {
    std::cerr << "MGSystem::read_par: no parameters.in file found" << std::endl;
    abort();
  }
  // set non dimensional groups 
  //  set_nondimgroups();
  return;
}
// =======================================================
/// This function prints all the  parameters
void MGSystem::print_par() const {
  std::cout << "\n =================================== ";
  std::cout << "\n MGSystem::print_par(): " << (int)parameters.size() << " parameters: \n";
  std::cout << " ----------------------------------- \n";
  std::map<std::string, double>::const_iterator pos = parameters.begin();
  std::map<std::string, double>::const_iterator pos_e = parameters.end();
  for (; pos != pos_e; pos++) {
    std::string name = pos->first;
    double value = pos->second;
    std::cout << "  " << std::left << std::setw(10) << name << " = " << std::setprecision(12) << value
              << std::endl;
  }
  return;
}
// ===========================================================
/// This function sets an  external  field B
double MGSystem::F_ext(     // <-
    const double /*xx*/[],  // coordinates->
    const int /* ivarq*/    // component to return ->
    ) const {               // 

  // field Bext
  double val_return = .1;
  return val_return;
}
// =============================================
/// Print xml attrib
void MGSystem::print_xml_attrib(std::ofstream& out /*file xdmf*/, int /*nelements*/, int nodes, std::string file_name) const {  
  std::ostringstream var_name;
  const int top_lev = _mgmesh._NoLevels;
  int n_elements = _mgmesh._NoElements[0][top_lev - 1];
  for (int i = 0; i < _n_data[0] + _n_data[1]; i++) {
    int nn = 0;
    if (stoi(_mgutils._dict_config["DynamicalTurbulence"]) != 0) {
      nn = n_elements * _mgmesh._n_subdom;
    } else {
      nn = (i < _n_data[0]) ? nodes : n_elements * _mgmesh._n_subdom;
    }

    var_name.str("");
    var_name << "DATA" << i;
    //       std::string var_name = var_names[ivar];
    out << "<Attribute Name=\"" << var_name.str() << "\" AttributeType=\"Scalar\" Center=";
    if (stoi(_mgutils._dict_config["DynamicalTurbulence"]) != 0) {
      out << "\"Cell\">\n";
    } else {
      if (i < _n_data[0]) {
        out << "\"Node\">\n";
      } else {
        out << "\"Cell\">\n";
      }
    }

    out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << nn << "  " << 1
        << "\" Format=\"HDF\">  \n";
    out << file_name << ":" << var_name.str() << "\n";
    out << "</DataItem>\n"
        << "</Attribute>";
  }

  return;
}
// =============================================
/// Print xml material
void MGSystem::print_xml_mat(std::ofstream& out /*file xdmf*/, int /*nodes*/, int nelems, std::string file_name) const {  

  out << "<Attribute Name=\"MAT\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << nelems << "  " << 1
      << "\" Format=\"HDF\">  \n";
  out << file_name
      //       femus_dir << "/" << output_dir << basesol << "."
      //       << setw(ndigits) << setfill('0') << t_step << ".h5"
      << ":"
      << "COLOR"
      << "\n";
  out << "</DataItem>\n"
      << "</Attribute>";
  //   }
  return;
}
// ======================================================
/// This function read data
void MGSystem::read_data(const std::string& file_name) {
  // set the file name
  const int dim_space = _mgmesh._dim;
  const int top_lev = _mgmesh._NoLevels;
  const int n_nodes = _mgmesh._NoNodes[top_lev - 1];
  const int n_elements = _mgmesh._NoElements[0][top_lev - 1];
  std::ifstream inf(file_name.c_str());

  if (!inf) {  // autogenerate data file
               //       double xm[3]; // cell center points
    for (int ivar = _n_data[0]; ivar < _n_data[0] + _n_data[1]; ivar++) {
      for (int iel = 0; iel < n_elements; iel++) {  // element loop
        // center coordinates
        //           const int k=_mgmesh._el_map[0][iel*NDOF_FEM+NDOF_FEM-1];
        for (unsigned int idim = 0; idim <  dim_space; idim++)
        // 	    xm[idim] =_mgmesh._xyz[k+idim*n_nodes];
        // component to return ->
        {
          _sys_data[ivar][iel] = 1.;  // Fc_ext(xm,ivar);
        }

      }  // iel
    }    // ivar
    std::cout << " " << _n_data[1] << " Data cell field generated from function Fc_ext(xm,ivar) \n";

    // point data from function  Fp_ext(xm,ivar)
    for (int ivar = 0; ivar < _n_data[0]; ivar++) {
      // point loop
      for (int inode = 0; inode < n_nodes; inode++) {
        for (int idim = 0; idim <  dim_space; idim++)
        // 	   xm[idim] =_mgmesh._xyz[inode+idim*n_nodes];
        {
          _sys_data[ivar][inode] = .1;  // Fp_ext(xm,ivar);     // component to return ->
        }
      }  // inode
    }    // ivar
    std::cout << " " << _n_data[0] << " Data points fields generated from function Fp_ext(xm,ivar) \n";
  }  
  else {
    // Open an existing file.
    std::ostringstream file_name;
    file_name << _mgutils._inout_dir << "/data_sys.h5";
    std::cout << " Reading mesh from= " << file_name.str() << std::endl;
    hid_t file_id = H5Fopen(file_name.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t status = 0;
    // Reading DFL

    //     status=H5Dread(H5Dopen(file_id, "/DFLS",
    //     H5P_DEFAULT),H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,topdata);

    for (int kc = 0; kc < _n_data[0] + _n_data[1]; kc++) {
      std::ostringstream Name("");
      Name << "DATA" << kc;
      status = H5Dread(
          H5Dopen(
              file_id, Name.str().c_str()
#if HDF5_VERSIONM != 1808
                           ,
              H5P_DEFAULT
#endif
              ),
          H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _sys_data[kc]);
      assert(status == 0);
    }
  }
#ifdef PRINT_INFO
  std::cout << " Read data(MGSys): pt " << _n_data[0] << " cell " << _n_data[1] << " field(s) \n";
#endif
  return;
}
// ====================================================
/// This function prints the distance from the wall
void MGSystem::print_data_view(
    const std::string& filename,  // filename
    int i_choice,
    std::string dir_name) const {
  const int top_lev = _mgmesh._NoLevels;
   const int dim = _mgmesh._dim;
  //   const int n_nodes=_mgmesh._NoNodes[top_lev-1];
  //   const int n_elements=_mgmesh._NoElements[0][top_lev-1];
  // setup ucoord 
  int n_elements = _mgmesh._NoElements[0][top_lev - 1];
  double* ucoord;
  ucoord = new double[n_elements * 4 * (dim - 1)];
  // storage DIST -> ucoord 
  int cel = 0;
  for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++) {
    const int nel_b = _mgmesh._off_el[0][top_lev - 1 + iproc * top_lev];
    for (int iel = 0; iel < _mgmesh._off_el[0][top_lev - 1 + iproc * top_lev + 1] - nel_b; iel++) {
      for (int is = 0; is < 4 * (dim - 1); is++) {
        ucoord[cel * 4 * (dim - 1) + is] = _sys_data[i_choice][iel + nel_b];
      }
      cel++;
    }
  }
  // print to hdf5 format 
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dimsf[2];
  dimsf[0] = n_elements * 4 * (dim - 1);
  dimsf[1] = 1;

  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(
      file_id, dir_name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
#if HDF5_VERSIONM != 1808
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  hid_t status = 0;
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ucoord);
  assert(status == 0);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Fclose(file_id);

  // clean 
  delete[] ucoord;

  return;
}
// ========================================================
/// This function prints meshdata
void MGSystem::print_data(const std::string& filename) const {
  const int top_lev = _mgmesh._NoLevels;
  const int n_nodes = _mgmesh._NoNodes[top_lev - 1];
  const int n_elements = _mgmesh._NoElements[0][top_lev - 1];
  // file
  std::ostringstream name1;
  name1.str("");
  name1 << filename.c_str() << ".h5";
  hid_t fileP = H5Fcreate(name1.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //  hid_t fileP = H5Fopen (filename.c_str(),H5F_ACC_RDWR, H5P_DEFAULT );
  // Point data field

  hsize_t dimsf[2];
  dimsf[1] = 1;
  std::ostringstream name;
  for (int i = 0; i < _n_data[0] + _n_data[1]; i++) {
    name.str("");
    name << "DATA" << i;
    dimsf[0] = (i < _n_data[0]) ? n_nodes : n_elements;
    hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
    hid_t dataset = H5Dcreate(fileP, name.str().c_str(), H5T_NATIVE_DOUBLE, dataspace,
                              H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _sys_data[i]);
    assert(status == 0);  H5Sclose(dataspace);    H5Dclose(dataset);
  }

  // clean
  H5Fclose(fileP);
P_INFO(std::cout << " Print data(MGSys): pt " << _n_data[0] << " cell " << _n_data[1] << " field(s) into "<< name1.str() << " \n";)
  return;
}

// #endif
