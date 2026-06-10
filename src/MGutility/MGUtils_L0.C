

// Ctrol femus --------------------->
#include "Printinfo_conf.h"  // petsc conf
#include "Solverlib_conf.h"  // petsc conf
// std libraries ----------------
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
// Femus Class includes --------------->
#include "MGUtils_L0.h"
// #include "MGFE.h"

// ===============================
//  MGUtils Class functions
// ===============================


// ================================================================================================
/// Constructor
// ================================================================================================
MGUtils::MGUtils(int mesh_number, std::string_view femus_dir, std::string_view app_dir, TurbUtils* turb_parameters):
  _ProbID(mesh_number),
  _femus_dir(femus_dir),
  _app_dir(app_dir),
  _DIMENSION(0)
// ================================================================================================
{TRACKING_FUN(printf(" --> MGUtils.C:  MGUtils::MGUtils (Constructor ******) \n");)
// ================================================================================================
  if (femus_dir == "") {
    char * val_femus = getenv("FEMUS_DIR");
    if(val_femus != nullptr) {_femus_dir = val_femus;} else {std::cout<<"$FEMUS_DIR ??"<<std::endl; abort();}
  }
  _contrib_dir = _femus_dir + "/contrib/"; // algebraic solver directory

  if (app_dir == "") {
    char * val_app = getenv("APP_PATH");
    if(val_app != nullptr) {_app_dir = val_app;} else {std::cout<<"$APP_PATH ??"<<std::endl; abort();}
  }

  std::string suffix = (_ProbID == -1) ? "" : std::to_string(_ProbID);
  std::string loc_par_mesh_file = _app_dir + "/DATA/param_files_msh" + suffix + ".in";
  read_param_files_msh(loc_par_mesh_file);       // from param_files_msh*.in

  // read parameters
  _inout_dir = _app_dir + "/" + get_file("INPUT_DIR");         //  RESU in/out directory
  _data_dir =  _app_dir + "/" + get_file("CONFIG_DIR");         // Data directory
  _mesh_dir=    _app_dir   + get_file("MESH_DIR");      // MESH directory
  StandardBuild();  // from all other *.in

  // interfaces
  if (_ProbID != -1) {
    std::string mesh_nameP = get_file("F_MESH_READ");  // name mesh
    int posP = mesh_nameP.find(".");  // position of "live" in str
    std::string filenameP = mesh_nameP.substr(0, posP) + "_MedToMg.med";
    _interface_mesh = filenameP;
  }

  // turbulence
  if (turb_parameters) set_turbulence_info(turb_parameters);
  TRACKING_FUN(printf(" <-- MGUtils.C:  MGUtils::MGUtils (Constructor ******) \n");)
}


// ================================================================================================
void MGUtils::StandardBuild()
// ================================================================================================
{ TRACKING_FUN(printf(" --> MGUtils.C:  MGUtils::StandardBuild \n"));
 // ===============================================================================================
 // ------------------------------------------------------------------------
//   FillSolverMap();
//   _mesh_dir = _app_dir + "/" + _dict_dir_files.find("MESH_DIR")->second;
//   if(NUM_MESH > 1) { _data_dir = _app_dir + "/DATA/DATA" + std::to_string(_ProbID) + "/"; }
  read_par();  // read parameters

  // clean and close --------------------------------------------------------------------------------
  P_INFO(print_param_files_msh();  print_par();)    // print in console parameters
  TRACKING_FUN(printf(" <-- MGUtils.C:  MGUtils::StandardBuild  \n"));

  return;
}



// ===============================================================================================
// Thisfunction READ PARAMETER FILES AND FILL RELATIVE MAPS
void MGUtils::read_par() {
// ================================================================================================
  TRACKING_FUN(printf(" --> MGUtils.C:  MGUtils::read_par \n");)
// ===============================================================================================
  // file names ------------------
  std::vector<std::string> maps;
  maps.resize(6);  // PARAMETER files
  maps[0] = get_file("BASEPARUTILS");  maps[1] = get_file("GEOM_PAR");
  maps[2] = get_file("SIM_CONFIG");    maps[3] = get_file("MAT_PROP");
  maps[4] = "Equations.in";
  maps[5] = "ControlConfig.in"; //  only for optimal  control
  double double_value;  std::string string_value;  std::string buf = "";  // read double, string, dummay

  for(int i = 1; i < 4; i++) {   // CYCLE ON FILES TO READ -------------------------------------------
    std::ostringstream filename;    filename << _data_dir << maps[i];  // file name
    std::ifstream fin;    fin.open(filename.str().c_str());  // stream file
    buf = "";
// #ifdef PRINT_INFO
//     if(fin.is_open()) { std::cout << "Init Reading = " << filename.str() << std::endl; }
// #endif
    if(fin.is_open()) {   // -------------------------------------------------------------
      while(buf != "/") {fin >> buf;}  // find "/" file start
      fin >> buf;  while(buf != "/") {
        if(buf == "#") { getline(fin, buf);}  // comment line
        else {
          if(i == 0) { fin >> double_value; set_par(buf, double_value); }
          else if(i == 1) { fin >> double_value;  set_geom_par(buf, double_value); }
          else if(i == 2 || i == 4) { fin >> string_value;  set_sim_par(buf, string_value); }
          else {fin >> double_value;  set_mat_par(buf, double_value);  }
        }
        fin >> buf;  // std::cerr <<buf.c_str() << "\n ";
      }
    }  // --------------------------------------------------------------------------------------
    else {std::cerr << "MGUtils::read_par: no parameter file found" << std::endl; abort();}
    fin.close();
    P_INFO(std::cout << "End Reading file " << filename.str() << std::endl;)
  }  // END CYCLE ON FILES TO READ ----------------------------------------------------------
  _DIMENSION=_dict_geom["DIMENSION"];
  
  // READING EQUATIONS.in FILE ---------------------------------->
  std::ostringstream filename;  filename << _data_dir << maps[4] ;  // file name
  std::ifstream fin; fin.open(filename.str().c_str());  // stream file
  buf = "";
// #ifdef PRINT_INFO
//   if(fin.is_open()) { std::cout << "Init Reading = " << filename.str() << std::endl; }
// #endif
  if(fin.is_open()) {   // -------------------------------------------------------------
    while(buf != "/") {fin >> buf;}  // find "/" file start
    fin >> buf;  while(buf != "/") {
      if(buf == "#") {getline(fin, buf);}  // comment line
      else {fin >> string_value; set_sim_par(buf, string_value); }
      fin >> buf;  // std::cerr <<buf.c_str() << "\n ";
    }
  }  // --------------------------------------------------------------------------------------
  
  int readControl =stoi(_dict_config["readControl"]); 
  if(readControl >0){ // READING ControlConfig.in FILE ---------------------------------->
  std::ostringstream filename;  filename << _data_dir << maps[5] ;  // file name
  std::ifstream fin; fin.open(filename.str().c_str());  // stream file
  buf = "";
// #ifdef PRINT_INFO
//   if(fin.is_open()) { std::cout << "Init Reading = " << filename.str() << std::endl; }
// #endif
  if(fin.is_open()) {   // -------------------------------------------------------------
    while(buf != "/") {fin >> buf;}  // find "/" file start
    fin >> buf;  while(buf != "/") {
      if(buf == "#") {getline(fin, buf);}  // comment line
      else {fin >> string_value; set_sim_par(buf, string_value);std::cerr <<buf.c_str() <<"  "<< string_value<< "\n "; }
      fin >> buf;   
    }
  } 
  } // --------------------- end ControlConfig.in-----------------------------------------------------------------
  // clear and close ------------------------------------------------------------------------------
  maps.clear();
  TRACKING_FUN(printf(" <-- MGUtils.C:  MGUtils::read_par \n");)
  return;
}


// ================================================================================================
/// This function reads file name parameter  and fill relative maps
//  file format:
//  / ----- start symbol
//  a      // data to read
//  #      comment line to skip
//  b      // data to read
// / ------ end symbol
// =================================================================================================
void MGUtils::read_temp(
const std::string& name_file_in  ///< (input) file to read
// ================================================================================================
) {TRACKING_FUN(printf(" --> MGUtils.C:  MGUtils::read_temp \n");)
// ================================================================================================
  // function setup ***************************************************************************** A
  std::ostringstream filename; filename << _app_dir << "/" << name_file_in;
  std::ifstream fin;  fin.open(filename.str().c_str());  // stream file
#ifdef PRINT_INFO
  if(fin.is_open()) { std::cout << "Init Reading = " << filename.str() << std::endl; }
#endif

  std::string buf = "";
  if(fin.is_open()) {   // -------------------------------------------------------------
    while(buf != "/") fin >> buf;  // find "/" file start
    fin >> buf;  while(buf != "/") {   //  while till to "/ " file end
      if(buf == "#") { getline(fin, buf);}  // # comment line (skip)
      else {std::string string_value; fin >> string_value; set_temp(buf, string_value); }
      fin >> buf;  // fill buf
    }
  }       // --------------------------------------------------------------------------------------
  else {std::cerr <<"MGUtils::read_par: no "<<name_file_in<<"file found"<<std::endl; abort();}

  // clean and close --------------------------------------------------------------------
  fin.close();
  P_INFO(std::cout << "End Reading file " << filename.str() << std::endl;)
  TRACKING_FUN(printf(" <-- MGUtils.C:  MGUtils::read_temp \n");)
  return;
}

// ================================================================================================
/// This function reads the file names
void MGUtils::read_param_files_msh(
  const std::string& name_file_in ///< filename where to read
)
// ================================================================================================
{ TRACKING_FUN(printf(" --> MGUtils.C:  MGUtils::read \n");)
// ================================================================================================
  // read femus dir from shell -----------------------
//   _femus_dir=getenv("FEMUS_DIR"); if(_femus_dir=="") {std::cout<<"$FEMUS_DIR ??"<<std::endl; abort();}
//   _myapp_name=getenv("FM_MYAPP"); if(_myapp_name=="") {std::cout<<"$MYAPP ??"<<std::endl; abort();}
//   std::ostringstream filename;  filename << _app_dir << "/" << name_file_in;
//   P_INFO(std::cout<<" femus_dir is="<<_femus_dir<<" MGUtils::read= "<< filename.str().c_str()<<"\n";)

//   std::ifstream fin(filename.str().c_str());  std::string buf = "";  std::string value;
     std::ifstream fin(name_file_in.c_str());  std::string buf = "";  std::string value;
  if(fin.is_open()) {
    while(!fin.eof()) {
      fin >> buf; if(buf == "#") { fin.ignore(200, '\n');}
      else { fin >> value; set_file(buf, value);}  // set new parameter
    }
  }
  else {std::cerr<<" MGFiles::read: "<<name_file_in<<" file not found"<< std::endl; abort();}
  // cleaning and check ---------------------------------------------------------------------------
  fin.close();
  //   check_dirs();
  TRACKING_FUN(printf(" <-- MGUtils.C:  MGUtils::read \n");)
  return;
}

// ============================================================
/// This file prints the file name map
void MGUtils::print_param_files_msh()
// ================================================================================================
{TRACKING_FUN(printf(" --> MGUtils.C:  MGUtils::print \n");)
 // ================================================================================================
  std::cout << "\n ================================================ ";
  std::cout << "\n Class MGFiles: " << (int)_dict_dir_files.size() << " file names: " << std::endl;
  std::map<std::string, std::string>::const_iterator pos = _dict_dir_files.begin();
  std::map<std::string, std::string>::const_iterator pos_e = _dict_dir_files.end();
  for(; pos != pos_e; pos++) {
    std::string name = pos->first;    // get the string
    std::string value = pos->second;  // get the name
    std::cout << " " << std::left << std::setw(15) << name << " = " << value << std::endl;
  }
TRACKING_FUN(printf(" <-- MGUtils.C:  MGUtils::print \n");)
  return;
}



// ================================================================================================
/// This function prints all the  parameters to stdout
void MGUtils::print_par() const
// ================================================================================================
{TRACKING_FUN(printf(" --> MGUtils.C:  MGUtils::print_par \n");)
// ================================================================================================
  std::cout << "\n ============================= ";
  std::cout << "\n   MGUtils: " << (int)_param_utils.size() << " parameters: \n";

  std::map<std::string, double>::const_iterator pos = _param_utils.begin();
  std::map<std::string, double>::const_iterator pos_e = _param_utils.end();
  for(; pos != pos_e; pos++) {
    std::string name = pos->first;  // get name
    double value = pos->second;     // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = " << std::setprecision(12) << value
              << std::endl;
  }

  std::cout << " ---------------- GEOMETRY PARAMETERS ---------------\n  ";
  std::map<std::string, double>::const_iterator pos1 = _dict_geom.begin();
  std::map<std::string, double>::const_iterator pos_e1 = _dict_geom.end();
  for(; pos1 != pos_e1; pos1++) {
    std::string name = pos1->first;  // get name
    double value = pos1->second;     // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = " << std::setprecision(12) << value
              << std::endl;
  }

  std::cout << " ---------------- SIMULATION PARAMETERS ---------------\n  ";
  std::map<std::string, std::string>::const_iterator pos2 = _dict_config.begin();
  std::map<std::string, std::string>::const_iterator pos_e2 = _dict_config.end();
  for(; pos2 != pos_e2; pos2++) {
    std::string name = pos2->first;    // get name
    std::string value = pos2->second;  // get value
    std::cout <<"  "<<std::setw(15)<<name<<" = "<< std::setprecision(12) << value<< std::endl;
  }
TRACKING_FUN(printf(" <-- MGUtils.C:  MGUtils::print_par \n");)
  return;
}

// ================================================================================================
// This function fills the problem vector
void MGUtils::FillFieldsVector(
  EquationsMap& map_str2fieldclass,
  std::vector<FIELDS>& myproblemP
// ================================================================================================
) {TRACKING_FUN(printf(" --> MGUtils.C:  MGUtils::FillFieldsVector \n");)
 // ================================================================================================
  myproblemP.clear();
  for(std::map<std::string, std::string>::iterator it = _dict_config.begin(); it != _dict_config.end(); ++it) {
    std::cout << it->first << " => " << it->second << '\n';
    if((it->first).compare(0, 3, "MG_") == 0) {
      if(stoi(it->second) != 0) {
        FIELDS ff = map_str2fieldclass._map_str2field[it->first];
        std::cout << it->first << " before " << map_str2fieldclass._map_str2field[it->first] << " "
                  << myproblemP.size() << '\n';
        myproblemP.push_back(ff);
        std::cout << it->first << " after " << map_str2fieldclass._map_str2field[it->first] << " "
                  << myproblemP.size() << '\n';
        if(it->first == "MG_DynamicalTurbulence") {
          std::cout << "Adding wall distance field \n";         myproblemP.push_back(DIST);
          std::cout << "Adding dynamical turbulence field \n";  myproblemP.push_back(MU_T);
        }
        if(it->first == "MG_ThermalTurbulence") {
          std::cout << "Adding thermal turbulence field \n";    myproblemP.push_back(ALPHA_T);
        }
      }
    }
  }

  std::cout << "\n ========= myproblemP vector ======================= \n\n";
  int sizep = myproblemP.size();
  for(int ii = 0; ii < sizep; ii++) { std::cout << "myproblemP[" << ii << "]=" << myproblemP[ii] << "\n"; }
  std::cout << "\n =============================== \n";
TRACKING_FUN(printf(" <-- MGUtils.C:  MGUtils::FillFieldsVector \n");)
  return;
}

// #if HDF5_VERSIONM == 188
// // =============================================================
// hid_t MGUtils::read_Dhdf5(hid_t file, const std::string& name, double data[]) {
//   hid_t dataset = H5Dopen(file, name.c_str());
//   hid_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
//   H5Dclose(dataset);
//   return status;
// }
//
// hid_t MGUtils::print_Dhdf5(hid_t file, const std::string& name, hsize_t dimsf[], double data[]) {
//   hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
//   hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
//   hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//   return status;
// }
//
// /// Print int data into dhdf5 file
//
// hid_t MGUtils::print_Ihdf5(hid_t file, const std::string& name, hsize_t dimsf[], int data[]) {
//   hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
//   hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
//   hid_t status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//   return status;
// }
//
// // ===========================================================================
// hid_t MGUtils::read_Ihdf5(hid_t file, const std::string& name, int data[]) {
//   hid_t dataset = H5Dopen(file, name.c_str());
//   hid_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
//   H5Dclose(dataset);
//   return status;
// }
//
// #else


// =============================================================
hid_t MGUtils::read_Dhdf5(
  hid_t file,
  const std::string& name,
  double data[]
) {// =============================================================================================
// ================================================================================================

  hid_t dataset = H5Dopen(file, name.c_str(),H5P_DEFAULT);
  hid_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  return status;
}

// ================================================================================================
hid_t MGUtils::print_Dhdf5(
  hid_t file,
  const std::string& name,
  hsize_t dimsf[],
  double data[])
{// ===============================================================================================
// ================================================================================================
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
                  H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}


// ================================================================================================
/// Print int data into dhdf5 file
hid_t MGUtils::print_Ihdf5(
  hid_t file,
  const std::string& name,
  hsize_t dimsf[],
  int data[]
) {//==============================================================================================
// ================================================================================================
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_INT, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}



// ================================================================================================
hid_t MGUtils::read_Ihdf5(
  hid_t file,
  const std::string& name,
  int data[]
) {// =============================================================================================
// ================================================================================================
  hid_t dataset = H5Dopen(file, name.c_str(),H5P_DEFAULT);
  hid_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  return status;
}
// #endif

// // ================================================================================================
// void MGUtils::FillSolverMap()
// // ================================================================================================
// {TRACKING_FUN(  printf(" --> MGUtils.C:  MGUtils::FillSolverMap \n");)
// // ================================================================================================
//
//
//
//
// //   _SolverTypeMap["CGNM"]    = CGNM;
// //   _SolverTypeMap["CGSM"]    = CGSM;
// //   _SolverTypeMap["CRM"]     = CRM;
// //   _SolverTypeMap["QMRM"]    = QMRM;
// //   _SolverTypeMap["TCQMRM"]  = TCQMRM;
// //   _SolverTypeMap["TFQMRM"]  = TFQMRM;
// //   _SolverTypeMap["BICGM"]   = BICGM;
// //   _SolverTypeMap["BICGSTABM"]   = BICGSTABM;
// //   _SolverTypeMap["MINRESM"]     = MINRESM;
// //   _SolverTypeMap["GMRESM"]      = GMRESM;
// //   _SolverTypeMap["VANKATM"]     = VANKATM;
// //   _SolverTypeMap["VANKANSM"]    = VANKANSM;
// //   _SolverTypeMap["LSQRM"]       = LSQRM;
// //   _SolverTypeMap["JACOBIM"]     = JACOBIM;
// //   _SolverTypeMap["SOR_FORWARDM"]  = SOR_FORWARDM;
// //   _SolverTypeMap["SOR_BACKWARDM"] = SOR_BACKWARDM;
// //   _SolverTypeMap["SSORM"]         = SSORM;
// //   _SolverTypeMap["RICHARDSONM"]   = RICHARDSONM;
// //   _SolverTypeMap["CHEBYSHEVM"]    = CHEBYSHEVM;
// //   _SolverTypeMap["LUMPM"]         = LUMPM;
// //   _SolverTypeMap["INVALID_SOLVERM"] = INVALID_SOLVERM;
// TRACKING_FUN( printf(" <-- MGUtils.C:  MGUtils::FillSolverMap \n");)
//   return;
// }
