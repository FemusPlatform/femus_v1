#pragma once


// std libraries ----------------------
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "hdf5.h"
// #ifdef TBK_EQUATIONS

// alg femus lib include -------------->
#include "Solvertype_enum.h"

// configure files -------------------->
// #include "Domain_conf.h"
// #include "MGSystem_L1.h"
// #include "MGSystem_L2.h"
#include "MGSystem_L2_EquationsMap.h"
#include "MGTurbModels_L0.h"


// #include "IbUtils.h"
// #endif
// Forwarding class -------------------
// class MGFiles;

class IbUtils;
// class TurbUtils;
// =========================================
//              MGUtils
// DATA directory
// Domain_conf.h  Equations_conf.h  MGSclass_conf.h Printinfo_conf.h Solverlib_conf.h MGFE_conf.h
// Equations.in         MaterialProperties.in    param_files_msh1.in  SimulationConfiguration.in
// GeometrySettings.in      NSproperties.in

// =========================================
/// This class contains information about parameters and directory path for
/// a problem defined by a number (_ProbID) and a name (_name)
class MGUtils {
    // ========================================================
    //                    DATA
    // ========================================================
protected:
    // Problem data ---------------------------------------
    std::map<std::string, double> _param_utils;  ///< Parameters map
    int _name;                                   ///< Problem name
    int _ProbID;                                 ///< Problem ID
   

public:
     int _DIMENSION;
    std::map<std::string, std::string> _dict_dir_files;  ///< Files map dictionary
    std::map<std::string, double>      _dict_geom;       ///< Geometry properties: bounding box generation, axisym, wall_dist
    std::map<std::string, double>      _dict_mat;        ///< Material properties: density, viscosity and other physical prop.
    std::map<std::string, std::string> _dict_config;     ///< Simulation properties: dt, print step, n. of steps, etc.
    std::map<std::string, std::string> _temp_info;       ///< Simulation properties: dt, print step, n. of steps, etc.

    std::map<std::string, SolverTypeM> _SolverTypeMap{
        {"CGNM",        CGNM},
        {"CGSM",        CGSM},
        {"CRM",         CRM},
        {"QMRM",        QMRM},
        {"TCQMRM",      TCQMRM},
        {"TFQMRM",      TFQMRM},
        {"BICGM",       BICGM},
        {"BICGSTABM",   BICGSTABM},
        {"MINRESM",    MINRESM},
        {"GMRESM",     GMRESM},
        {"VANKATM",    VANKATM},
        {"VANKANSM",   VANKANSM},
        {"LSQRM",      LSQRM},
        {"JACOBIM",    JACOBIM},
        {"SOR_FORWARDM",  SOR_FORWARDM},
        {"SOR_BACKWARDM", SOR_BACKWARDM},
        {"SSORM",         SSORM},
        {"RICHARDSONM",   RICHARDSONM},
        {"CHEBYSHEVM",    CHEBYSHEVM},
        {"LUMPM",         LUMPM},
        {"INVALID_SOLVERM",  INVALID_SOLVERM}
    };
    double _factor = 1.;
    double _press = 0.;
    // femus directories (from platform environment) ------------------
    std::string _femus_dir;    ///< Femus directory
    std::string _app_dir;      ///< Application directory
    std::string _myapp_name;   ///< Application name
    // application directories (DATA,RESU,MESH, contrib) ---------------
    std::string _data_dir;     ///< Data directory
    std::string _inout_dir;    ///< RESU directory
    std::string _mesh_dir;     ///< MESH directory
    std::string _contrib_dir;  ///< FEM directory
    std::string _interface_mesh;
    //------------------------------------------------------------------------
    // #ifdef TBK_EQUATIONS
    TurbUtils* _TurbParameters = NULL;
    IbUtils* _IBParameter = NULL;
    EquationsMap* _FieldClass;
    // #endif
    ///@}
    //   MGFiles&           _files;     ///< MGFiles class pointer

    // Constructor-Destructor --------------------
    ///@{ \name CONSTRUCTOR-DESTRUCTOR
    explicit MGUtils(
      int mesh_number = -1,
      std::string_view femus_dir = "",
      std::string_view app_dir = "",
      TurbUtils* turb_parameters = nullptr);

    void StandardBuild();
    void StandardBuild(int a);
    
    /// Destructor
    ~MGUtils() {   clean(); }
    void clean() {
        _param_utils.clear();
        _dict_dir_files.clear();
        // if(_TurbParameters != NULL) { _TurbParameters->~TurbUtils(); delete _TurbParameters; }
        //     if(_IBParameter != NULL) { _IBParameter->~IbUtils(); delete _IBParameter;  }
    }
    ///@}
    //--------------------------------------------------------------------------------------
    ///@{ \name GET/SET

    //  _dict_files dictionary
    inline std::string get_file(const std::string& name) const {
        const std::string tmp=_dict_dir_files.find(name)->second;
        if(_dict_dir_files.end()->second == tmp) {std::cout << "MGUtils: Get par failed "; abort(); }
        return tmp;
    }  ///< Return file name from the corresponding map
    inline void set_file(const std::string& name, std::string& value) {
        _dict_dir_files.insert(make_pair(name, value)); }  ///< Set file name into the corresponding map

    //  _param_utils dictionary
    inline void set_par(const std::string& name, double value) {
        _param_utils[name] = value;
    }
    inline double get_par(const std::string& name) const {
        const  double tmp=_param_utils.find(name)->second;
        if(_param_utils.end()->second == tmp) {std::cout << "MGUtils: Get par failed "; abort(); }
        return tmp;
    }  ///< Return a parameter from the corresponding map

    inline void set_name(const int value) { _name = value; }
    inline int get_name() {return _name;  }

    /// Set a parameter in the corresponding map
    inline void set_geom_par(const std::string& name, double value) { _dict_geom[name] = value; }
    inline void set_mat_par(const std::string& name, double value) {  _dict_mat[name] = value;  }

    inline void set_sim_par(const std::string& name, std::string value) {_dict_config[name] = value; }
    inline void set_temp(const std::string& name, std::string value) {_temp_info[name] = value;}


    // #ifdef TBK_EQUATIONS
    inline void set_turbulence_info(TurbUtils* Parameters) {_TurbParameters = Parameters; }
    inline void AddFieldClass(EquationsMap* FIELDclass) { _FieldClass = FIELDclass;}

    inline void set_IB_info(IbUtils* Parameters) {_IBParameter = Parameters;}
    inline void ClearTemp() {
        _temp_info.clear();
    }
    // #endif
    ///@}
    //-----------------------------------------------------------------------------------------
    ///@{ \name READ-PRINT (IN CONSTRUCTOR)
    void read_par();                                              ///< Read parameters from file
    void read_param_files_msh(const std::string& name = "/DATA/param_files.in");  ///< Read file names from file
    void read_temp(const std::string& name);                      ///< Read file names from file

    void print_param_files_msh();                                                 ///< Print in console the file names
    void print_par() const;   /// Print in console the parameters read


    // hdf5 ------------------------------------
    /// READ-PRINT HDF5
    // ----------------------------------------
    /// file   HDF5 file name (hid_t format)
    /// name   Directory name inside the HDF5 file (std::string)
    /// dimsf[]  Dimension of the vector in the directory
    /// data   Vector where put or get data

    //-----------------------------------------------------------------------------------------
    hid_t print_Dhdf5(hid_t file, const std::string& name, hsize_t dimsf[], double data[]);  ///< Print vector of double
    hid_t print_Ihdf5(hid_t file, const std::string& name, hsize_t dimsf[], int data[]);  ///< Print vector of integer
    hid_t read_Dhdf5(hid_t file, const std::string& name, double data[]);   ///< Read vector of double
    hid_t read_Ihdf5(hid_t file, const std::string& name, int data[]); ///< Read vector of integer
    ///@}

    //  ------------------------------------
    /// Field vector
    // ----------------------------------------
//     void FillSolverMap();
    void FillFieldsVector(EquationsMap& map_str2fieldclass, std::vector<FIELDS>& myproblemP);
    //  ------------------------------------
    /// BASIC ALGEBRA
    // ----------------------------------------
    void cross(const double* a, const double* b, double* res);
    double dot(const double* a, const double* b, double ab);

};

inline void MGUtils::cross(const double* a, const double* b, double* axb) {
    // a,b,res are 3D vectors
    for(int i=0; i<3; i++) axb[i]= (a[(i+1)%3]*b[(i+2)%3]-a[(i+2)%3]*b[(i+1)%3]);
    return;
}
inline double MGUtils::dot(const double* a, const double* b, double ab) {
    ab=0.;
    for(int i=0; i<3; i++) ab += a[i]*b[i]; // a,b,res are 3D vectors
    return ab;
}


