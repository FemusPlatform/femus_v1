#ifndef __mgequationsmap_h__
#define __mgequationsmap_h__


// This class contains the system map
using namespace std;
#include <map>
#include <string>

class MGUtils;
class MGMesh;
class MGFEMap;

// Femus lib include
#include "MGSolver_L2_DA.h" // MGsolver at level L2 DAsolver(QLK)
#include "MGSystem_L0.h"    // basic system class level L0

#ifdef TWO_PHASE
class MGSolCC;
#endif


// ================================================================================================
class MGEquationsSystem : public MGSystem {
/// The MGEquationsSystem  class has two maps:
/// _num_equations[string=class_name="NS0"] -> # number eq in the table femus/include/EquationMaps.h
/// _equations[string=class_name="NS0"]     -> pointer to class MGSolDA*
/// and one FEM class:    
/// _mgfemap
// ================================================================================================    
 protected:
  // data --------------------------------
 int _NLevels;
 
 public:
  map<std::string, int> _num_equations;     ///< system map
  map<std::string, MGSolDA*> _equations;    ///< system map
  MGFEMap& _mgfemap;                        ///< MGFEMap class  pointer

  // Constructor / Destructor -----------------------------------
  MGEquationsSystem(
      MGUtils& mgutils_in,MGMesh& mgmesh_in, MGFEMap& mgfemap_in, 
      int np_data, int ncell_data
  );

  ~MGEquationsSystem();
  virtual void init(const std::vector<FIELDS>& pbName);  ///< Initialize
  /// Clean all substructures
  void clean();
  //---------------------------------------------------------------------------------------------
  ///@{ \name EQUATIONS GET/SET
  template <class SolverClass>
  void AddSolver(
      std::string SystemName, int nSys, int nPieceWise, int nLinear, int nQuadratic, std::string VarName);

  template <class SolverClass>
  void AddSolver(
      std::string SystemName, int nSys, int nPieceWise, int nLinear, int nQuadratic, std::string VarName,
      std::vector<FIELDS> PBname);

  inline void set_num_eqs(std::string name, int num) { _num_equations.insert(make_pair(name, num));}
  inline void set_eqs(MGSolDA* value)                { _equations.insert(make_pair(value->_eqname_D, value));}
  inline MGSolDA* get_eqs(const std::string& name)   { return _equations.find(name)->second; }
  inline const MGSolDA* get_eqs(const std::string& name) const { return _equations.find(name)->second; }
  void get_eqs_names(std::vector<std::string>& FieldsNames);
  ///@}
  //-----------------------------------------------------------------------------------------------
  
  ///@{ \name ITERATORS FOR EQUATION MAP
  typedef std::map<std::string, MGSolDA*>::iterator iterator;
  typedef std::map<std::string, MGSolDA*>::const_iterator const_iterator;

  inline iterator begin() { return _equations.begin(); }
  inline iterator end() { return _equations.end(); }
  inline const_iterator begin() const { return _equations.begin(); }
  inline const_iterator end() const { return _equations.end(); }
  ///@}
#ifdef TWO_PHASE
  void set_mgcc(MGSolCC& cc);
#endif
  void setDofBcOpIc();
  // void eqnmap_ctrl_domain(
      // const double xMin, const double xMax, const double yMin, const double yMax, const double zMin,
      // const double zMax);
  void eqnmap_ctrl_domain(const std::vector<std::vector<double>> coord);
    // void eqnmap_controlled_domain(
    //   const double xMin, const double xMax, const double yMin, const double yMax, const double zMin,
    //   const double zMax);
  void eqnmap_controlled_domain(const std::vector<std::vector<double>> coord);
   void eqnmap_lift_domain(const std::vector<std::vector<double>> coord);
  /// This function returns a value from the systems
  double GetValue(const int& ff, int flag);
  /// This function sets a value in the systems
  void SetValue(const int& ff, double value);
  void SetValueVector(const int& ff, std::vector<double> value);

  void eqnmap_steady_loop(
      const int& nmax_step,          ///< number max of steps
      const double& toll,            ///< tolerance
      const double delta_t_step_in,  //   (in)
      const int& eq_min,             ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
      const int& eq_max              ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
  );
  void set_uooold(
      const int& vec_from,           ///< source vector to be copied     (in)
      const int& vec_to,             ///< target vector                  (in)
      const double& toll,            ///< tolerance                      (in)
      const double delta_t_step_in,  //                         (in)
      const int& eq_min,             ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
      const int& eq_max              ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
  );
  int eqnmap_timestep_loop(
      const double time, const int delta_t_step_in,
      const int& eq_min,  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max   ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );
  void eqnmap_timestep_loop_control(
      const int& nmax_step,           ///< number max of steps                                (in)
      const int& it,                  ///< iteration number                                   (in)
      const double& delta_t_step_in,  ///< delta t timestep                                   (in)
      const int& eq_min,              ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
      const int& eq_max,              ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
      bool& converged                 ///< check if the solution converged (1->converged)     (out)
  );

  void eqnmap_timestep_loop_control(
      const int& nmax_step,            ///< number max of steps                                        (in)
      const int& it,                   ///< iteration number                                           (in)
      const double& delta_t_step_in,   ///< delta t timestep                                           (in)
      const int& eq_min,               ///< eq min to solve -> enum  FIELDS (equations_conf.h)         (in)
      const int& eq_max,               ///< eq max to solve -> enum  FIELDS (equations_conf.h)         (in)
      std::vector<int> controlled_eq,  ///< equations to solve and to control convergence              (in)
      bool& converged,                 ///< check if the solution converged (1->converged)             (out)
      const double& toll               ///< tolerance                                                  (in)
  );
void eqnmap_timestep_loop_and_update(
    const double time,          // real time
    const int delta_t_step_in,  // integer time
    const int& eq_min,          ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int& eq_max           ///< eq max to solve -> enum  FIELDS (equations_conf.h)
); 


  double eqnmap_timestep_nonl_solve(
      const double time,         ///< time                  (in)
      const int delta_t_step_in,   ///< time step                   (in)
      double toll,            ///< tolerance                   (in)
      int max_iter,           ///< max non linear iterations   (in)
      const int& eq_min,  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max   ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );
 void eqnmap_timestep_update(
      const double time,
      const int delta_t_step_in,
      double toll,            ///< tolerance                   (in)
      int max_iter,           ///< max non linear iterations   (in)
      const int& eq_min,  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max   ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );
 
  void eqnmap_timestep_loop_underrelaxed(
      const int& it,                      ///< tolerance                                          (in)
      const double& delta_t_step_in,      ///< delta t timestep (in)
      const int& eq_min,                  ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
      const int& eq_max,                  ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
      std::vector<double> controlled_eq,  ///< vector with the number of convergence-controlled equation  (in)
      bool& converged,                    ///< check if the solution converged (1->converged) (out)
      const double& toll                  ///< tolerance (in)
  );

  //----------------------------------------------------------------------------------------------
  ///@{ \name READ-PRINT FUNCTIONS
  void print_soln(const int t_step);         ///< Print solution
  void print_mesh_data(double vect_data[]);  ///< Print data from to mesh class (2-mesh code)
  void print_case(const int t_init);         ///< Print ic and bc

  /// Read solution
  void read_soln(const int t_step);
  ///@}
#ifdef TWO_PHASE
  void readCC(const int t_init);
#endif
  double System_functional(
      const int& ff,  ///<  eq_system
      double parameter,
      double& control  ///< parameter
  );
  double System_functional(
      const int& ff  ///<  eq_system
  );

  void movemesh();  /// Displace the mesh according to a given field
  //-------------------------------------------------------------------------------------------------
 private:
  ///@{ \name PRINT XMF/H5 FUNCTIONS
  void print_soln_xmf(const int t_step, int n_l_out, int n_c_out);
  void print_soln_h5(const int t_flag);

  void print_case_xmf(const int t_init, const int n_lines, const int n_lin);
  void print_case_h5(const int t_init);
  ///@}
  //  #ifdef TWO_PHASE
  //  /// Print xmf file CC solution
  //   void print_xmfCC(std::ofstream & out,const int t_init, const uint n_lines,const uint n_cells);
  //  /// Print CC solution in hdf5 format
  //   void print_h5CC(hid_t file,const uint flag_print,int *n_l_out,uint *n_c_out);
  // #endif
};

#include "MGSystem_L1_II.h"

#endif
