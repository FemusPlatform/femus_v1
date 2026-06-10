// std libraries --------------->
#include <math.h>
#include <iomanip>
#include <sstream>
// control Femus ---------->
#include "Printinfo_conf.h"


// Femus alg lib --------------->
#include "numeric_vectorM.h"
// femus lib class   ------------>
#include "MGSystem_L1.h"
// local inlcudes
// #include "MGFEMap.h"
#include "MGMesh_L0.h"
// #include "MGSolverDA.h"
// #include "MGUtils.h"

// conf files
// #include "Domain_conf.h"
// #include "MGFE_conf.h"


#include "MGGeom_L1_El.h"


// classes included in the map
// ====================================================
/// This function constructs all the MGSystems
MGEquationsSystem::MGEquationsSystem (
    MGUtils& mgutils_in,  // MGUtils pointer
    //   MGSystem& mgphys_in,// MGSystem pointer
    MGMesh& mgmesh_in,    // MGMesh pointer
    MGFEMap& mgfemap_in,  // MGFEMap pointer
    int np_data, int ncell_data )
    : MGSystem ( mgutils_in, mgmesh_in, np_data, ncell_data ),
      //   _mgutils(mgutils_in),
      //   _mgphys(mgphys_in),
      //   _mgmesh(mgmesh_in),
      _mgfemap ( mgfemap_in )
{
#ifdef TRACKING_FUN
    printf ( " --> MGEquationsSystem.C: MGEquationsSystem::MGEquationsSystem(Constructor *****) \n" );
#endif
    _NLevels=mgmesh_in._NoLevels;
#ifdef TRACKING_FUN
    printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::MGEquationsSystem(Constructor *****) \n" );
#endif
}

// ====================================================
MGEquationsSystem::~MGEquationsSystem() {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::MGEquationsSystem( Destructor *****) \n" );)
    clean();  // deallocates the map of equations
    _equations.clear();
    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::MGEquationsSystem( Destructor *****) \n" );)
    return;
}
// ====================================================
/// This function destroys all the MGSystems
void MGEquationsSystem::clean() {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::clean \n" );)
//     int icount=0;
//     for ( MGEquationsSystem::iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
//         icount++;
//      if(eqn->second != NULL)   delete eqn->second;
//     }
//     _equations.clear();
    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::clean \n" );)
    return;
}
// ====================================================
/// This sets dof initial and boundary conditions and sets the operators
void MGEquationsSystem::setDofBcOpIc() {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::setDofBcOpIc \n" );)
    // Reading operators
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;  // get the pointer
        mgsol->MGSol_init_B();           // init dof, GenBc, ReadOperators
    }
    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::setDofBcOpIc \n" );)
    return;
}
// ==========================================================================================
/// This function performes all the MGSystem time step routines
void MGEquationsSystem::eqnmap_steady_loop (
    const int& nmax_step,          ///< number max of steps
    const double& toll,            ///< tolerance
    const double delta_t_step_in,  ///<   (in)
    const int& eq_min,             ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int& eq_max              ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
)
{
    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_steady_loop \n" );)
    // Loop for time steps
    int NoLevels = _NLevels;
    double norm_new = 1.e-20;
    double norm_old = 1.e-20;
    double diff_norm = 1.e-20;
    double diff_norm_old = 10000.;

    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
            norm_old += mgsol->x_old[0][NoLevels - 1]->l2_norm();
        }
    }
    double time_step = delta_t_step_in;
    double time = 0;
    for ( int istep = 1; istep <= nmax_step; istep++ ) {
        // equation loop
        // -----------------------------------------------------------------------
        norm_new = 1.e-20;
        for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
            MGSolDA* mgsol = eqn->second;
            //       if(_num_equations[eqn->first] <flag_state) {
            if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
//                 NoLevels = mgsol->_NoLevels;
                mgsol->set_dt ( time_step );
                mgsol->MGTimeStep ( time, delta_t_step_in,0,2);
                norm_new += mgsol->x_old[0][NoLevels - 1]->l2_norm();
            }
        }
        diff_norm = fabs ( norm_old - norm_new );
        // ---------------------------------------------------------------------------------------
        std::cout << "\n step " << istep << ": old norm=" << norm_old << "; new norm=" << norm_new
                  << "; err  =" << diff_norm / ( norm_new ) << std::endl;
        if ( diff_norm / norm_old < toll ) { // diff_norm/norm_old < toll
            std::cout << "*** Steady state found on  n =" << istep << " ** Time step= " << time_step << " ***"
                      << std::endl;
            break;
        } else {  //   diff_norm/norm_old > toll
            std::cout << "\n  *** Steady state NOT found: relative norm difference is " << diff_norm;
            if ( diff_norm < diff_norm_old ) { // diff_norm< diff_norm_old -> step ok
                // ----------------------
                time_step *= 1.25;
                std::cout << ", increasing  step to  " << time_step << std::endl;
            }  // --------------
            else {
                time_step *= 0.95;
                std::cout << ", reduce step  to  " << time_step << std::endl;
            }
            diff_norm_old = diff_norm;
            norm_old = norm_new;
            time += time_step;
        }  // end else  diff_norm/norm_old < toll
    }
    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_steady_loop \n" );)
    return;
}
// ==========================================================================================
/// This function performes all the MGSystem time step routines
void MGEquationsSystem::set_uooold (
    const int& vec_from,           ///< source vector to be copied       (in)
    const int& vec_to,             ///< target vector                    (in)
    const double& /*toll*/,        ///< tolerance                       (in)
    const double /*delta_t_step_in*/,  ///                             (in)
    const int& eq_min,             ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int& eq_max              ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
)
{
#ifdef TRACKING_FUN
    printf ( " --> MGEquationsSystem.C: MGEquationsSystem::set_uooold \n" );
#endif
    // loop for time steps
//     int NoLevels = _NLevels;

    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
//         NoLevels = mgsol->_NoLevels;
        if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
            if ( vec_from == 3 && vec_to == 0 ) {
                mgsol->set_xooold2x();
            } else {
                mgsol->set_cp_vector (mgsol->_varst_D, vec_from, vec_to );
            }
        }
    }
#ifdef TRACKING_FUN
    printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::set_uooold \n" );
#endif
    return;
}
// ==========================================================================================
/// This function performes all the MGSystem time step routines for control
/// problem
void MGEquationsSystem::eqnmap_timestep_loop_control (
    const int& nmax_step,           ///< number max of steps (in)
    const int& it,                  ///< tolerance                                          (in)
    const double& delta_t_step_in,  ///< delta t timestep (in)
    const int& eq_min,              ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int& eq_max,              ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
    bool& converged                 ///< check if the solution converged (1->converged)     (in)
)
{
#ifdef TRACKING_FUN
    printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop_control \n" );
#endif
    // Loop for time steps
    int NoLevels = _NLevels;
    double norm_new = 1.e-20;
    double norm_old = 1.e-20;
    double diff_norm = 1.e-20;
//     double diff_norm_old = 10000.;
    double toll = 1.e-5;

    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;

        if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
//             NoLevels = mgsol->_NoLevels;
            if ( _num_equations[eqn->first] == 0 ) {
                mgsol->x_old[0][NoLevels - 1]->close();
                norm_old += mgsol->x_old[0][NoLevels - 1]->l2_norm();
            }
            //       mgsol->x_ooold[NoLevels-1]=mgsol->x_old[0][NoLevels-1];
        }
    }
    double time_step = delta_t_step_in;
    double time = time_step * it;

    for ( int istep = 1; istep <= nmax_step; istep++ ) {
        //;*(istep);
        std::cout << "\n*** Solving steady iteration n =" << istep << " ** Time= " << time <<" ***"<< std::endl;

        // equation loop
        // -----------------------------------------------------------------------
        norm_new = 1.e-20;
        for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
            MGSolDA* mgsol = eqn->second;
//             NoLevels = mgsol->_NoLevels;
            if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
                if(eqn->first=="NSOX")  {
                    std::cout << eqn->first;
                }
                mgsol->MGTimeStep ( time, delta_t_step_in,0,2);
                if ( _num_equations[eqn->first] == 0 ) {
                    norm_new += mgsol->x_old[0][NoLevels - 1]->l2_norm();
                }
            }
        }
        diff_norm = fabs ( norm_old - norm_new );
        // ---------------------------------------------------------------------------------------
        std::cout << "\n step " << istep << ": old norm=" << norm_old << "; new norm=" << norm_new
                  << "; err  =" << diff_norm / ( norm_new ) << std::endl;
        if ( diff_norm / norm_old < toll ) { // diff_norm/norm_old < toll
            std::cout <<"*** Steady state found on  n =" << istep <<" ** Time step= "<< time <<" ***"<< std::endl;
            converged = true;
            break;
        } else {  //   diff_norm/norm_old > toll
            std::cout << "\n  *** Steady state NOT found: relative norm difference is " << diff_norm / norm_old;
            if ( norm_new > 1.e+8 || std::isnan ( norm_new ) ) {
                std::cout << "\n*** Steady state NOT found ABORTING and reducing control" << endl;
                break;
            }
            //       if(diff_norm< diff_norm_old) {  // diff_norm< diff_norm_old ->
            //       step ok ----------------------
            //         time_step  *=1.25;   std::cout<<", increasing  step to  "<<
            //         time_step<<std::endl;
            //       } // --------------
            //       else {
            //         time_step  *=0.95; std::cout<<", reduce step  to  "<< time_step
            //         << std::endl;
            //       }
//             diff_norm_old = diff_norm;
            norm_old = norm_new;
            time += time_step;

        }  // end else  diff_norm/norm_old < toll
    }
#ifdef TRACKING_FUN
    printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop_control \n" );
#endif
    return;
}
//===============================================================================================
/// This function performes all the MGSystem time step routines for control
/// problem
void MGEquationsSystem::eqnmap_timestep_loop_control (
    const int& nmax_step,            ///< number max of steps (in)
    const int& it,                   ///< tolerance                                          (in)
    const double& delta_t_step_in,   ///< delta t timestep (in)
    const int& /*eq_min*/,               ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int& /*eq_max*/,               ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
    std::vector<int> controlled_eq,  ///< equations to solve and to control convergence      (in)
    bool& converged,                 ///< check if the solution converged (1->converged)     (in)
    const double& toll               ///< tolerance (in)
)
{
#ifdef TRACKING_FUN
    printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop_control \n" );
#endif
    unsigned conv_dim = controlled_eq.size();
    // Loop for time steps
    int NoLevels = _NLevels;
    double *norm_new;
    norm_new=new double [conv_dim];
    double *norm_old;
    norm_old=new double [conv_dim];
    double *diff_norm;
    diff_norm=new double[conv_dim];
    double *err_rel;
    err_rel=new double[conv_dim];
    for ( unsigned i = 0; i < conv_dim; i++ ) {
        norm_new[i] = norm_old[i] = diff_norm[i] = err_rel[i] = 1.e-20;
        for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
            MGSolDA* mgsol = eqn->second;
            if ( _num_equations[eqn->first] == controlled_eq[i] ) {
                mgsol->x_old[0][NoLevels - 1]->close();
                norm_old[i] += mgsol->x_old[0][NoLevels - 1]->l2_norm();
                // mgsol->x_nonl[NoLevels - 1]->close();
                // norm_new[i] += mgsol->x_nonl[NoLevels - 1]->l2_norm();
            }
        }
    }
    double time_step = delta_t_step_in;
    double time = time_step * it;

    for ( int istep = 1; istep <= nmax_step; istep++ ) {
        std::cout<< "\n*** Solving steady iteration n =" << istep << " ** Time= " << time <<" ***"<< std::endl;

        // equation loop
        // -----------------------------------------------------------------------
        for ( unsigned i = 0; i < conv_dim; i++ ) {
            norm_new[i] = 1.e-20;
            for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
                MGSolDA* mgsol = eqn->second;
                if ( _num_equations[eqn->first] == controlled_eq[i] ) {
                    mgsol->MGTimeStep(time,delta_t_step_in,2,2); //TODO CHECK
                    norm_new[i] += mgsol->x_old[0][NoLevels - 1]->l2_norm();
                    // norm_new[i] += mgsol->x_nonl[NoLevels - 1]->l2_norm();
                    diff_norm[i] = fabs ( norm_old[i] - norm_new[i] );
                    err_rel[i] = diff_norm[i] / norm_old[i];
                }
            }
        }
        // diff_norm=fabs(norm_old-norm_new);
        printf ( "\033[01;33m" );
        printf ( "\nConvergence control step %d:\n", istep );
        printf ( "\033[0m" );
        // ---------------------------------------------------------------------------------------
        for ( unsigned i = 0; i < conv_dim; i++ ) {
            printf (
                "Equation: %2d; Old norm: %15.8f; New norm: %15.8f; error: %15.8f\n", controlled_eq[i], norm_old[i],
                norm_new[i], err_rel[i] );
        }
        // ---------------------------------------------------------------------------------------
        bool flag = true;
        for ( unsigned i = 0; i < conv_dim; i++ ) {
            if ( err_rel[i] > toll ) {
                flag = false;
//                 diff_norm_old[i] = diff_norm[i];
                norm_old[i] = norm_new[i];
                time += time_step;
            }
            if ( norm_new[i] > 1.e+8 ) {
                std::cout << "\n  *** Steady state NOT found ABORTING and reducing control " << endl;
                break;
            }
        }
        if ( flag == true ) {
            std::cout<<"*** Steady state found on  n =" << istep << " ** Time step= " << time <<" ***"<< std::endl;
            converged = true;
            break;
        }
    }
#ifdef TRACKING_FUN
    printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop_control \n" );
#endif
    return;
}
//=================================================================================
/// This function performes all the requested MGSystem time step routines for
/// under relaxed problems. It has no loop inside, only a feedback on the
/// reached convergence of the solution norm.
void MGEquationsSystem::eqnmap_timestep_loop_underrelaxed (
    const int& it,                      ///< tolerance                                          (in)
    const double& delta_t_step_in,      ///< delta t timestep (in)
    const int& eq_min,                  ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int& eq_max,                  ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
    std::vector<double> controlled_eq,  ///< vector whose compontents are the
    ///< flags of the equations that have to
    ///< converge
    bool& converged,                    ///< check if the solution converged (1->converged) (out)
    const double& toll                  ///< tolerance (in)
)
{
#ifdef TRACKING_FUN
    printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop_underrelaxed \n" );
#endif
    unsigned conv_dim = controlled_eq.size();
    // Loop for time steps
    int NoLevels = _NLevels;
    double *norm_new;
    norm_new=new double [conv_dim];
    double *norm_old;
    norm_old=new double [conv_dim];
    double *diff_norm;
    diff_norm=new double[conv_dim];
//     double diff_norm_old[conv_dim];
    double *err_rel;
    err_rel=new double[conv_dim];
    for ( unsigned i = 0; i < conv_dim; i++ ) {
        norm_new[i] = norm_old[i] = diff_norm[i] = err_rel[i] = 1.e-20;
    }
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
//             NoLevels = mgsol->_NoLevels;
            for ( unsigned i = 0; i < conv_dim; i++ )
                if ( _num_equations[eqn->first] == controlled_eq[i] ) {
                    mgsol->x_old[0][NoLevels - 1]->close();
                    norm_old[i] += mgsol->x_old[0][NoLevels - 1]->l2_norm();
                }
        }
    }
    double time_step = delta_t_step_in;
    double time = time_step * it;

    // equation loop
    // -----------------------------------------------------------------------
    for ( unsigned i = 0; i < conv_dim; i++ ) {
        norm_new[i] = 1.e-20;
    }
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
//             NoLevels = mgsol->_NoLevels;
            mgsol->MGTimeStep ( time, delta_t_step_in,0,2 );
            // norm_new[_num_equations[eqn->first]]+= mgsol
            // ->x_old[0][NoLevels-1]->l2_norm();
            for ( unsigned i = 0; i < conv_dim; i++ )
                if ( _num_equations[eqn->first] == controlled_eq[i] ) {
                    norm_new[i] += mgsol->x_old[0][NoLevels - 1]->l2_norm();
                    diff_norm[i] = fabs ( norm_old[i] - norm_new[i] );
                    err_rel[i] = diff_norm[i] / norm_old[i];
                }
        }
    }
    // ---------------------------------------------------------------------------------------
    converged = true;
    for ( unsigned i = 0; i < conv_dim; i++ ) {
        std::cout << "\n equation " << controlled_eq[i] << ": old norm=" << norm_old[i]
                  << "; new norm=" << norm_new[i] << "; err  =" << err_rel[i] << std::endl;

        if ( diff_norm[i] > toll ) {
            converged = false;
            norm_old[i] = norm_new[i];
        }
    }
    if ( converged == true ) {
        std::cout << "*** Solution with under relaxation found" << std::endl;
    }
#ifdef TRACKING_FUN
    printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop_underrelaxed \n" );
#endif
    return;
}
// ==========================================================================================
/// This function performes all the MGSystem time step routines
/// mgsol->MGTimeStep(time, delta_t_step_in);
int MGEquationsSystem::eqnmap_timestep_loop (
    const double time,          // real time
    const int delta_t_step_in,  // integer time
    const int& eq_min,          ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int& eq_max           ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)
{
    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop \n -> mgsol->MGUpdateStep  \n" );)

    int err= 0;
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
            std::cout<< std::endl << "\033[038;5;" << _num_equations[eqn->first] + 50 << ";1m "
            << "--------------------------------------------------- \n\t"
            << " SOLUTION OF " << (eqn->first).c_str()
            << "\n ---------------------------------------------------\n\033[0m";

            mgsol->MGTimeStep ( time, delta_t_step_in,0,2 );
//           if(err1<0) return err1;
//           if(err1==1) err=1;
        }
    }
    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop \n" );)
    return err;
}
// ==========================================================================================
/// This function performes all the MGSystem time step routines
void MGEquationsSystem::eqnmap_timestep_loop_and_update (
    const double time,          // real time
    const int /*delta_t_step_in*/,  // integer time
    const int& eq_min,          ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int& eq_max           ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)
{
    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop_and_update \n -> mgsol->MGUpdateStep  \n" );)
    // solving  -> MGTimeStep_no_up() --------------------------------------------
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max )
            mgsol->MGTimeStep ( time,0,0,0);
    }
    // solution update -> MGUpdateStep() ------------------------------------------
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max )
            mgsol->MGTimeStep(time,0,0,1);
    }
    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop_and_update \n" );)
    return;
}


// ===========================================================================
double MGEquationsSystem::eqnmap_timestep_nonl_solve (
    const double time,         ///< time                  (in)
    const int delta_t_step_in,   ///< time step                   (in)
    double /*toll*/,            ///< tolerance                   (in)
    int iter,           ///< max non linear iterations   (in)
    const int& eq_min,  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int& eq_max   ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)
{
    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_solve \n -> mgsol->MGUpdateStep_no_up  \n" );)

    // loop for time steps
//     double norm_diff=2.*toll;
//     double norm_old;//[3+1];
    double norm_new=0;//[3+1];
    const int NoLevels = ( int ) ( _mgutils._dict_geom["nolevels"] );


//     for ( int kdir=0; kdir<=_mgmesh._dim ; kdir++ ) norm_old[kdir] =0.;

//     int kiter=0;
//     int err=0;
    // test tolerance  ------------------------------------------
//     while ( norm_diff>toll && kiter<max_iter ) {
//         norm_diff =0.;
//         int dir=0;
        // solution update -> MGUpdateStep() ------------------------------------------
        for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
            MGSolDA* mgsol = eqn->second;
            const int eq_pos_tab=_num_equations[eqn->first];
            int eq_field=mgsol->_field_D;
            if (eq_field >= eq_min && eq_field <= eq_max) {
                /*err=*/  mgsol->MGTimeStep( time,delta_t_step_in,iter,0);
                mgsol->x_nonl[NoLevels - 1]->close();
                norm_new/*[dir]*/ += mgsol->x_nonl[NoLevels - 1]->l2_norm();
//                 norm_diff +=fabs ( norm_new[dir]-norm_old[dir] );
                 std::cout << "\n ** non linear loop: num of equation " <<  eq_pos_tab << " time= " << time << "  eq-nsys= " << mgsol->_field_D  <<" iter=" << iter << " newnorm = " << norm_new <<  "\n";
//                 std::cout << "\n num of equation " <<  eq_pos_tab << "\n Iteration= " << kiter << " time= " << time << "  eq-nsys= " << dir
//                           << " oldnorm = " << norm_old[dir] << " newnorm = " << norm_new[dir]<< " \n non linear norm error = " << norm_diff<< "\n";
//                 norm_old[dir] = norm_new[dir];
//                 dir++;
            }
//         }

//         kiter++;
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_solve \n" );)
    return norm_new;
}
// ===============================================================
//   Update only
void MGEquationsSystem::eqnmap_timestep_update (
    const double time,           ///< time
    const int /*delta_t_step_in*/,   ///< time step
    double /*toll*/,            ///< tolerance                   (in)
    int /*max_iter*/,           ///< max non linear iterations   (in)
    const int& eq_min,  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int& eq_max   ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)
{
    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_update \n -> mgsol->MGUpdateStep  \n" );)

    // solution update -> MGUpdateStep() ---------------------------------------------------
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        int eq_field= mgsol->_field_D;
        if(eq_field>=eq_min && eq_field<= eq_max) {
            mgsol->MGTimeStep(time,0,0,1);
        }
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_timestep_loop_and_update \n" );)
    return;
}
// ==========================================================================================
/// This function sets the controlled domain
// void MGEquationsSystem::eqnmap_ctrl_domain (
//     const double xMin, const double xMax,
//     const double yMin, const double yMax,
//     const double zMin, const double zMax ) {
// 
//     TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_ctrl_domain \n" );)
// 
//     // loop for time steps
//     for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
//         MGSolDA* mgsol = eqn->second;
//         mgsol->set_ctrl_dom_B ( xMin, xMax, yMin, yMax, zMin, zMax );
//     }
// 
//     TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_ctrl_domain \n" );)
// 
//     return;
// }
void MGEquationsSystem::eqnmap_ctrl_domain (const std::vector<std::vector<double>> coord) {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_ctrl_domain \n" );)

    // loop for time steps
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        mgsol->set_ctrl_dom_B (coord);
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_ctrl_domain \n" );)

    return;
}
// ==========================================================================================
/// This function sets the controlled domain
// void MGEquationsSystem::eqnmap_controlled_domain (
//     const double xMin, const double xMax,
//     const double yMin, const double yMax,
//     const double zMin, const double zMax ) {
// 
//     TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_controlled_domain \n" );)
// 
//     // loop for time steps
//     for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
//         MGSolDA* mgsol = eqn->second;
//         mgsol->set_controlled_dom_B ( xMin, xMax, yMin, yMax, zMin, zMax );
//     }
// 
//     TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_controlled_domain \n" );)
// 
//     return;
// }
void MGEquationsSystem::eqnmap_controlled_domain (const std::vector<std::vector<double>> coord) {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_controlled_domain \n" );)

    // loop for time steps
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        mgsol->set_controlled_dom_B (coord);
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_controlled_domain \n" );)

    return;
}

void MGEquationsSystem::eqnmap_lift_domain (const std::vector<std::vector<double>> coord) {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::eqnmap_lift_domain \n" );)

    // loop for time steps
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        mgsol->set_lift_dom_B (coord);
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::eqnmap_lift_domain \n" );)

    return;
}
// ==========================================================================================
/// This function returns a value from the systems
double MGEquationsSystem::GetValue (const int& ff, int flag) {
    // loop for time steps
    double val;
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] == ff ) {
            val = mgsol->GetValue ( flag );
        }
    }
    return val;
}
// ==========================================================================================
/// This function returns a value from the systems
void MGEquationsSystem::SetValue ( const int& ff, double value ) {
    // loop for time steps
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] == ff ) {
            mgsol->SetValue ( value );
        }
    }
    return;
}
// ==========================================================================================
/// This function returns a set of values from the systems
void MGEquationsSystem::SetValueVector ( const int& ff, std::vector<double> value ) {
    // loop for time steps
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] == ff ) {
            mgsol->SetValueVector ( value );
        }
    }
    return;
}
// ==========================================================================================
/// This function prints xdmf and hdf5 file
void MGEquationsSystem::print_soln ( const int t_step) { // time step

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::print_soln \n  " );)

    const int iproc = _mgmesh._iproc;
    if ( iproc == 0 ) { // print only one processor

        print_soln_h5 ( t_step ); // print sol h5
        int n_lines = 0, n_cells = 0;
        // #ifdef TWO_PHASE  // --------- cc ----------------
        //     print_h5CC(file,t_flag,n_lines,n_lines); // print CC
        // #endif // ----------------- cc --------------------
        print_soln_xmf ( t_step, n_lines, n_cells ); // print xdmf file
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::print_soln \n" );)
    return;
}
// =================================================================
/// This function prints the attributes into the corresponding hdf5 file
void MGEquationsSystem::print_soln_h5 ( const int t_flag ) {// time flag

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::print_soln_h5 \n  " );)

    const int NoLevels = ( int ) ( _mgutils._dict_geom["nolevels"] );
    const int ndigits = stoi ( _mgutils._dict_config["ndigits"] );

    // file  ---------------------------------------------
    // file name
    std::ostringstream filename;
    filename << _mgutils._inout_dir << _mgutils.get_file ( "BASESOL" ) << "." << setw ( ndigits ) << setfill ( '0' )
             << t_flag << ".h5";
    // open file for hf5 storage
    hid_t file = H5Fcreate ( filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    H5Fclose ( file );

    // print all systems ---------------------------
    MGEquationsSystem::const_iterator pos = _equations.begin();
    MGEquationsSystem::const_iterator pos_e = _equations.end();
    for ( ; pos != pos_e; pos++ ) {
        MGSolDA* mgsol = pos->second;
        mgsol->print_u_D ( filename.str(), NoLevels - 1 );
#ifdef HAVE_MED
        if ( _mgutils._dict_config["MG_ControlTemperature"] != "" )
            if ( stoi ( _mgutils._dict_config["MG_ControlTemperature"] ) != 0 ) {
                mgsol->print_weight_med_D ( "mesh_sol.med", NoLevels - 1 );
            }
#endif
    }
    // printing cell system data (MGSystem-> printdata)
    for ( int idata = 0; idata < _n_data[0] + _n_data[1]; idata++ ) {
        std::ostringstream dir_name;
        dir_name << "DATA" << idata;
        print_data_view ( filename.str(), idata, dir_name.str() );
    }
    if ( _mgutils._dict_config["MG_ImmersedBoundary"] != "" )
        if ( stoi ( _mgutils._dict_config["MG_ImmersedBoundary"] ) != 0 ) {
            _mgmesh.print_VolFrac_hf5 ( filename.str(), "Piece_VolFrac" );
        }

    if ( _mgutils._dict_config["MG_DynamicalTurbulence"] != "" )
        if ( stoi ( _mgutils._dict_config["MG_DynamicalTurbulence"] ) != 0 ) {
            _mgmesh.print_dist_hf5 ( filename.str(), _mgmesh._yplus, "YPLUS" );
        }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::print_soln_h5 \n  " );)
    return;
}
// void MGEquationsSystem::print_soln_h5(const int t_flag) {
//
//   const int NoLevels = (int)(_mgutils._dict_geom["nolevels"]);
//   const int ndigits = stoi(_mgutils._dict_config["ndigits"]);
//
//   // file  ---------------------------------------------
//   // file name
//   std::ostringstream filename;
//   filename << _mgutils._inout_dir << _mgutils.get_file("BASESOL") << "." << setw(ndigits) << setfill('0')
//            << t_flag << ".h5";
//   // open file for hf5 storage
//   hid_t file = H5Fcreate(filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//   H5Fclose(file);
//
//   // print all systems ---------------------------
//   MGEquationsSystem::const_iterator pos = _equations.begin();
//   MGEquationsSystem::const_iterator pos_e = _equations.end();
//   for (; pos != pos_e; pos++) {
//     MGSolDA* mgsol = pos->second;
//     mgsol->print_u_D(filename.str(), NoLevels - 1);
//     if (_mgutils._dict_config["MG_ControlTemperature"] != "")
//       if (stoi(_mgutils._dict_config["MG_ControlTemperature"]) != 0)
//         mgsol->print_weight_med_D("mesh_sol.med", NoLevels - 1);
//   }
//   // printing cell system data (MGSystem-> printdata)
//   for (int idata = 0; idata < _n_data[0] + _n_data[1]; idata++) {
//     std::ostringstream dir_name;
//     dir_name << "DATA" << idata;
//     print_data_view(filename.str(), idata, dir_name.str());
//   }
//   if (_mgutils._dict_config["MG_ImmersedBoundary"] != "")
//     if (stoi(_mgutils._dict_config["MG_ImmersedBoundary"]) != 0)
//       _mgmesh.print_VolFrac_hf5(filename.str(), "Piece_VolFrac");
//
//   if (_mgutils._dict_config["MG_DynamicalTurbulence"] != "")
//     if (stoi(_mgutils._dict_config["MG_DynamicalTurbulence"]) != 0)
//       _mgmesh.print_dist_hf5(filename.str(), _mgmesh._yplus, "YPLUS");
//
//   return;
// }

// ===================================================================
/// It prints the attributes in  Xdmf format for one time step
void MGEquationsSystem::print_soln_xmf ( const int t_step, int /*n_lines*/, int /*n_cells*/ ) {
    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::print_soln_xmf \n  " );)

    const int NoLevels = ( int ) ( _mgutils._dict_geom["nolevels"] );
    const int ndigits = stoi ( _mgutils._dict_config["ndigits"] );

    //  Mesh ----------------------
    const MGMesh& mgmesh = _mgmesh;
    int n_nodes = mgmesh._NoNodes[NoLevels - 1];
    int n_elements = mgmesh._NoElements[0][NoLevels - 1];
    const int n_subdom = mgmesh._n_subdom;

    // get parameters
    std::string inout_dir = _mgutils._inout_dir;
    std::string basesol = _mgutils.get_file ( "BASESOL" );
    std::string basemesh = _mgutils.get_file ( "BASEMESH" );
    std::string contrib_dir = _mgutils.get_file ( "CONTRIB_DIR" );
    //   std::string aux_xdmf   = _mgutils.get_file("AUX_XDMF");
    std::string connlin = _mgutils.get_file ( "CONNLIN" );
    std::string basecase = _mgutils.get_file ( "BASECASE" );
    const double dt = stod ( _mgutils._dict_config["dt"] );

    // files
    //   std::ostringstream conn_file; // connectivity file (mesh_conn_lin.h5)
    //   conn_file /* <<femus_dir  <<  "/" << appl_dir << "/"  << myapp << "/" <<
    //   input_dir */ << basemesh;
    std::ostringstream topol_file;  // topology file file (mesh_conn_lin.h5)
    topol_file << basemesh << connlin << ".h5";
    //   conn_file << ".h5";
    std::ostringstream coord_time_file;  // connectivity file (mesh_conn_lin.h5)
    coord_time_file /*<<femus_dir  << "/"<< output_dir << "/"*/ << basemesh << "." << std::setw ( ndigits )
            << std::setfill ( '0' ) << t_step << ".h5";

    std::ostringstream filename;  //  solution file xmf
    filename << inout_dir << basesol << "." << setw ( ndigits ) << setfill ( '0' ) << t_step;
    std::ostringstream casefilename;  //  solution file xmf
    casefilename << basecase << "." << setw ( ndigits ) << setfill ( '0' ) << 0 << ".h5";

    std::ostringstream attr_file;  //  solution file h5
    /* attr_file << filename.str()<< ".h5"; */
    attr_file << basesol << "." << setw ( ndigits ) << setfill ( '0' ) << t_step << ".h5";
    filename << ".xmf";

    // solution file  xmf
    std::ofstream out ( filename.str().c_str() );
    std::cout << " Solution written to= " << filename.str().c_str() << std::endl;
    //  ++++++++++++ Header ++++++++++++++
    out << "<?xml version=\"1.0\" ?> \n";
    //   out << "<!DOCTYPE Xdmf SYSTEM ";
    //   out <<  "\"" << aux_xdmf << "\" > \n";
    out << "<Xdmf> \n"
        << "<Domain> \n"
        << "<Grid Name=\"Mesh\"> \n";
    // time
    const double restartime =
        ( stoi ( _mgutils._dict_config["restart"] ) != 0 ) ? stod ( _mgutils._dict_config["restartime"] ) : 0.0;
    const int t_in = stoi ( _mgutils._dict_config["restart"] );
    out << "<Time Value =\"" << restartime + ( t_step - t_in ) * dt << "\" /> \n";
    // +++++ Topology ++++++++++
    _mgmesh.print_xmf_topology ( out, topol_file.str(), NoLevels - 1, 0 );
    // +++++++  Geometry +++++++++++++++++
    _mgmesh.print_xmf_geometry ( out, coord_time_file /*conn_file*/.str(), NoLevels - 1, 0 );
    // ++++  Attributes ++++++++++++
    MGEquationsSystem::const_iterator pos1 = _equations.begin();
    MGEquationsSystem::const_iterator pos1_e = _equations.end();
    for ( ; pos1 != pos1_e; pos1++ ) {
        MGSolDA* mgsol = pos1->second;
        mgsol->print_u_xdmf_D ( out, n_nodes, n_elements * n_subdom, attr_file.str() );
    }
    //   print_xml_attrib(out,n_elements,n_nodes,attr_file.str());
    //   printining cell attributes
    if ( _n_data[0] + _n_data[1] > 0 ) {
        print_xml_attrib ( out, n_elements, n_nodes, attr_file.str() );
    }

    if ( _mgutils._dict_config["MG_FluidStructure"] != "" )
        if ( stoi ( _mgutils._dict_config["MG_FluidStructure"] ) != 0 ) {
            print_xml_mat ( out, n_nodes, n_elements * n_subdom, casefilename.str() );
        }
    // #ifdef TWO_PHASE
    //   // print of CC
    //   print_xmfCC(out,t_step,n_lines,n_cells);
    // #endif
    if ( _mgutils._dict_config["MG_ImmersedBoundary"] != "" )
        if ( stoi ( _mgutils._dict_config["MG_ImmersedBoundary"] ) != 0 ) {
            // ----------------------------------------------------------
            out << "<Attribute Name=\"Piece_VolFrac\" AttributeType=\"Scalar\" "
                "Center=\"Cell\">\n";
            out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_elements * n_subdom << "  "
                << 1 << "\" Format=\"HDF\">  \n";
            out << attr_file.str()
                //       femus_dir << "/" << output_dir << basesol << "."
                //       << setw(ndigits) << setfill('0') << t_step << ".h5"
                << ":Piece_VolFrac\n";
            out << "</DataItem>\n"
                << "</Attribute>";
        }

    if ( _mgutils._dict_config["MG_DynamicalTurbulence"] != "" )
        if ( stoi ( _mgutils._dict_config["MG_DynamicalTurbulence"] ) != 0 ) {
            out << "<Attribute Name=\"YPLUS\" AttributeType=\"Scalar\" "
                "Center=\"Cell\">\n";
            out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_elements * n_subdom << "  "
                << 1 << "\" Format=\"HDF\">  \n";
            out << attr_file.str()
                //       femus_dir << "/" << output_dir << basesol << "."
                //       << setw(ndigits) << setfill('0') << t_step << ".h5"
                << ":YPLUS\n";
            out << "</DataItem>\n"
                << "</Attribute>";
        }

    out << "</Grid>\n"
        << "</Domain> \n"
        << "</Xdmf> \n";
    out.close();

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::print_soln_xmf \n  " );)
    return;
}
// ========================================================================
/// This function read the solution for all the system (restart)
void MGEquationsSystem::read_soln ( const int t_step ) {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::read_soln \n  " );)

    // --------------------------------------------------------------
    const int ndigits = stoi ( _mgutils._dict_config["ndigits"] );
    //   const int ndigits  = _mgutils.get_par("ndigits");
    const int restart_lev_flag = stoi ( _mgutils._dict_config["restart_lev"] );

    // open file -----------------------------
    std::ostringstream namefile;
    namefile << _mgutils._inout_dir << _mgutils.get_file ( "BASESOL" ) << "." << setw ( ndigits ) << setfill ( '0' )
             << t_step << ".xmf";

    std::cout<< namefile.str().c_str() <<std::endl;

    P_INFO(std::cout << "\n MGEquationsSystem::read_soln: Reading time  from " << namefile.str().c_str();)
    std::ifstream in;
    in.open ( namefile.str().c_str() ); // associate the file stream with the name of the file
    if ( !in.is_open() ) {
        std::cout << " MGCase: restart .xmf file not found " << std::endl;
        abort();
    }

    // reading time from xmf file --------------
    std::string buf = "";
    while ( buf != "<Time" ) {
        in >> buf;
    }
    in >> buf >> buf;
    buf = buf.substr ( 2, buf.size() - 3 );
    // create an istringstream from a string
    std::istringstream buffer ( buf );
    double restart_time;
    buffer >> restart_time;

    // add parameter to system
    //   _mgutils.set_par("restartime",restart_time);
    _mgutils.set_sim_par ( "restartime", std::to_string ( restart_time ) );
    // ---------------------------------------------------
    // reading data from  sol.N.h5
    // ---------------------------------------------------
    // file name -----------------------------------------
    namefile.str ( "" ); // empty string
    namefile << _mgutils._inout_dir << _mgutils.get_file ( "BASESOL" ) << "." << setw ( ndigits ) << setfill ( '0' )
             << t_step << ".h5";

    P_INFO(std::cout << "\n MGEquationsSystem::read_soln: Reading from file " << namefile.str().c_str() << std::endl;)
    // loop reading over the variables ---------------------
    for ( MGEquationsSystem::const_iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        mgsol->read_u_D ( namefile.str(), restart_lev_flag );
    }  //  loop --------------------------------------------------------

    // #ifdef TWO_PHASE
    //   readCC(t_step);
    // //   _mgsys.get_mgcc()->read_fine ( t_init,_mgsys.get_mgcc()->_Nlev-1 );
    // //   readCC ( t_step );
    // #endif

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::read_soln \n  " );)
    return;
}

// ========================================================================
/// This function print data from single class equation to the mesh system  on
/// vect_data
void MGEquationsSystem::print_mesh_data ( double vect_data[] ) {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::print_mesh_data` \n  " );)

    // loop reading/printing over the equation ---------------------
    int count = 0;
    //   eqn=_equations.begin();
    for ( MGEquationsSystem::const_iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        if ( _mgutils.get_file ( "MESHNUMBER" ) == "mesh1" ) {
            if ( count == 1 ) {
                MGSolDA* mgsol = eqn->second;
                mgsol->print_ext_data ( vect_data ); // compute from a system/mesh to vect_data
            }
            count++;
        }
        if ( _mgutils.get_file ( "MESHNUMBER" ) == "mesh2" ) {
            if ( count == 0 ) {
                MGSolDA* mgsol = eqn->second;
                mgsol->print_ext_data ( vect_data ); // compute from a system/mesh to vect_data
            }
            count++;
        }
    }  //  loop --------------------------------------------------------

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::print_mesh_data \n  " );)
    return;
}
// ================================================================================================
/// This function prints initial and boundary data in xdmf+hdf5 format
void MGEquationsSystem::print_case ( const int t_init ) {

    TRACKING_FUN(printf ( "--> MGEquationsSystem.C: MGEquationsSystem::print_case \n  " );)

    const int iproc = _mgmesh._iproc;

    if ( iproc == 0 ) {       // print only one processor
        print_case_h5 ( t_init ); // ic+bc print format h5

        int n_lines = 0, n_cells = 0;  // for VOF
        // #ifdef TWO_PHASE
        //     print_h5CC(hid_t file,t_init,&n_lines,&n_cells);
        // #endif
        print_case_xmf ( t_init, n_lines, n_cells ); // xml format
    }

    TRACKING_FUN(printf ( "<-- MGEquationsSystem.C: MGEquationsSystem::print_case \n  " );)
    return;
}
// ================================================================================================
/// This function prints initial and boundary data in hdf5 fromat
/// in the file case.h5
void MGEquationsSystem::print_case_h5 ( const int t_init ) {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::print_case_h5 \n  " );)

    const int NoLevels = ( int ) ( _mgutils._dict_geom["nolevels"] );
    const int ndigits = stoi ( _mgutils._dict_config["ndigits"] );
    std::string output_dir = _mgutils._inout_dir;
    std::string basecase = _mgutils.get_file ( "BASECASE" );
    //  Mesh ---- ---------------------------------------------
    const MGMesh& mgmesh = _mgmesh;
    // file ---------------------------------------
    std::ostringstream filename;  // file name
    filename << output_dir << basecase << "." << setw ( ndigits ) << setfill ( '0' ) << t_init << ".h5";
    hid_t file = H5Fcreate ( filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    mgmesh.print_subdom_hf5 ( filename.str() ); // PID (n processor)
    if ( _mgutils._dict_config["MG_DynamicalTurbulence"] != "" )
        if ( stoi ( _mgutils._dict_config["MG_DynamicalTurbulence"] ) != 0 ) {
            mgmesh.print_dist_hf5 ( filename.str(), mgmesh._dist,"DIST" ); // PID (n processor)
            mgmesh.print_dist_hf5 ( filename.str(), mgmesh._yplus,"YPLUS" ); // PID (n processor)
        }
    H5Fclose ( file );
    // loop over all systems ---------------------------
    MGEquationsSystem::const_iterator pos = _equations.begin();
    MGEquationsSystem::const_iterator pos_e = _equations.end();
    for ( ; pos != pos_e; pos++ ) {
        MGSolDA* mgsol = pos->second;
        mgsol->print_u_D ( filename.str(), NoLevels - 1 ); // initial solution
        mgsol->print_bc_D ( filename.str(), NoLevels - 1 ); // boundary condition
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::print_case_h5 \n  " );)
    return;
}
// ================================================================================================
/// It prints the Xdmf file to read the initial and boundary conditions
void MGEquationsSystem::print_case_xmf ( const int t_init, const int /*n_lines*/, const int /*n_cells*/ ) {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::print_case_xmf \n  " );)

    // file ----------------------------------------------
    std::string inout_dir = _mgutils._inout_dir;
    std::string basecase = _mgutils.get_file ( "BASECASE" );
    std::string basemesh = _mgutils.get_file ( "BASEMESH" );
    std::string contrib_dir = _mgutils.get_file ( "CONTRIB_DIR" );
    std::string aux_xdmf = _mgutils.get_file ( "AUX_XDMF" );
    std::string connlin = _mgutils.get_file ( "CONNLIN" );
    const int NoLevels = ( int ) ( _mgutils._dict_geom["nolevels"] );
    const int ndigits = stoi ( _mgutils._dict_config["ndigits"] );

    // files
    std::ostringstream conn_file;  // connectivity file (mesh_conn_lin.h5)
    conn_file << basemesh;
    std::ostringstream topol_file;  // topology file file (mesh_conn_lin.h5)
    topol_file << conn_file.str() << connlin << ".h5";
    conn_file << ".h5";
    std::ostringstream filename;  //  solution file xmf
    filename << inout_dir << basecase << "." << setw ( ndigits ) << setfill ( '0' ) << t_init << ".xmf";

    //  Mesh ---- ---------------------------------------------
    const MGMesh& mgmesh = _mgmesh;
    int n_nodes = mgmesh._NoNodes[NoLevels - 1];
    int n_elements = mgmesh._NoElements[0][NoLevels - 1];
    std::string var_name[3];
    std::string var_type[3];

    //   File xdmf -------------------------------
    std::ofstream out ( filename.str().c_str() );
    out << "<?xml version=\"1.0\" ?> \n";
    out << "<!DOCTYPE Xdmf SYSTEM ";
    out << "\"" << aux_xdmf << "\" \n";
    //    out << " [ <!ENTITY HeavyData \"\"> ] ";
    out << ">\n";
    out << "<Xdmf> \n"
        << "<Domain> \n"
        << "<Grid Name=\"Mesh\"> \n";
    // +++++ Topology ++++++++++
    _mgmesh.print_xmf_topology ( out, topol_file.str(), NoLevels - 1, 0 );
    // +++++++  Geometry +++++++++++++++++
    _mgmesh.print_xmf_geometry ( out, conn_file.str(), NoLevels - 1, 0 );

    // ++++  Attributes for each system +++++++++++++++++++
    MGEquationsSystem::const_iterator pos1 = _equations.begin();
    MGEquationsSystem::const_iterator pos1_e = _equations.end();
    for ( ; pos1 != pos1_e; pos1++ ) {
        MGSolDA* mgsol = pos1->second;
//         int ivar=0;
       for ( int ivar = 0; ivar < mgsol->_varst_D[0] + mgsol->_varst_D[1] +  mgsol->_varst_D[2]; ivar++ ) {
            // Volume and boundary conditions
            var_name[0] = mgsol->_var_names_D[ivar];
            var_name[1] = "bc_"+var_name[0];
//             var_name[2] = var_name[0] + "vl";
            var_type[0] = "Float";
            var_type[1] = "Int";
            int n_ibvar=1; if(ivar==0) n_ibvar=2; // only the first variable has bc
//             var_type[2] = "Int";
//             for ( int ibvar = 0; ibvar < 3; ibvar++ ) {
                for ( int ibvar = 0; ibvar <n_ibvar; ibvar++ ) {
                out << "<Attribute Name=\"" << var_name[ibvar] << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
                out << "<DataItem  DataType=\"" << var_type[ibvar].c_str() << "\" Precision=\"8\" Dimensions=\""
                    << n_nodes << "  " << 1 << "\" Format=\"HDF\">  \n";
                out << /* femus_dir << "/" << output_dir <<*/ basecase << "." << setw ( ndigits ) << setfill ( '0' )
                    << t_init << ".h5"
                    << ":" << var_name[ibvar].c_str() << "\n";
                out << "</DataItem>\n"
                    << "</Attribute>\n";
            }
        }
    }

    // ----------------------------------------------------------
    out << "<Attribute Name=\"PID\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_elements * _mgmesh._GeomEl.n_se[0]
        << "  " << 1 << "\" Format=\"HDF\">  \n";
    out << basecase << "." << setw ( ndigits ) << setfill ( '0' ) << t_init << ".h5"
        << ":PID\n";
    out << "</DataItem>\n"
        << "</Attribute>\n";
    if ( _mgutils._dict_config["MG_DynamicalTurbulence"] != "" ) {
        if ( stoi ( _mgutils._dict_config["MG_DynamicalTurbulence"] ) != 0 ) {
            // ----------------------------------------------------------
            out << "<Attribute Name=\"DIST\" AttributeType=\"Scalar\" "
                "Center=\"Cell\">\n";
            out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
                << n_elements * _mgmesh._GeomEl.n_se[0] << "  " << 1 << "\" Format=\"HDF\">  \n";
            out << basecase << "." << setw ( ndigits ) << setfill ( '0' ) << t_init << ".h5"
                << ":DIST\n";

            out << "</DataItem>\n"
                << "</Attribute>\n";
            out << "<Attribute Name=\"YPLUS\" AttributeType=\"Scalar\" "
                "Center=\"Cell\">\n";
            out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
                << n_elements * _mgmesh._GeomEl.n_se[0] << "  " << 1 << "\" Format=\"HDF\">  \n";
            out << basecase << "." << setw ( ndigits ) << setfill ( '0' ) << t_init << ".h5"
                << ":YPLUS\n";
            out << "</DataItem>\n"
                << "</Attribute>\n";
        }
    }
    // #ifdef TWO_PHASE
    //   // print of CC
    //   print_xmfCC(out, t_init,n_lines,n_cells);
    // #endif
    out << "</Grid>\n"
        << "</Domain> \n"
        << "</Xdmf> \n";
    out.close();

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::print_case_xmf \n  " );)
    return;
}
// ============================================================================
double MGEquationsSystem::System_functional (
    const int& ff,  ///< initial time iteration
    double parameter,
    double& control  ///< step time
)
{
    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::System_functional \n  " );)

    double value = 0.;
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] == ff ) {
//             const int NoLevels = mgsol->_NoLevels;
            value += mgsol->MGFunctional_B ( parameter, control );
            std::cout << "  functional " << _num_equations[eqn->first] << "\n -----";
        }
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::System_functional \n  " );)
    return value;
}
// =============================================================================
double MGEquationsSystem::System_functional (
    const int& ff  ///< initial time iteration
)
{
    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::System_functional \n  " );)

    double value = 0.;
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;
        if ( _num_equations[eqn->first] == ff ) {
            // ************************************
            value +=0.;// mgsol->MGFunctional_B ();
        // *******************************************
            std::cout << "  functional " << _num_equations[eqn->first] << "\n -----";
        }
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::System_functional \n  " );)
    return value;
}
// =============================================================================
void MGEquationsSystem::get_eqs_names ( std::vector<string>& FieldsNames ) {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::get_eqs_names \n  " );)

    for ( MGEquationsSystem::iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {

        FieldsNames.push_back ( eqn->first );
        //     std::cout<<" MGEquationsSystem::get_eqs_names "<< eqn->first
        //     <<std::endl;
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::get_eqs_names \n  " );)
    return;
}
// =============================================================================
/// This function prints initial and boundary data in hdf5 fromat
/// in the file case.h5
void MGEquationsSystem::movemesh() {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::movemesh \n  " );)

    std::cout << "Now we try to move mesh" << std::endl;
    MGSolDA* mgsoldsx = get_eqs ( "SDSX" );
    MGSolDA* mgsoldsy = get_eqs ( "SDSY" );
    MGSolDA* mgsoldsz =nullptr;
if(_mgmesh._dim == 3) mgsoldsz = get_eqs ( "SDSZ" );

    //   (mgsolt->x_old[0][NoLevels-1])->localize(*mgsolt ->x_oold[NoLevels-1]);
    const int flag_moving_mesh = ( int ) ( _mgutils._dict_geom["moving_mesh"] );
    const int NoLevels = ( int ) ( _mgutils._dict_geom["nolevels"] );
//     const int n_nodes = _mgmesh._NoNodes[NoLevels - 1];
//     const int n_elem = _mgmesh._NoElements[0][NoLevels - 1];

    if ( flag_moving_mesh ) {
        /// E) mesh update
        const int n_nodes = _mgmesh._NoNodes[NoLevels - 1];
        int offsetp = 0 * n_nodes;
        for ( int inode = 0; inode < n_nodes; inode++ ) {
            offsetp = 0 * n_nodes;
            double disp = ( *mgsoldsx->x_old[0][NoLevels - 1] ) ( inode ) - ( *mgsoldsx->x_old[1][NoLevels - 1] ) ( inode );
            _mgmesh._xyz[inode + offsetp] =
                _mgmesh._xyzo[inode + offsetp] + ( *mgsoldsx->x_old[0][NoLevels - 1] ) ( inode );
            _mgmesh._dxdydz[inode + offsetp] = disp;
            offsetp = 1 * n_nodes;
            disp = ( *mgsoldsy->x_old[0][NoLevels - 1] ) ( inode ) - ( *mgsoldsy->x_old[1][NoLevels - 1] ) ( inode );
            //             _mgmesh._xyz[inode+offsetp] += disp; cerroni
            _mgmesh._dxdydz[inode + offsetp] = disp;
            _mgmesh._xyz[inode + offsetp] =
                _mgmesh._xyzo[inode + offsetp] + ( *mgsoldsy->x_old[0][NoLevels - 1] ) ( inode );
if(_mgmesh._dim == 3){
            offsetp = 2 * n_nodes;
            disp = ( *mgsoldsz->x_old[0][NoLevels - 1] ) ( inode ) - ( *mgsoldsz->x_old[1][NoLevels - 1] ) ( inode );
            _mgmesh._dxdydz[inode + offsetp] = disp;
            _mgmesh._xyz[inode + offsetp] =
                _mgmesh._xyzo[inode + offsetp] + ( *mgsoldsz->x_old[0][NoLevels - 1] ) ( inode );
}
        }
    }
}
// =============================================================================
void MGEquationsSystem::init ( const std::vector<FIELDS>& /*pbName*/ ) {

    TRACKING_FUN(printf ( " --> MGEquationsSystem.C: MGEquationsSystem::init \n  " );)

    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
//         MGSolDA* mgsol = eqn->second;  // get the pointer
        //     mgsol -> set_ext_fields(pbName);          // init ext fields
    }

    TRACKING_FUN(printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::init \n  " );)
    return;
}
// =============================================================================
#ifdef TWO_PHASE
void MGEquationsSystem::set_mgcc ( MGSolCC& cc ) {

    TRACKING_FUN( printf ( " --> MGEquationsSystem.C: MGEquationsSystem::set_mgcc \n" );)
    // Reading operators
    for ( iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++ ) {
        MGSolDA* mgsol = eqn->second;  // get the pointer
        mgsol->set_mgcc ( cc );          // set mgcc
    }

    TRACKING_FUN( printf ( " <-- MGEquationsSystem.C: MGEquationsSystem::set_mgcc \n" ); )
    return;
}
#endif
// kate: indent-mode cstyle; replace-tabs on;



