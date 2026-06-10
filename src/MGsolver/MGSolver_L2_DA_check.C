// std lib
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

// Femus class files ----------------->
#include "MGSolver_L2_DA.h"
// #include "MGGeom_L1_El.h"
// #include "MGMesh_L1_Extended.h"
// #include "MGSystem_L2.h"
// // #include "MGFE.h"
// #include "MGFEMap_L0.h"
// #include "MGGraph_L0.h"
// #include "MGUtils_L0.h"

void MGSolDA::check_multiphysics_D(int /*mode*/) {
    
    std::fstream fs ("check_multiphysics.txt",std::fstream::out);
    fs << " System =" << _eqname_D << " \n" << " ======================================= \n";
    fs << " Equations 0-> " << _DAdata_eq_D[2].max_neqs;
        
    for ( int deg=0; deg<3; deg++ ) {
        fs << "  \n =================================" << deg << " \n";
        fs << "  \n Coupling with =" << _DAdata_eq_D[deg].n_eqs << " equations with degree =" << deg << " \n";
        if(_DAdata_eq_D[deg].n_eqs>0) {
            fs << "_DAdata_eq_D[" << deg<<   "]:tab_eqs \n";
            for(int keq=0;keq<_DAdata_eq_D[deg].max_neqs;keq++) if(_DAdata_eq_D[deg].tab_eqs[keq]>-1) { 
                    fs <<  setw(4)<<keq<<"("<< _DAdata_eq_D[deg].tab_eqs[keq]<<")";}fs << " \n";
        }
            
        for ( int eq=0; eq<_DAdata_eq_D[deg].n_eqs; eq++ ) {
            fs << "  \n Coupling with system " << eq << "  " << _DAdata_eq_D[deg].mg_eqs[eq]->_eqname_D << ": \n";
            fs << "Pointer DASolver =" << _DAdata_eq_D[deg].mg_eqs[eq] << "\n";
            fs <<  "Index ub-data from= " << _DAdata_eq_D[deg].indx_ub[eq]*_el_dof_D[0]<<  "  to " << _DAdata_eq_D[deg].indx_ub[eq+1]*_el_dof_D[0]-1<< "\n";
        }
    }
    fs.close();
  
  return;
}
// =======================================================================
void MGSolDA::check_mtrx_rhs_D(int mode) {
      
    std::fstream fs ("check_mtrx.txt",std::fstream::out);
    if(mode<10) fs << " \n mtrx \n " << _KeM_D;  
    if(mode<20) fs << " \n rhs \n "  <<_FeM_D;      

    return;
}
// =======================================================================      
void MGSolDA::check_ass_nan_inf_D(int /*mode*/) {
    
    for(int i=0;i<_varst_D[0]*_el_dof_D[0];i++) {
        double x=_FeM_D(i);
        if(x !=x) std::cout <<  " Nan FeM ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! line=" << i; 
      
        for(int j=0;j<_varst_D[0]*_el_dof_D[0];j++) {
            double x=_KeM_D(i,j);
            if(x !=x) std::cout <<  " Nan Kem ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! line=" << i << "  " << j; 
        }
    }
    return;
}
// =======================================================================  
void MGSolDA::check_bc_D(int /*mode*/){

    std::fstream infile( "check_bc_0.txt",std::fstream::out);
    infile  << " inode  " << " " << " _bc[0][inode]  " << " " << " _bc[1][inode] " << "\n";
    for(int inode=0;inode<_Ndof_lev[_NoLevels-1];inode++){
// // printf(in" %d %d \n ",inode,_bc[1][inode]);
        infile  << inode << " " << _bc[0][inode] << " " << _bc[1][inode] << "\n";
    }
    infile.close();
 
 return;
}
