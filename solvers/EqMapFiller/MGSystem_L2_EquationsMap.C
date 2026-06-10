#include <fstream>
#include <vector>
#include <string_view>
#include <functional>
// #include "EquationsMap.h"
// #include "Equations_conf.h"
// #include "Solverlib_conf.h"

#include "MGSystem_L2.h"
#include "MGSolver_L0_Base.h"
#include "MGSolver_L2_DA.h"


// =================================================================================================
//  EQUATION TABLE
// =================================================================================================
//  FIELDS -> (*)_F -> init(*)_F(EqSys)   
//   init(*)_F(EqSys) mus be defined in (*)_FEquationMap.C file in PLAT_USER/application
// Initialization of problem classes
 // ===========================================================================
 //   Initialization  dummy functions ------------------------------------------------------------------>
 // ===========================================================================
// ================================== NS_EQUATIONS ================================
// if(FIELDS == NS_F || FIELDS == NSX_F || FIELDS == NSY_F || FIELDS == NSZ_F)    
// initNavierStokes(EqSys)  --> constructor;
// couplingNavierStokes  ProbObj<-- add equation to the problem
void EquationsMap::initNavierStokes(EquationSystemsExtendedM& EqSys) {std::cout<<"initNavierStokes is not defined!!!";  return; }
void EquationsMap::couplingNavierStokes(MGSolDA*& ProbObj,int EqnName ){std::cout<<"couplingNavierStokes is not defined!!!!!";  return;}; /// Initialization of Navier Stokes system 
  
// ================================== FSI_EQUATIONS ================================
void EquationsMap::initFluidStructure(EquationSystemsExtendedM& EqSys) {std::cout<<"initFluidStructure is not defined";  return;  }
    // if(FIELDS == FS_F)                                               initFluidStructure(EqSys);
// ================================== DS_EQUATIONS ================================
  void EquationsMap::initDisplacements(EquationSystemsExtendedM& EqSys) {std::cout<<"initDisplacements is not defined";  return; }    
    // if(FIELDS == SDS_F)                                              initDisplacements(EqSys);
// ================================== B_EQUATIONS ================================
  void EquationsMap::initBoundary(EquationSystemsExtendedM& EqSys) {std::cout<<"MGSolverDA not defined";   return; }    
    // if(FIELDS == B_F || FIELDS == BX_F || FIELDS == BY_F ||  FIELDS == BZ_F)        initBoundary(EqSys);
// ================================== T_EQUATIONS ================================
  void EquationsMap::initTemperature(EquationSystemsExtendedM& EqSys) {   return; }    
    // if(FIELDS == T_F)                                                initTemperature(EqSys);
// ================================== DA_EQUATIONS ================================
void EquationsMap::initGenericDA(EquationSystemsExtendedM& EqSys) {std::cout<<"MGSolverDA not defined";return;}   
    // if(FIELDS == DA_F)                                               initGenericDA(EqSys);
// ================================== TBK_EQUATIONS ================================
  void EquationsMap::initDynamicTurbulence(EquationSystemsExtendedM& EqSys) {  return; }    
    // if(FIELDS == K_F || FIELDS == EW_F)                                   initDynamicTurbulence(EqSys);
// ================================== TTBK_EQUATIONS ================================
  void EquationsMap::initThermalTurbulence(EquationSystemsExtendedM& EqSys) {   return;  }    
    // if(FIELDS == KTT_F || FIELDS == EWTT_F)      initThermalTurbulence(EqSys);
    
// TWO-PHASE FLOW ----------------------------------------------------->    
    
// ================================== COLOR_EQUATIONS ================================
void EquationsMap::initColor(EquationSystemsExtendedM& EqSys) {   return;  }   
    // if(FIELDS == CO_F)                                               initColor(EqSys); 
// ================================== CURVATURE_EQUATIONS ================================
  void EquationsMap::initCurvature(EquationSystemsExtendedM& EqSys) {   return;  }
    // if(FIELDS == KK_F)                                               initCurvature(EqSys);
 
//  AdJOINT CONTROL  OPTIMAL ------------------------------------------------->
// ================================== ADJ_EQUATIONS ================================
  void EquationsMap::initAdjointFSI(EquationSystemsExtendedM& EqSys) { return; }   
    // if(FIELDS == NSA_F || FIELDS == NSAX_F || FIELDS == NSAY_F || FIELDS == NSAZ_F) initAdjointNavierStokes(EqSys);
// ================================== NSA_EQUATIONS ================================
  void EquationsMap::initAdjointNavierStokes(EquationSystemsExtendedM& EqSys) {  return; }    
    // if(FIELDS == FSA_F)                                              initAdjointFSI(EqSys);
 // ================================== TA_EQUATIONS ================================
  void EquationsMap::initAdjointTemperature(EquationSystemsExtendedM& EqSys) { return;  }
    // if(FIELDS == TA_F)                                               initAdjointTemperature(EqSys);
 // ================================== TBKA_EQUATIONS ================================
  void EquationsMap::initAdjointDynamicTurbulence(EquationSystemsExtendedM& EqSys) {  return; }    
    // if(FIELDS == KA_F || FIELDS == EWA_F)                                 initAdjointDynamicTurbulence(EqSys);
  // ================================== CTRL_EQUATIONS ================================
  void EquationsMap::initControlTemperature(EquationSystemsExtendedM& EqSys) {  return;  }  
    // if(FIELDS == CTRL_F)                                             initControlTemperature(EqSys);

// POROUS MEDIA EQS ------------------------------------------------>    
  void EquationsMap::initImmersedBoundary(EquationSystemsExtendedM& EqSys) {    return;  }
  // if(FIELDS == IB_F)                                               initImmersedBoundary(EqSys);
  
// COUPLING EQ ----------------------------------------------->    
    // ================================== T_COUPLED_EQUATIONS ================================
  void EquationsMap::initCoupledTemperature(EquationSystemsExtendedM& EqSys) {  return; } 
    // if(FIELDS == TCOUP_F)                                            initCoupledTemperature(EqSys);
    
// AUXILIARY EQS ------------------------------------------>    
void EquationsMap::initAreaMono(EquationSystemsExtendedM& EqSys) {   return; }
    // if(FIELDS == A_F)                                                initAreaMono(EqSys);
void EquationsMap::initFlowrateMono(EquationSystemsExtendedM& EqSys) {    return;  } 
    // if(FIELDS == Q_F)                                                initFlowrateMono(EqSys);
void EquationsMap::initReynoldsStressTensor(EquationSystemsExtendedM& EqSys) {    return;  }
    // if(FIELDS == TAU_F || FIELDS==TAUXX_F || FIELDS==TAUXY_F || FIELDS==TAUYY_F)    initReynoldsStressTensor(EqSys);
void EquationsMap::initReynoldsHeatFlux(EquationSystemsExtendedM& EqSys) {    return;  }        
 // if(FIELDS == THF_F || FIELDS==THFX_F || FIELDS==THFY_F || FIELDS == THFZ_F)     initReynoldsHeatFlux(EqSys);
 
 // ===========================================================================
 // Activate dummy functions ------------------------------------------------------------------>
 // ===========================================================================
    
  
  void EquationsMap::couplingStructuralMechanics(MGSolDA*& ProbObj,int EqnName ){std::cout<<"couplingStructuralMechanics is not defined!!!!!";  return;}; /// Initializaiton of Navier Stokes system  
  void EquationsMap::couplingFluidStructure(MGSolDA*& ProbObj,int EqnName ){std::cout<<"couplingFluidStructure is not defined!!!!!";  return;}; /// Initializaiton of Fluid Structure Interaction system
  void EquationsMap::couplingBoundary(MGSolDA*& ProbObj,int EqnName ){std::cout<<"couplingBoundary is not defined!!!!!";  return;}; /// Initializaiton of boundary system
  void EquationsMap::couplingDisplacement(MGSolDA*& ProbObj,int EqnName ){std::cout<<"couplingDisplacement is not defined!!!!!";  return;}; /// Initializaiton of solid displacement system
  void EquationsMap::couplingTemperature(MGSolDA*& ProbObj,int EqnName ){std::cout<<"couplingTemperature is not defined!!!!!";  return;};/// Initializaiton of Temperature equation
  void EquationsMap::couplingGenericDA(MGSolDA*& ProbObj,int EqnName ){std::cout<<"couplingGenericDA is not defined!!!!!";  return;};  /// Initializaiton of DA equation

 
  void EquationsMap::couplingDynamicTurbulence(MGSolDA*& ProbObj,int EqnName ){}; /// Initializaiton of dynamical turbulence system
  void EquationsMap::couplingThermalTurbulence(MGSolDA*& ProbObj,int EqnName ){};/// Initializaiton of thermal turbulence

 // two-phase
  void EquationsMap::couplingColor(MGSolDA*& ProbObj,int EqnName ){}; /// Initializaiton of color-curvature system
  void EquationsMap::couplingCurvature(MGSolDA*& ProbObj,int EqnName ){};/// Initializaiton of color-curvature system

  // Control optimal
  void EquationsMap::couplingAdjointNavierStokes(MGSolDA*& ProbObj,int EqnName ){};/// Initializaiton of Adjoint Navier Stokes system
  void EquationsMap::couplingAdjointFSI(MGSolDA*& ProbObj,int EqnName ){}; /// Initializaiton of Adjoint Fluid Structure Interaction system
  void EquationsMap::couplingAdjointTemperature(MGSolDA*& ProbObj,int EqnName ){};  /// Initializaiton of Adjoint temperature equation
  void EquationsMap::couplingAdjointDynamicTurbulence(MGSolDA*& ProbObj,int EqnName ){}; /// Initializaiton of Adjoint dynamical turbulence system

  void EquationsMap::couplingControlTemperature(MGSolDA*& ProbObj,int EqnName ){};  /// Initializaiton of Control temperature system
  void EquationsMap::couplingCoupledTemperature(MGSolDA*& ProbObj,int EqnName ){}; /// Initializaiton of Temperature control optimality system
  void EquationsMap::couplingImmersedBoundary(MGSolDA*& ProbObj,int EqnName ){}; /// Initializaiton of Immersed Boundary system

  void EquationsMap::couplingAreaMono(MGSolDA*& ProbObj,int EqnName ){};  /// Initializaiton of one-dimensional Area system
  void EquationsMap::couplingFlowrateMono(MGSolDA*& ProbObj,int EqnName ){}; /// Initializaiton of one-dimensional Flowrate system
  void EquationsMap::couplingReynoldsStressTensor(MGSolDA*& ProbObj,int EqnName ){}; /// Initialization of Reynolds Stress Tensor
  void EquationsMap::couplingReynoldsHeatFlux(MGSolDA*& ProbObj,int EqnName ){};
// =================================================================================================
//  EQUATION TABLE
// =================================================================================================
















  




//======================================================================
/// EquationsMap standard constructor
EquationsMap::EquationsMap(std::string_view app_dir) {
  TRACKING_FUN(printf(" --> EquationsMap.C:EquationsMap(constructor) \n");)
  Fill_FIELD_map();
  ReadEquationsToAdd( app_dir);
  Fill_pbName();
  TRACKING_FUN(printf(" <--- EquationsMap.C:EquationsMap(constructor) \n");)
}


//======================================================================
/// Function that handles #EquationsMap::_myproblemP vector filling
void EquationsMap::Fill_pbName() {
  TRACKING_FUN(printf(" -> EquationsMap.C:Fill_pbName \n");)
  _myproblemP.clear();
#ifdef _TURBULENCE_
  ReadTurbulenceInfo();
#endif

  int i = 0;
  std::cout << "------------------------------------------------------------------------------ \n";
  std::cout << " EquationMap::Fill_pbName() \n";

  for(std::map<std::string, int>::iterator it = _EquationsToAdd.begin(); it != _EquationsToAdd.end(); ++it) {
    if((it->first).compare(0, 3, "MG_") == 0) {
      if(it->second != 0) {
        FIELDS ff = _map_str2field[it->first];
        _myproblemP.push_back(ff);
        std::cout << " Adding field " << it->first << " with Equation tab number " << _myproblemP[i] << "\n";
        i++;
        // -------------------------------------------------------------------------------------------
#ifdef _TURBULENCE_
        if(it->first == "MG_DynamicalTurbulence") {
          _myproblemP.push_back(DIST);
          std::cout << " Adding wall distance field with Equation tab number " << _myproblemP[i] << "\n";
          i++;
          _myproblemP.push_back(MU_T);
          std::cout << " Adding dynamical turbulence field with Equation tab number " << _myproblemP[i]
                    << "\n";
          i++;
        }
        if(it->first == "MG_ThermalTurbulence") {
          _myproblemP.push_back(ALPHA_T);
          std::cout << " Adding thermal turbulence with Equation tab number " << _myproblemP[i] << "\n";
          i++;
        }
#endif
        // -------------------------------------------------------------------------------------------
      }
    }
  }
  TRACKING_FUN(printf(" <- EquationsMap.C:Fill_pbName \n");)
// ==================================================

  return;
}
 // void EquationsMap::SolverMapNS_F(){};
// void initSolverNS(EquationSystemsExtendedM &EqMap){}; /// Initializaiton of Navier Stokes system 
// ==============================================================================================
/// This function adds equations to #FEMUS::_mg_equations_map
void EquationsMap::FillEquationMap(
EquationSystemsExtendedM&
EqSys
) {  // ==========================================================================
// TRACKING_FUN(printf(" -> EquationsMap.C: EquationsMap::FillEquationMap \n");)
// ==============================================================================================
  int n_equations = _myproblemP.size();
   // void initSolverNS
  // FredMemFn p = &Fred::f;
 
  // InitSolver p = &std::function<void(EquationSystemsExtendedM &EqMap)>(initSolverNS);
  // InitSolverMap[NS_F]= p;
  // InitSolver InitSolverMap[NS_F]=&EquationsMap::initSolverNS;
// SolverMap[NS_F]=SolverMapNS_F;
  // Initialization of problem classes
  for(int iname = 0; iname < n_equations; iname++) {
    FIELDS p = _myproblemP[iname];
    if(p == NS_F || p == NSX_F || p == NSY_F || p == NSZ_F)    initNavierStokes(EqSys);
    if(p == FS_F)                                               initFluidStructure(EqSys);
    if(p == SDS_F)                                              initDisplacements(EqSys);
    if(p == B_F || p == BX_F || p == BY_F ||  p == BZ_F)        initBoundary(EqSys);
    if(p == T_F)                                                initTemperature(EqSys);
    if(p == DA_F)                                               initGenericDA(EqSys);
    if(p == K_F || p == EW_F)                                   initDynamicTurbulence(EqSys);
    if(p == KTT_F || p == EWTT_F)                               initThermalTurbulence(EqSys);
    if(p == CO_F)                                               initColor(EqSys); 
    if(p == KK_F)                                               initCurvature(EqSys);
    if(p == NSA_F || p == NSAX_F || p == NSAY_F || p == NSAZ_F) initAdjointNavierStokes(EqSys);
    if(p == FSA_F)                                              initAdjointFSI(EqSys);
    if(p == TA_F)                                               initAdjointTemperature(EqSys);
    if(p == KA_F || p == EWA_F)                                 initAdjointDynamicTurbulence(EqSys);
    if(p == CTRL_F)                                             initControlTemperature(EqSys);
    if(p == IB_F)                                               initImmersedBoundary(EqSys);
    if(p == TCOUP_F)                                            initCoupledTemperature(EqSys);
    if(p == A_F)                                                initAreaMono(EqSys);
    if(p == Q_F)                                                initFlowrateMono(EqSys);
    if(p == TAU_F || p==TAUXX_F || p==TAUXY_F || p==TAUYY_F)    initReynoldsStressTensor(EqSys);
    if(p == THF_F || p==THFX_F || p==THFY_F || p == THFZ_F)     initReynoldsHeatFlux(EqSys);
    
  /*      switch(p){   
    case NS_F: case NSX_F: case NSY_F: case NSZ_F:    initNavierStokes(EqSys); break;   
    case FS_F:                                               initFluidStructure(EqSys);break;   
    case SDS_F:                                              initDisplacements(EqSys); break;   
    case B_F : case BX_F: case BY_F:  case BZ_F:        initBoundary(EqSys);   break; }
    case T_F:                                                initTemperature(EqSys);  break;
    case DA_F:                                               initGenericDA(EqSys);  break; // This
    case K_F: case EW_F:                                   initDynamicTurbulence(EqSys);break; for
    case KTT_F: case EWTT_F:                               initThermalTurbulence(EqSys); break; */  
    
  }

  // This Function calls the MGSolDA::init_ext_fields()
  for(auto eqn = EqSys._equations.begin(); eqn != EqSys._equations.end(); eqn++) {
    MGSolDA* mgsol = eqn->second;  // get the pointer
    mgsol->setUpExtFieldData();
    setProblems(mgsol);
  }
  // ==============================================================================================
  TRACKING_FUN(printf(" <- EquationsMap.C:FillEquationMap \n");)

  return;
}

/// This function fills #EquationsMap::_map_str2field map
void EquationsMap::Fill_FIELD_map() {  // map equation   map_str2field
  TRACKING_FUN(printf(" -> EquationsMap.C:Fill_FIELD_map \n");)
// ==============================================================================================
  _map_str2field["MG_NavierStokes"] = NS_F;
  _map_str2field["NS_F"] = NS_F;    // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["NSX_F"] = NSX_F;  // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["NSY_F"] = NSY_F;  // [1] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["NSZ_F"] = NSZ_F;  // [2] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["MG_FluidStructure"] = FS_F;
  _map_str2field["FS_F"] = FS_F;    // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["FSX_F"] = FSX_F;  // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["FSY_F"] = FSY_F;  // [1] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["FSZ_F"] = FSZ_F;  // [2] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["MG_StructuralMechanics"] = SM_F;
  _map_str2field["SM_F"] = SM_F;    // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["SMX_F"] = SMX_F;  // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["SMY_F"] = SMY_F;  // [1] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["SMZ_F"] = SMZ_F;  // [2] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
  _map_str2field["MG_Pressure"] = P_F;
  _map_str2field["P_F"] = P_F;  // [3] -> Pressure (linear (1);NS_EQUATIONS==0 or 2)
  _map_str2field["MG_Temperature"] = T_F;
  _map_str2field["T_F"] = T_F;  // [4] -> Temperature   (quadratic (2);T_EQUATIONS)
  _map_str2field["MG_DynamicalTurbulence"] = K_F;
  _map_str2field["K_F"] = K_F;    // [5] -> Turbulence K  (quadratic (2);TB K_EQUATIONS)
  _map_str2field["EW_F"] = EW_F;  // [6] -> Turbulence W  (quadratic (2);TB W_EQUATIONS)
  _map_str2field["MG_ThermalTurbulence"] = KTT_F;
  _map_str2field["KTT_F"] = KTT_F;    // [7] -> Turbulence K  (quadratic (2);TB K_EQUATIONS)
  _map_str2field["EWTT_F"] = EWTT_F;  // [8] -> Turbulence W  (quadratic (2);TB W_EQUATIONS)
  _map_str2field["MG_Displacement"] = SDS_F;
  _map_str2field["SDS_F"] = SDS_F;    // [9] -> Displacement (quadratic (2); DS_EQUATIONS)
  _map_str2field["SDSX_F"] = SDSX_F;  // [9] -> Displacement (quadratic (2); DS_EQUATIONS)
  _map_str2field["SDSY_F"] = SDSY_F;  // [10]-> Displacement (quadratic (2); DS_EQUATIONS)
  _map_str2field["SDSZ_F"] = SDSZ_F;  // [11]-> Displacement (quadratic (2); DS_EQUATIONS)
  _map_str2field["MG_Boundary"] = B_F;
  _map_str2field["B_F"] = B_F;    // [9] -> Displacement (quadratic (2); DS_EQUATIONS)
  _map_str2field["BX_F"] = BX_F;  // [9] -> Displacement (quadratic (2); DS_EQUATIONS)
  _map_str2field["BY_F"] = BY_F;  // [10]-> Displacement (quadratic (2); DS_EQUATIONS)
  _map_str2field["BZ_F"] = BZ_F;  // [11]-> Displacement (quadratic (2); DS_EQUATIONS)
  _map_str2field["BP_F"] = BP_F;  // [11]-> Displacement (quadratic (2); DS_EQUATIONS)
  _map_str2field["MG_DA"] = DA_F;
  _map_str2field["DA_F"] = DA_F;  // [12]-> DA solver (quadratic (2); DA_EQUATIONS)
  // two-phase
  _map_str2field["MG_ColorFunction"] = CO_F;
  _map_str2field["CO_F"] = CO_F;
  _map_str2field["MG_Curvature"] = KK_F;
  _map_str2field["KK_F"]  = KK_F;
  _map_str2field["MG_ExtendedCurvature"] = EK_F;
  _map_str2field["EK_F"]  = EK_F;
  _map_str2field["MG_TwoFluidModel"] = TFM_F;
  _map_str2field["TFM_F"] = TFM_F;
  _map_str2field["NS1_F"] = NS1_F;
  _map_str2field["NS1X_F"] = NS1X_F;
  _map_str2field["NS1Y_F"] = NS1Y_F;
  _map_str2field["NS1Z_F"] = NS1Z_F;
  _map_str2field["NS2_F"] = NS2_F;
  _map_str2field["NS2X_F"] = NS2X_F;
  _map_str2field["NS2Y_F"] = NS2Y_F;
  _map_str2field["NS2Z_F"] = NS2Z_F;
  _map_str2field["P_TFM_F"] = PTF_F;
  // adjoint
  _map_str2field["MG_AdjointNavierStokes"] = NSA_F;
  _map_str2field["NSA_F"] = NSA_F;    // [15]-> adjoint NS or FSI equations (Dimension)
  _map_str2field["NSAX_F"] = NSAX_F;  // [15]-> adjoint NS or FSI equations (Dimension)
  _map_str2field["NSAY_F"] = NSAY_F;  // [15]-> adjoint NS or FSI equations (Dimension)
  _map_str2field["NSAZ_F"] = NSAZ_F;  // [15]-> adjoint NS or FSI equations (Dimension)
  _map_str2field["MG_AdjointFluidStructure"] = FSA_F;
  _map_str2field["FSA_F"] = FSA_F;    // [15]-> adjoint NS or FSI equations (Dimension)
  _map_str2field["FSAX_F"] = FSAX_F;  // [15]-> adjoint NS or FSI equations (Dimension)
  _map_str2field["FSAY_F"] = FSAY_F;  // [15]-> adjoint NS or FSI equations (Dimension)
  _map_str2field["FSAZ_F"] = FSAZ_F;  // [15]-> adjoint NS or FSI equations (Dimension)
  _map_str2field["MG_AdjointDA"] = DA_P;
  _map_str2field["DA_P"] = DA_P;  // [13]-> DA solver (piecewise; DA_EQUATIONS)
  _map_str2field["MG_AdjointTemperature"] = TA_F;
  _map_str2field["TA_F"] = TA_F;  // [14]-> Temp adjoint
  _map_str2field["MG_AdjointTurbulence"] = KA_F;
  _map_str2field["KA_F"] = KA_F;    // [18] + 19 Adjoint turbulence
  _map_str2field["EWA_F"] = EWA_F;  // [18] + 19 Adjoint turbulence
  // control
  _map_str2field["MG_Laplacian"] = CO_F;
  _map_str2field["CO_F"] = CO_F;  // [16]-> Color function for FSI equations
  _map_str2field["MG_ControlTemperature"] = CTRL_F;
  _map_str2field["CTRL_F"] = CTRL_F;                  // [16]-> Color function for FSI equations
  _map_str2field["CTRLX_F"] = CTRLX_F;                // [16]-> Color function for FSI equations
  _map_str2field["CTRLY_F"] = CTRLY_F;                // [16]-> Color function for FSI equations
  _map_str2field["CTRLZ_F"] = CTRLZ_F;                // [16]-> Color function for FSI
  _map_str2field["MG_ImmersedBoundary"] = IB_F;       // [16]-> Color function for FSI
  _map_str2field["MG_CoupledTemperature"] = TCOUP_F;  // [29]-> Coupled Temperature optimality system
  _map_str2field["TCOUP_F"] = TCOUP_F;                // [4] -> Temperature   (quadratic (2);T_EQUATIONS)
  _map_str2field["MG_MonoArea"] = A_F;                // [30]-> Area for monodimensional code
  _map_str2field["MG_MonoFlowrate"] = Q_F;            // [31]-> Flowrate for monodimensional code
  _map_str2field["MG_ReynoldsStressTensor"] = TAU_F;  // [32] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
  _map_str2field["TAU_F"] = TAU_F;                    // [32] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
  _map_str2field["TAUXX_F"] = TAUXX_F;                // [32] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
  _map_str2field["TAUXY_F"] = TAUXY_F;                // [33] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
  _map_str2field["TAUYY_F"] = TAUYY_F;                // [34] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
  _map_str2field["THF_F"] = THF_F;                    // [32] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
  _map_str2field["MG_ReynoldsHeatFlux"] = THF_F;      // [32] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
  _map_str2field["THFX_F"] = THFX_F;                  // [32] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
  _map_str2field["THFY_F"] = THFY_F;                  // [33] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
  _map_str2field["THFZ_F"] = THFZ_F;                  // [34] -> Anisotropy (quadratic (2); TAU_EQUATIONS)
// ==============================================================================================
  TRACKING_FUN(printf(" <- EquationsMap.C:Fill_FIELD_map \n");)
  return;
}


// ====================================================================================
void EquationsMap::ReadEquationsToAdd(std::string_view app_dir) {
  TRACKING_FUN(printf(" -> EquationsMap.C:ReadEquationsToAdd \n");)
  std::ostringstream filename; filename << getenv("APP_PATH") << "/DATA/Equations.in";
  std::ifstream fin; fin.open(filename.str().c_str());  // stream file

  std::string buf = "";
  int value;
  if(fin.is_open()) {   // ==================================================-------
    while(buf != "/") {
      fin >> buf;  // find "/" file start
    }
    fin >> buf;
    while(buf != "/") {
      if(buf == "#") {
        getline(fin, buf);  // comment line
      }
      else {
        fin >> value;
        _EquationsToAdd.insert(std::pair<std::string, int>(buf, value));
      }
      fin >> buf;
    }
  }  // ==================================================--------------------------------
  TRACKING_FUN(printf(" <- EquationsMap.C:ReadEquationsToAdd \n");)
  return;
}

// ==============================================================================================
/// Solution sharing between couplingd equations
void EquationsMap::setProblems(
MGSolDA*& ProbObj
) { 
  TRACKING_FUN(printf(" -> EquationsMap.setProblems \n");)
// ==============================================================================================
  _PieceEq = _LinearEq = _QuadEq = 0;
  std::cout << "\n------------------------------------------------\n";
  std::cout << "  Activating fields for problem " << ProbObj->_eqname_D;

  for(std::map<std::string, int>::iterator eq = _EquationsToAdd.begin(); eq != _EquationsToAdd.end(); ++eq) {
    std::string EqnName = eq->first;
    int EqnLabel = eq->second;

    if(EqnLabel > 0) {
      
      // DA --------------------------------------------------------------------------------------
      if(EqnName == "MG_DA") {  couplingGenericDA(ProbObj,EqnLabel);}
      // NAVIER-STOKES ----------------------------------------------------------------------------
      if(EqnName == "MG_NavierStokes") {  couplingNavierStokes(ProbObj,EqnLabel);}
      // } // ----------------------------------------------------------------------------------------
      if(EqnName == "MG_AdjointNavierStokes"){couplingAdjointNavierStokes(ProbObj,EqnLabel); } 
      // ADJOINT NAVIER STOKES ------------------------
        // ----------------------------------------------------------------------------------------
      // FLUID STRUCTURE
      if(EqnName == "MG_FluidStructure") {  couplingFluidStructure(ProbObj,EqnLabel);}
      // ADJOINT FLUID STRUCTURE
      if(EqnName == "MG_AdjointFluidStructure") { couplingAdjointFSI(ProbObj,EqnLabel);}
      // STRUCTURAL MECHANICS
      if(EqnName == "MG_StructuralMechanics") {couplingStructuralMechanics(ProbObj,EqnLabel);}
      // SOLID DISPLACEMENTS
      if(EqnName == "MG_Displacement") { couplingDisplacement(ProbObj,EqnLabel);}
      // BOUNDARY
      if(EqnName == "MG_Boundary") { couplingBoundary(ProbObj,EqnLabel);}
      // TEMPERATURE
      if(EqnName == "MG_Temperature")  { couplingTemperature(ProbObj,EqnLabel);}
      // CONTROL TEMPERATURE
      if(EqnName == "MG_ControlTemperature") { couplingControlTemperature(ProbObj,EqnLabel);}
      // ADJOINT TEMPERATURE
      if(EqnName == "MG_AdjointTemperature") { couplingAdjointTemperature(ProbObj,EqnLabel);}
      // DYNAMICAL TURBULENCE
      if(EqnName == "MG_DynamicalTurbulence") { couplingDynamicTurbulence(ProbObj,EqnLabel);}
      // THERMAL TURBULENCE
      if(EqnName == "MG_ThermalTurbulence") { couplingThermalTurbulence(ProbObj,EqnLabel);}
      // COLOR FUNCTION
      if(EqnName == "MG_ColorFunction") { couplingColor(ProbObj,EqnLabel);}
      if(EqnName == "MG_Curvature")  { couplingCurvature(ProbObj,EqnLabel);}
      // IMMERSED BOUNDARY
      if(EqnName == "MG_ImmersedBoundary") { couplingImmersedBoundary(ProbObj,EqnLabel);}
      // ADJOINT TURBULENCE
      if(EqnName == "MG_AdjointTurbulence") { couplingAdjointDynamicTurbulence(ProbObj,EqnLabel);}
      // COUPLED OPTIMAL CONTROL TEMPERATURE
      if(EqnName == "MG_CoupledTemperature") { couplingCoupledTemperature(ProbObj,EqnLabel);}
      // ONE-DIMENSIONAL CODE
      if(EqnName == "MG_MonoArea") { couplingAreaMono(ProbObj,EqnLabel);}
      // REYNOLDS STRESS TENSOR
      if(EqnName == "MG_ReynoldsStressTensor")  { couplingReynoldsStressTensor(ProbObj,EqnLabel);}
      if(EqnName == "MG_ReynoldsHeatFlux")   { couplingReynoldsHeatFlux(ProbObj,EqnLabel);}
      }
  }
  // ==============================================================================================
  TRACKING_FUN(printf(" <- EquationsMap.setProblems \n");)
  return;
}


// ================================================================================================
void EquationsMap::ReadTurbulenceInfo() {
  TRACKING_FUN(printf(" -> EquationsMap.C:ReadTurbulenceInfo \n");)
// ================================================================================================
  std::ostringstream filename;  filename << getenv("APP_PATH") << "/DATA/Turbulence.in";
  std::ifstream fin; fin.open(filename.str().c_str());  // stream file
  std::string buf = "";  std::string value;
  if(fin.is_open()) {   // ==================================================-------
    while(buf != "/") fin >> buf;   // find "/" file start
    fin >> buf;
    while(buf != "/") {
      if(buf == "#") getline(fin, buf);   // comment line
      else {fin >> value; _TurbulenceModel.insert(std::pair<std::string, std::string>(buf, value));}
      fin >> buf;
    }
  }  // ==================================================--------------------------------
// ================================================================================================
  TRACKING_FUN(printf(" <- EquationsMap.C:ReadTurbulenceInfo \n");)
  return;
}

// // ==============================================================================================
// /// Solution sharing between couplingd equations
// void EquationsMap::setProblems(
// MGSolDA*& ProbObj
// ) { 
//   TRACKING_FUN(printf(" -> EquationsMap.setProblems \n");)
// // ==============================================================================================
//   _PieceEq = _LinearEq = _QuadEq = 0;
//   std::cout << "\n------------------------------------------------\n";
//   std::cout << "  Activating fields for problem " << ProbObj->_eqname_D;
// 
//   for(std::map<std::string, int>::iterator eq = _EquationsToAdd.begin(); eq != _EquationsToAdd.end(); ++eq) {
//     std::string EqnName = eq->first;
//     int EqnLabel = eq->second;
// 
//     if(EqnLabel > 0) {
//       
//       // DA --------------------------------------------------------------------------------------
//       if(EqnName == "MG_DA") {  couplingGenericDA(ProbObj,EqnLabel);}
//       //   // setSolverMG_DA(int EqnName){
//       //   if(EqnLabel%10 >0)   ProbObj->ActivateDA(2, DA_F, "DA1", _QuadEq);
//       //   if(EqnLabel/10 >0)    ProbObj->ActivateScalar(1, DA_F, "DA1", _LinearEq);
//       //   if(EqnLabel/100 >0)    ProbObj->ActivateScalar(0, DA_F, "DA1", _PieceEq);
//       //   // }
//       // }
//       
//       
//       // NAVIER-STOKES ----------------------------------------------------------------------------
//       if(EqnName == "MG_NavierStokes") {  couplingNavierStokes(ProbObj,EqnLabel);}
// //         if(EqnLabel == 1) {
// //           ProbObj->ActivateVectField(2, NS_F, "NS0", _QuadEq,1);
// //           ProbObj->ActivateScalar(1, P_F, "NS0", _LinearEq);
// //         }
// //         if(EqnLabel == 2) {
// //           ProbObj->ActivateVectField(2, NS_F, "NS0", _QuadEq,0);
// //           ProbObj->ActivateScalar(1, P_F, "NS2P", _LinearEq);
// //         }
// //         if(EqnLabel == 3) {
// //           ProbObj->ActivateVectField(2, NS_F, "NS0", _QuadEq,0);
// //           ProbObj->ActivateScalar(1, P_F, "NS2P", _LinearEq);
// //           ProbObj->ActivateScalar(0, P_F, "NS2P", _PieceEq);
// //           
// //         }
// //  if(EqnLabel == 4) {
// //           ProbObj->ActivateVectField(2, NS_F, "NS0", _QuadEq,1);
// //           ProbObj->ActivateScalar(1, P_F, "NS0", _LinearEq);
// //           ProbObj->ActivateScalar(1, PENR_F, "NS0", _LinearEq);
// //         }
// // 
// // //          if (EqnLabel == 3) { ProbObj->ActivateScalar(2, NS_F, "NS0", _QuadEq); }
//       // } // ----------------------------------------------------------------------------------------
//       if(EqnName == "MG_AdjointNavierStokes") {   // ADJOINT NAVIER STOKES ------------------------
//          couplingAdjointNavierStokes(ProbObj,EqnLabel);
//         // if(EqnLabel <= 2) {
//         //   int coupled = (EqnLabel == 1) ? 1 : 0;
//         //   ProbObj->ActivateVectField(2, NSA_F, "NSA0", _QuadEq, coupled);
//         //   if(coupled == 0) { ProbObj->ActivateScalar(1, P_F, "NSAP", _LinearEq); }
//         // }
//       }// ----------------------------------------------------------------------------------------
//       // FLUID STRUCTURE
//       if(EqnName == "MG_FluidStructure") {  couplingFluidStructure(ProbObj,EqnLabel);}
//         // if(EqnLabel <= 2) {
//         //   int coupled = (EqnLabel == 1) ? 1 : 0;
//         //   ProbObj->ActivateVectField(2, FS_F, "FSI0", _QuadEq, coupled);
//         //   if(coupled == 0) { ProbObj->ActivateScalar(1, P_F, "FSIP", _LinearEq); }
//         // }
//       // }// ----------------------------------------------------------------------------------------
//       // ADJOINT FLUID STRUCTURE
//       if(EqnName == "MG_AdjointFluidStructure") { couplingAdjointFSI(ProbObj,EqnLabel);}
//       //   if(EqnLabel <= 2) {
//       //     int coupled = (EqnLabel == 1) ? 1 : 0;
//       //     ProbObj->ActivateVectField(2, NSA_F, "FSIA0", _QuadEq, coupled);
//       //     if(coupled == 0) { ProbObj->ActivateScalar(1, P_F, "FSIAP", _LinearEq); }
//       //   }
//       // }// ----------------------------------------------------------------------------------------
//       // STRUCTURAL MECHANICS
//       if(EqnName == "MG_StructuralMechanics") {couplingStructuralMechanics(ProbObj,EqnLabel);}
//       //   if(EqnLabel <= 2) {
//       //     int coupled = (EqnLabel == 1) ? 1 : 0;
//       //     ProbObj->ActivateVectField(2, SM_F, "SM0", _QuadEq, coupled);
//       //   }
//       // }// ----------------------------------------------------------------------------------------
//       // SOLID DISPLACEMENTS
//       if(EqnName == "MG_Displacement") { couplingDisplacement(ProbObj,EqnLabel);}
//       //   if(EqnLabel <= 2) {
//       //     int coupled = (EqnLabel == 1) ? 1 : 0;
//       //     ProbObj->ActivateVectField(2, SDS_F, "SDS", _QuadEq, coupled);
//       //   }
//       // }// ----------------------------------------------------------------------------------------
//       // BOUNDARY
//       if(EqnName == "MG_Boundary") { couplingBoundary(ProbObj,EqnLabel);}
// //         if(EqnLabel <= 2) {
// //           int coupled = (EqnLabel == 1) ? 1 : 0;
// //           ProbObj->ActivateVectField(2, B_F, "B0", _QuadEq, coupled);
// // 
// //           if(coupled == 0) {
// //             // if(NDOF_K == 1) {
// //             //   ProbObj->ActivateScalar(0, BP_F, "B2P", _PieceEq);
// //             // }
// //             // else {
// //               ProbObj->ActivateScalar(1, BP_F, "B2P", _LinearEq);
// //             // }
// //           }
// //         }
// // 
// // //         if (EqnLabel == 3) { ProbObj->ActivateScalar(2, NS_F, "NS0", _QuadEq); }
// //       }// ----------------------------------------------------------------------------------------
// 
// 
// 
// 
//       // TEMPERATURE
//       if(EqnName == "MG_Temperature")  { couplingTemperature(ProbObj,EqnLabel);}
//         // ProbObj->ActivateScalar(2, T_F, "T", _QuadEq);
//       // }// ----------------------------------------------------------------------------------------
// 
//       // CONTROL TEMPERATURE
//       if(EqnName == "MG_ControlTemperature") { couplingControlTemperature(ProbObj,EqnLabel);}
//       //   int vector;
//       //   vector = (_EquationsToAdd["MG_ControlTemperature"] == 2) ? 1 : 0;
//       //   ProbObj->ActivateControl(2, CTRL_F, "CTRL", _QuadEq, vector, 2);
//       // }// ----------------------------------------------------------------------------------------
// 
//       // ADJOINT TEMPERATURE
//       if(EqnName == "MG_AdjointTemperature") { couplingAdjointTemperature(ProbObj,EqnLabel);}
//       //   ProbObj->ActivateScalar(2, TA_F, "TA", _QuadEq);
//       // }// ----------------------------------------------------------------------------------------
// 
//       // DYNAMICAL TURBULENCE
//       if(EqnName == "MG_DynamicalTurbulence") { couplingDynamicTurbulence(ProbObj,EqnLabel);}
// //         ProbObj->ActivateCoupled(2, K_F, "K2K", _QuadEq, "K1W");
// // #ifdef _TURBULENCE_
// //         ProbObj->ActivateScalar(2, DIST, "DIST", _QuadEq);
// //         ProbObj->ActivateScalar(2, MU_T, "MU_T", _QuadEq);
// // #endif
//       // }// ----------------------------------------------------------------------------------------
// 
//       // THERMAL TURBULENCE
//       if(EqnName == "MG_ThermalTurbulence") { couplingThermalTurbulence(ProbObj,EqnLabel);}
// //         ProbObj->ActivateCoupled(2, KTT_F, "TK", _QuadEq, "TK2");
// // #ifdef _TURBULENCE_
// //         ProbObj->ActivateScalar(2, ALPHA_T, "ALPHA_T", _QuadEq);
// // #endif
// //       }// ----------------------------------------------------------------------------------------
// 
//       // COLOR FUNCTION
//       if(EqnName == "MG_ColorFunction") { couplingColor(ProbObj,EqnLabel);}
//       //   ProbObj->ActivateCoupled(2, CO_F, "C", _QuadEq, "CK"); 
//       // }
//       // CURVATURE FUNCTION
//       if(EqnName == "MG_Curvature")  { couplingCurvature(ProbObj,EqnLabel);}
//       //   ProbObj->ActivateScalar(1, KK_F, "K", _LinearEq);
//       //   ProbObj->ActivateScalar(1, EK_F, "EK", _LinearEq); 
//       // }
//       // IMMERSED BOUNDARY
//       if(EqnName == "MG_ImmersedBoundary") { couplingImmersedBoundary(ProbObj,EqnLabel);}
//       // { ProbObj->ActivateCoupled(2, IB_F, "IB1", _QuadEq, "IB2"); }
//       // ADJOINT TURBULENCE
//       if(EqnName == "MG_AdjointTurbulence") { couplingAdjointDynamicTurbulence(ProbObj,EqnLabel);}
//       // { ProbObj->ActivateCoupled(2, KA_F, "K2KA", _QuadEq, "K1WA"); }
//       // COUPLED OPTIMAL CONTROL TEMPERATURE
//       if(EqnName == "MG_CoupledTemperature") { couplingCoupledTemperature(ProbObj,EqnLabel);}
//       
// //      {
//   //      int coupled = 1; ProbObj->ActivateVectField(2, TCOUP_F, "TCOUP0", _QuadEq, coupled);
//     //    // if(coupled==0){ ProbObj->ActivateScalar ( 1, P_F, "NSAP", _LinearEq );}
//     //  }
// 
//       //
// 
//       // ONE-DIMENSIONAL CODE
//       if(EqnName == "MG_MonoArea") { couplingAreaMono(ProbObj,EqnLabel);}
//       //   ProbObj->ActivateScalar(1, A_F, "A", _LinearEq);
//       //   ProbObj->ActivateScalar(1, Q_F, "Q", _LinearEq);
//       // }
//       // REYNOLDS STRESS TENSOR
//       if(EqnName == "MG_ReynoldsStressTensor")  { couplingReynoldsStressTensor(ProbObj,EqnLabel);}
//       //   if(EqnLabel <= 2) {
//       //     int coupled = (EqnLabel == 1) ? 1 : 0;
//       //     ProbObj->ActivateVectField(2, TAU_F, "TAU", _QuadEq, coupled);
//       //   }
//       // }
//       if(EqnName == "MG_ReynoldsHeatFlux")   { couplingReynoldsHeatFlux(ProbObj,EqnLabel);}
//         // if(EqnLabel <= 2) {
//         //   int coupled = (EqnLabel == 1) ? 1 : 0;
//         //   ProbObj->ActivateVectField(2, THF_F, "THF", _QuadEq, coupled);
//         // }
//       }
//     // }
//   }
//   // ==============================================================================================
//   TRACKING_FUN(printf(" <- EquationsMap.setProblems \n");)
//   return;
// }





