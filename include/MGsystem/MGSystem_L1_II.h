#ifndef __mgequationsmap_add_solver_h__
#define __mgequationsmap_add_solver_h__

// -----------------------------------------------------------------------------------
// ================================================================================================
template <class SolverClass>
void MGEquationsSystem::AddSolver(
    std::string SystemName, ///< system name
    int nSys, 
    int nPieceWise,    int nLinear, int nQuadratic, 
    std::string VarName     ///< variable name
) {// =============================================================================================
  int nvars_in[3]={nPieceWise,nLinear,nQuadratic};
  SolverClass* mgs = new SolverClass(*this, nvars_in, SystemName, VarName);
  set_eqs(mgs);
  set_num_eqs(mgs->_eqname_D, nSys);
  return;
}
// ================================================================================================
template <class SolverClass>
void MGEquationsSystem::AddSolver(
    std::string SystemName, ///< system name
    int nSys, 
    int nPieceWise, int nLinear, int nQuadratic, 
    std::string VarName,          ///< variable name
    std::vector<FIELDS> PBname
) {// ==============================================================================================
  int nvars_in[3]={nPieceWise,nLinear,nQuadratic};
  SolverClass* mgs = new SolverClass(*this, nvars_in, SystemName, VarName);
  set_eqs(mgs);
  set_num_eqs(mgs->_eqname, nSys);
  mgs->set_ext_fields(PBname);
  return;
}

#endif
