#ifndef __FEMUS0__
#define __FEMUS0__

#include "MGUtils_L0.h"
#include "Solverlib_conf.h"
#include "mpi.h"
#include <vector>
#include "MGGeom_L0_El.h"


#define TRACKING_FUN(x)
#ifndef P_PROC
  #define P_PROC(x)
#endif


class MeshExtended;
class MGFemusInit;
class EquationSystemsExtendedM;
class MGSystem;
template<unsigned int DIM,GEOFamily T>
class MGGeomEl;
class MGFEMap;
class MGTimeLoop;
class BoundInterp;
class TurbUtils;
class IbUtils;

#ifdef HAVE_MED
namespace MEDCoupling {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
class DataArrayIdType;
class DataArrayDouble;
}  // namespace MEDCoupling
#endif

class FEMUS0 {

protected:
    // data communication (defined in Constructor)
//    std::unique_ptr<MGFemusInit>  _start;                         // start function

     // processors 
    bool _local_MPI_Init;  // initial mpi flag   
    
public:
    // Constructor-Destructor 
    FEMUS0() = default;
  
   virtual ~FEMUS0() = default;

};

#endif
