#ifndef _printinfo_cfg
#define _printinfo_cfg

// ===========================

//  #define DEBUG2


#ifdef DEBUG2
#define DEBUG2(x) x
#else 
#define DEBUG2(x)
#endif

// ---------------------------------------------------------------------------
#define PRINT_INFO (1)
// PRINT_INFO (option=0) routine info
#ifdef PRINT_INFO
#define P_INFO(x) x
#else 
#define P_INFO(x)
#endif
// ---------------------------------------------------------------------------

#define PRINT_TIME (1)
// PRINT_INFO (option=1) time info
// PRINT_INFO (option=0) no time info
#ifdef PRINT_TIME
#define TIME(x) x
#else 
#define TIME(x)
#endif

// ---------------------------------------------------------------------------

  #define PRINT_CONV (1)
// PRINT_INFO (option=1) time info
#ifdef PRINT_CONV
#define P_CONV(x)  x
#else 
#define P_CONV(x)
#endif
// ---------------------------------------------------------------------------

#define PRINT_PROC (1)
#ifdef PRINT_PROC
#define P_PROC(x) x
#else 
#define P_PROC(x)
#endif

// ---------------------------------------------------------------------------
// #define PRINT_MED (1)
// print the field and interfaces passed between problems with MED interface

// ---------------------------------------------------------------------------
//    #define TRACK_FUN

#ifdef TRACK_FUN
#define TRACKING_FUN(x) x
#else 
#define TRACKING_FUN(x)
#endif
// ========= No user   ==================



// ============== no user ==================
#if PRINT_TIME==1
#include <ctime>
#endif



//*****************
#ifdef LM_REFCOUNT

//refcount needs init (and also debug mode, by the way)
    #ifndef LM_INIT   
      #define LM_INIT
    #endif

//also,if refcount starts then LM_REAL must also start, otherwise
                      //you get ambiguous references to Real
    #ifndef LM_REAL 
    #define LM_REAL
    #endif
    
  #include "reference_counted_object.h"

#endif
 //******************


#endif // ------end file ---------------------
