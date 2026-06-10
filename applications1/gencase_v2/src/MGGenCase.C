// std libraries ------------------------------------------
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// LibMesh includes ---------------------------------------
// #include "MGFE_conf.h"
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

// configure includes -------------------------------------
// #include "Domain_conf.h"    //  domain dimensions
#include "gencase_conf.h"      //
#include "MGGenCase.h"      //
#include "Printinfo_conf.h" //
#include "Solverlib_conf.h" //

// MED includes -------------------------------------------
#ifdef HAVE_MED
#include "MEDCoupling.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoader.hxx"
#endif

// local includes -----------------------------------------
#include "MGGenCase_conf.h"
#include "MGGeomEl.h"
#include "MGUtils.h"

// =======================================
/// Constructor
MGGenCase::MGGenCase(const Parallel::Communicator& comm_in, MGUtils& mgutils_in, MGGeomEl& geomel_in)
    : _comm(comm_in), _geomel(geomel_in), _mgutils(mgutils_in) {
  // =================================
  const int NoLevels = (int)_mgutils._geometry["nolevels"];
  _bcmat = 0;
  _dim = DIMENSION;
  _n_levels = NoLevels;
  _N_CHILD = ((DIMENSION > 1) ? 4 * (DIMENSION - 1) : 2);
  _dcl_nel = ((DIMENSION > 1) ? 1 + NDOF_FEM + 3 + 1 + _N_CHILD : 1 + NDOF_FEM + 3 + 1 + _N_CHILD);
  _N_NDV = ((DIMENSION > 1) ? 4 : 2);
  _dclb_nel = NDOF_FEMB + 5;
  _n_subdomains = libMesh::global_n_processors();

  _ElPerProcPerLevel = new int*[_n_subdomains];
  for (int i = 0; i < _n_subdomains; i++) {
    _ElPerProcPerLevel[i] = new int[_n_levels];
    for (int j = 0; j < _n_levels; j++) _ElPerProcPerLevel[i][j] = 0;
  }

  return;
}

MGGenCase::~MGGenCase() {
  for (int i = 0; i < _n_subdomains; i++) delete[] _ElPerProcPerLevel[i];
  delete[] _ElPerProcPerLevel;
}

// =======================================================
/// This function is the main mesh generation function
void MGGenCase::GenCase() {
  // =======================================================

  int dim = DIMENSION;
  std::cout << " Dimension= " << dim << "  " << _dim << "\n";
#ifdef PRINT_TIME  // ------------------------------------------------
  std::cout << " \n ================================================ \n";
  std::clock_t start_timeA = std::clock();
#endif  // -----------------------------------------------------------

  // *************************************************************
  //                   Coarse mesh generator
  // *************************************************************
  // coarse mesh (msh0) generation case 0) from file 1) Internal
  Mesh* msh0 = new Mesh(_comm, _dim);
  // fine mesh (msht) generation case   0) from file 1) Internal
  Mesh* msht = new Mesh(_comm, dim);
  const int libmesh_gen = (int)_mgutils._geometry["libmesh_gen"];  // gen param flag
  switch (libmesh_gen) {
      // mesh generator cases 0) File 1) Internal

    case 0: {
      // mesh from file at level 0 -----------------------
      std::cout << " Reading coarse Mesh File at level 0 \n";
      std::ostringstream mesh_infile;
      mesh_infile << _mgutils._mesh_dir << _mgutils.get_file("F_MESH_READ");
      std::cout << mesh_infile.str() << std::endl;
      msh0->allow_renumbering(false);
      msh0->read(mesh_infile.str().c_str());
      msht->allow_renumbering(false);
      msht->read(mesh_infile.str().c_str());
      // From HEX8 or HEX20 to HEX27 etc..
      if ((int)_mgutils._geometry["second_order"] == 1) {
        msh0->all_second_order(true);
        msht->all_second_order(true);
      }
      break;
    }  // ----------------------------------------------------------

    case 1: {
      // Internal mesh generator at level 0 ---------------
      std::cout << " Internal mesh generator at level 0 \n";
      libMesh::ElemType libmname;  // libmesh name elem type
      int nintervx = 0;
      int nintervy = 0;
      int nintervz = 0;  // rectangular dimension

      const double LXB = _mgutils._geometry["LXB"];
      const double LXE = _mgutils._geometry["LXE"];
      const double LYB = _mgutils._geometry["LYB"];
      const double LYE = _mgutils._geometry["LYE"];
      const double LZB = _mgutils._geometry["LZB"];
      const double LZE = _mgutils._geometry["LZE"];

      switch (dim) {
          // dimension cases 1,2,3 dim

        case 3: {
          // dimension 3 --------------------------------------
          nintervx = (int)_mgutils._geometry["nintervx"];
          nintervy = (int)_mgutils._geometry["nintervy"];
          nintervz = (int)_mgutils._geometry["nintervz"];
          if (_geomel.name[0] == "Hex_27") {
            libmname = HEX27;
          } else if (_geomel.name[0] == "Tet_10") {
            libmname = TET10;
          }
          // else if (_geomel.name[0] == "Wedge_18")  libmname = PRISM18;
          MeshTools::Generation::build_cube(
              *msh0, nintervx, nintervy, nintervz,  // number of intervals
              LXB, LXE, LYB, LYE, LZB, LZE,         // rectangular dimensions
              libmname);                            // element type
          MeshTools::Generation::build_cube(
              *msht, nintervx, nintervy, nintervz,  // number of intervals
              LXB, LXE, LYB, LYE, LZB, LZE,         // rectangular dimensions
              libmname);                            // element type
          break;                                    // end case 3 ----------------------------------------
        }
        case 2: {
          // dimension 2 --------------------------------------
          nintervx = (int)_mgutils._geometry["nintervx"];
          nintervy = (int)_mgutils._geometry["nintervy"];
          if (_geomel.name[0] == "Quad_9") {
            libmname = QUAD9;
          } else if (_geomel.name[0] == "Tri_6") {
            libmname = TRI6;
          }
          MeshTools::Generation::build_square(
              *msh0, nintervx, nintervy, LXB, LXE, LYB, LYE, libmname);  // AX0,AX1,AY0,AY1
          MeshTools::Generation::build_square(
              *msht, nintervx, nintervy, LXB, LXE, LYB, LYE, libmname);  // AX0,AX1,AY0,AY1
          break;  //  end case 2 ---------------------------------------
        }
        default: {
          // dimension 1 -------------------------------------
          nintervx = (int)_mgutils._geometry["nintervx"];
          libmname = EDGE3;
          MeshTools::Generation::build_square(*msh0, nintervx, 0, LXB, LXE, 0., 0., libmname);
          MeshTools::Generation::build_square(*msht, nintervx, 0, LXB, LXE, 0., 0., libmname);
        }     // end default ----------------------------------------------
      }       // switch dim
      break;  // ------------------------------------------------------
    }
    default: std::cout << "MGGenCase: Error mesh coarse generation" << std::endl; abort();

  }  //--------------------------------------------------------------------------
  std::cout << "\n ========================================================== \n"
            << "  ================= Coarse Mesh information ============= \n";
  msh0->print_info();  // print the mesh at coarse level

#ifdef PRINT_TIME  // ------------------------------------------------
  std::clock_t end_timeA = std::clock();
  std::cout << " *+* Generation/Reading coarse mesh time ="
            << double(end_timeA - start_timeA) / CLOCKS_PER_SEC << std::endl;
  std::cout << " \n ================================================= \n ";
  std::clock_t start_timeB = std::clock();
#endif  // -------------------------------------------------------------

  // ********************************************************
  /// Refine the coarse mesh
  // ********************************************************
  const int NoLevels = (int)_mgutils._geometry["nolevels"];  // param
  // type mesh map -----------------------------------
  int* ttype_FEM;
  ttype_FEM = new int[2];
  ttype_FEM[0] = _geomel.n_q[0];
  ttype_FEM[1] = _geomel.n_q[1];

  // mesh_map_in ----------------------------------
  int* mesh_map_in;
  mesh_map_in = new int[2 * NoLevels];
  for (int itp = 0; itp < 2; itp++) {
    for (int ilev = 0; ilev < NoLevels; ilev++) { mesh_map_in[ilev + itp * NoLevels] = ttype_FEM[itp]; }
  }

  // LibMesh refinement msht[0]-> msht[NoLevels-1]-----------
  if (NoLevels > 1) {
    std::cout << "\n LibMesh Mesh Refinement ---------  \n";
    MeshRefinement mesh_refinement(*msht);
    mesh_refinement.uniformly_refine(NoLevels - 1);
  }

  // *********************************************************
  //  Generating Boundary Mesh from top level
  // *********************************************************
  std::cout << " LibMesh BOUNDARY generation --------- \n";
  // boundary of the coarse  mesh
  BoundaryMesh* bd_msh0 = new BoundaryMesh(_comm, msh0->mesh_dimension() - 1);
  msh0->boundary_info->sync(*bd_msh0);
  BoundaryMesh* bd_msht = new BoundaryMesh(_comm, msht->mesh_dimension() - 1);
  msht->boundary_info->sync(*bd_msht);

#ifdef PRINT_TIME  // ------------------------------------------
  std::clock_t end_timeB = std::clock();
  std::cout << " *+* Generation refined mesh time =" << double(end_timeB - start_timeB) / CLOCKS_PER_SEC
            << std::endl;
  std::cout << " \n ==================================================  \n ";
  std::clock_t start_timeC = std::clock();
#endif  // -------------------------------------------------------
  std::cout << "\n =============== Fine Mesh info ==================== \n";
  bd_msht->print_info();
  std::cout << "\n =============== Fine Mesh info ==================== \n";
  msht->print_info();

  // ********************************************************
  ///              PrintMesh
  // ********************************************************
  if (libMesh::global_processor_id() == 0) { printMesh(*bd_msht, *msht, *bd_msh0, *msh0, mesh_map_in); }

#ifdef PRINT_TIME  // ---------------------------------------------------------
  std::clock_t end_timeC = std::clock();
  std::cout << " Print and Operators time =" << double(end_timeC - start_timeC) / CLOCKS_PER_SEC << std::endl;
#endif  // -----------------------------------------------------------------------

  // ********************************************************
  /// clean and stop
  // ********************************************************
  delete[] ttype_FEM;
  delete[] mesh_map_in;
  delete bd_msht;
  delete msht;
  delete bd_msh0;
  delete msh0;

#ifdef PRINT_TIME  // --------------------------------------------------------------------------------------------
  std::cout << " \n ====================================   \n";
  std::cout << " +*+*+* Total time =" << double(end_timeC - start_timeA) / CLOCKS_PER_SEC << std::endl;
#endif  // --------------------------------------------------------------------------------------------

  std::cout << "\n +++++++++ End GenCase ++++++++++++++++ " << std::endl;
  return;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
