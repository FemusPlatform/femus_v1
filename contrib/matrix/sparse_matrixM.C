// ========================================
//            SparseMatrix class
// ========================================
// std lib includes -------------------------
#include <fstream>  // std::ofsream
#include <sstream>  // std::ostringstream
#include "hdf5.h"   // hdf5

// this class ----------------------------
#include "sparse_matrixM.h"  // definition header

// configuration files -------------------
#include "Solverlib_conf.h"  // in ./config

// local includes ------------------------
#include "laspack_matrixM.h"  // laspack matrix (son of this class)
#include "numeric_vectorM.h"  // numerical base vector class
#include "petsc_matrixM.h"    // petsc matrix   (son of this class)
                              // #include "trilinos_epetra_matrix.h"  // trilinos matrix

// ====================================================================
//            SparseMatrix Methods: Constructor/Destructor/init/build
// ===================================================================

// =====================================================================================
/// Constructor;  before usage call init(...).
SparseMatrixM::SparseMatrixM(const ParallelM::Communicator& comm  // parallel communicator <-
                             )
    : ParallelObjectM(comm), _is_initialized(false) {}

// =====================================================================================
/// This function builds a  SparseMatrixM using the linear solver
/// package specified by  solver_package
std::unique_ptr<SparseMatrixM> SparseMatrixM::build(  // -----
    const ParallelM::Communicator& comm,
    const SolverPackageM solver_package  //  solver_package
) {  // =================================================================================
  // Build the appropriate vector
  switch (solver_package) {
#ifdef HAVE_LASPACKM  // ----------------------------
    case LASPACK_SOLVERSM: {
      std::unique_ptr<SparseMatrixM> ap(new LaspackMatrixM(comm));
      return ap;
    }
#endif
#ifdef HAVE_PETSCM  // ------------------------------
    case PETSC_SOLVERSM: {
      std::unique_ptr<SparseMatrixM> ap(new PetscMatrixM(comm));
      return ap;
    }
#endif
#ifdef HAVE_TRILINOSM  // ----------------------------
    case TRILINOS_SOLVERSM: {
      std::unique_ptr<SparseMatrixM> ap(new EpetraMatrix<double>);
      return ap;
    }
#endif
    default: std::cerr << "SolverPackageM solver_package:  Unrecognized: " << solver_package; abort();
  }
  std::cerr << "SolverPackageM solver_package:  Unrecognized: " << solver_package;
  //   std::unique_ptr<SparseMatrixM > ap(NULL);
  //   return ap;
}

// =================================================
//            SparseMatrix Methods: Add/mult
// =================================================

// =========================================================
void SparseMatrixM::vector_mult(NumericVectorM& dest, const NumericVectorM& arg) const {
  dest.zero();
  this->vector_mult_add(dest, arg);
}

// ===========================================================
void SparseMatrixM::vector_mult_add(NumericVectorM& dest, const NumericVectorM& arg) const {
  /* This functionality is actually implemented in the \p NumericVector class.  */
  dest.add_vector(arg, *this);
}

// =================================================
//            SparseMatrix Methods: setting/return
// =================================================

// ==========================================================
/// This function sets all row entries to 0 then puts diag_value in the diagonal entry
void SparseMatrixM::zero_rows(std::vector<int>&, double) {
  std::cerr << " SparseMatrixM::zero_rows: Not implememnted ";
  abort();
}

// ===============================================
//            SparseMatrix Methods: Print/Read
// ==============================================

// =====================================================
void SparseMatrixM::print(const std::string& name) const {
  std::ofstream out(name.c_str());
  print(out);
}



// =====================================================
/// This function reads len sparse matrix structures
void SparseMatrixM::read_len_hdf5(
    const std::string namefile,  // file name
    const int i,
    const  int j,
//     const std::string mode,      // linear or quadratic
    int len_row[],               // row lengths
    int len_off_row[]            // row off entry lengths
) {
  hid_t file_id = H5Fopen(namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  // row lengths -------------------------------------------------------------
//   std::ostringstream name_dst;
//   if (mode == "0" || mode == "1" || mode == "2")   name_dst << "LEN0" << mode;
//   else  name_dst << "LEN" << mode;

  hid_t dataset = H5Dopen(file_id, ("LEN"+std::to_string(i)+std::to_string(j)).c_str(), H5P_DEFAULT);
  hid_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, len_row);
  assert(status == 0);  H5Dclose(dataset);
  // matrix row off entry lengths --------------------------------------------
//   name_dst.str("");
//   if (mode == "0" || mode == "1" || mode == "2")   name_dst << "OFFLEN0" << mode;
//   else  name_dst << "OFFLEN" << mode;

  dataset = H5Dopen(file_id,("OFFLEN"+std::to_string(i)+std::to_string(j)).c_str(), H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, len_off_row);
  // closing 
  assert(status == 0);  H5Dclose(dataset);  H5Fclose(file_id);
  return;
}

// =====================================================
/// This function reads pos sparse matrix structure
void SparseMatrixM::read_pos_hdf5(
    const std::string namefile,  // file name
    const int i,
    const int j,
//     const int mode,              // linear or quadratic
    int pos_row[]                // compressed row positions
) {
  hid_t file_id = H5Fopen(namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

//   std::ostringstream name_dst;  name_dst.str("");  name_dst << "POS" << mode;
  hid_t dataset = H5Dopen(file_id, ("POS"+std::to_string(i)+std::to_string(j)).c_str(),H5P_DEFAULT);
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_row);
  // closing
  H5Dclose(dataset);  H5Fclose(file_id);
  return;
}
// =====================================================
/// This function reads quad-linear matrix dimensions
int SparseMatrixM::read_dim_hdf5(
    const std::string namefile,  // file name
    const int i,
    const int j,
//     const std::string mode,      // linear or quadratic
    int ldim[]                   // dimensions
) {
  hid_t file_id = H5Fopen(namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
//   std::ostringstream name_label;
//   if (mode == "0" || mode == "1" || mode == "2")   name_label << "DIM0" << mode;
//   else    name_label << "DIM" << mode;

  hid_t dataset = H5Dopen( file_id, ("DIM"+std::to_string(i)+std::to_string(j)).c_str(), H5P_DEFAULT);
  hid_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ldim);
  // closing 
  assert(status == 0);  H5Dclose(dataset);  H5Fclose(file_id);
  return ((int)status);
}
// =====================================================
/// This function reads quad-linear matrix dimensions
int SparseMatrixM::read_hdf5(
    const std::string namefile,  // file name
    const std::string mode,      // linear or quadratic
    int data[]                   // dimensions
) {
  hid_t file_id = H5Fopen(namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
//   std::ostringstream name_label;
//   if (mode == "0" || mode == "1" || mode == "2")   name_label << "DIM0" << mode;
//   else    name_label << "DIM" << mode;

  hid_t dataset = H5Dopen( file_id,mode.c_str(), H5P_DEFAULT);
  hid_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  // closing 
  assert(status == 0);  
  H5Dclose(dataset);  
  H5Fclose(file_id);
  return ((int)status);
}
//------------------------------------------------------------------
