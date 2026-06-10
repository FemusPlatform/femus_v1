#pragma once

#include <string> 
#include <iostream> 
#include <unistd.h>
// using namespace std --------------->
#include <map>
#include <memory>
// class MGFEbase ----------------->

#include "MGFE_L1_base.h"
// #include "MGFE_L2_template.h"
#include "MGquadrature_L1_gauss.h"
#include "MGGeom_L0_enum.h"


// ============================================================================
// ============================ MGFEMap =======================================
// ============================================================================
class MGFEMap {
  // ===============  data ====================================================
 private:
  // map<int, MGFE1*> FE_map;  ///< Map of the FEMs
// map<int, MGFE<2,LAGRANGE>*> FE_map;  ///< Map of the FEMs
  int dim_;
  GEOFamily geo_fam_;  // FEFamily::LAGRANGE
  ElemType maxElemType_;//   ElemType::HEX27 
  Order  Q_order_;      // Order::FIFTH quadrature
  

  std::string filename_femap_;

  
  std::unique_ptr<QGauss> qrule_;
  std::unique_ptr<QGauss> face_qrule_;
  
  std::unique_ptr<FEType> fetq_;
  std::unique_ptr<FEType> fetl_;
  std::unique_ptr<FEType> fetk_;
  
  std::map<int, std::unique_ptr<MGFEbase>> FE_map;  ///< Map of the FEMs

 public:
  
  // =======================  Constructor ===================================
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  /// Constructor
  // char* buff;//automatically cleaned when it exits scope
  //    return std::string(getcwd(buff,255));
 MGFEMap(int dim_in,GEOFamily geo_fam_in=LAGRANGE27, Order  qorder_in=FIFTH);

   // MGFEMap(std::string filename_femap_in){filename_femap_=filename_femap_in;}
  /// Destructor
  //   ~MGFEMap();
  ///@}
  // void init(Order  qorder);
  // ====================  functions ==========================================
  ///< This function inserts in the map a FEM
//  void set_FE(MGFE* mgfe) {  FE_map.insert(make_pair(mgfe->_order, mgfe));} 

// ================================================================================================
/// This function destroys the MGFE class
~MGFEMap() {  // ==============================
  // for (std::map<int, MGFE1*>::iterator itr = FE_map.begin(); itr != FE_map.end(); itr++) {
    // (itr->second)->~MGFE1();
  // }
   // for (std::map<int, MGFE<2,LAGRANGE>*>::iterator itr = FE_map.begin(); itr != FE_map.end(); itr++) {
    // (itr->second)->~MGFE();
  // }
  FE_map.clear();
  
//  P_INFO (  std::cout << " ~MGFEMap() called \n";)
}
  
  
// void set_FE(std::unique_ptr<MGFEbase> & mgfe,int order) {  FE_map.insert(std::make_pair(order, mgfe.release()));}  
  ///@{ \name FE GET/SET
  
  
 void attach_qrule(std::unique_ptr<QGauss> q){ qrule_.reset(q.release());     }    //{FE_map.insert(make_pair(mgfe._order,mgfe));} ///< Insert in the map a FEM
// void set_FE(MGFE1* mgfe) {  FE_map.insert(make_pair(mgfe._order+mgfe._FamType*3, mgfe));} 
// void set_FE(MGFE<2,LAGRANGE>* mgfe,int order) {  FE_map.insert(make_pair(order, mgfe));}
void set_FE(std::unique_ptr<MGFEbase> & mgfe,int order) {  FE_map.insert(std::make_pair(order, mgfe.release()));}
  /// Get a FEM from the map
  // inline MGFE<2,LAGRANGE>* get_FE(const int label) const { return FE_map.find(label).second; }
    inline MGFEbase* get_FE(const int label) const { return FE_map.find(label)->second.get(); }
  ///@}

  ///@{ \name ITERATORS FOR FE MAP
  // typedef std::map<int, MGFE<2,LAGRANGE>*>::const_iterator const_giterator;  ///< Iterator for the map
  // const_giterator gc_begin() const { return FE_map.begin(); }    ///< Return begin of the map
 // typedef std::map<int, MGFEbase*>::const_iterator const_giterator;  ///< Iterator for the map
  // const_giterator gc_begin() const { return FE_map.begin(); }    ///< Return begin of the map
  // ///< Return end of the map
  // const_giterator gc_end() const { return FE_map.end(); }
  std::string get_filename_femap(){return filename_femap_; }
  ///@}
};
