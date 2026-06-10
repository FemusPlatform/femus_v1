#pragma once



// using namespace std --------------->
#include <map>
// class MGGeomElBase ----------------->

#include "MGGeom_L1_El.h"
// #include "MGFE_L2_template.h"



// ============================================================================
// ============================ MGGeomMap =======================================
// ============================================================================
class MGGeomMap {
  // ===============  data ====================================================
 private:
  // map<int, MGFE1*> Geom_map;  ///< Map of the FEMs
// map<int, MGFE<2,LAGRANGE>*> Geom_map;  ///< Map of the FEMs
std::map<int, MGGeomElBase *> Geom_map;  ///< Map of the FEMs
 public:
  // =======================  Constructor ===================================
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  /// Constructor
  MGGeomMap(){};
  /// Destructor
//   ~MGGeomMap();
  ///@}

  // ====================  functions ==========================================
  ///< This function inserts in the map a FEM
//  void set_FE(MGFE* mgfe) {  Geom_map.insert(make_pair(mgfe->_order, mgfe));} 

// ================================================================================================
/// This function destroys the MGFE class
~MGGeomMap() {  // ==============================
  // for (std::map<int, MGFE1*>::iterator itr = Geom_map.begin(); itr != Geom_map.end(); itr++) {
    // (itr->second)->~MGFE1();
  // }
   // for (std::map<int, MGFE<2,LAGRANGE>*>::iterator itr = Geom_map.begin(); itr != Geom_map.end(); itr++) {
    // (itr->second)->~MGFE();
  // }
   for (std::map<int, MGGeomElBase*>::iterator itr = Geom_map.begin(); itr != Geom_map.end(); itr++) {
    (itr->second)->~MGGeomElBase();
  }
  Geom_map.clear();
//  P_INFO (  std::cout << " ~MGGeomMap() called \n";)
}
  
  
  
  ///@{ \name FE GET/SET
//   void set_FE(MGFE* mgfe);  //{Geom_map.insert(make_pair(mgfe->_order,mgfe));} ///< Insert in the map a FEM
// void set_FE(MGFE1* mgfe) {  Geom_map.insert(make_pair(mgfe->_order+mgfe->_FamType*3, mgfe));} 
// void set_FE(MGFE<2,LAGRANGE>* mgfe,int order) {  Geom_map.insert(make_pair(order, mgfe));}
void set_Geom(MGGeomElBase* mgfe,int order) {  Geom_map.insert(std::make_pair(order, mgfe));}
  /// Get a FEM from the map
  // inline MGFE<2,LAGRANGE>* get_FE(const int label) const { return Geom_map.find(label)->second; }
    inline MGGeomElBase* get_Geom(const int label) const { return Geom_map.find(label)->second; }
  ///@}

  ///@{ \name ITERATORS FOR FE MAP
  // typedef std::map<int, MGFE<2,LAGRANGE>*>::const_iterator const_giterator;  ///< Iterator for the map
  // const_giterator gc_begin() const { return Geom_map.begin(); }    ///< Return begin of the map
 typedef std::map<int, MGGeomElBase*>::const_iterator const_giterator;  ///< Iterator for the map
  const_giterator gc_begin() const { return Geom_map.begin(); }    ///< Return begin of the map
  ///< Return end of the map
  const_giterator gc_end() const { return Geom_map.end(); }
  ///@}
};
