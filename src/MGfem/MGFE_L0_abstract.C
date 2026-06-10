// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// libmesh includes

// #include "libmesh/libmesh_logging.h"
// #include "libmesh/enum_elem_type.h"
// #include "libmesh/boundary_info.h"
// #include "libmesh/mesh_base.h"
// #include "libmesh/dense_matrix.h"
// #include "libmesh/dense_vector.h"
// #include "libmesh/dof_map.h"
// #include "libmesh/elem.h"
// #include "libmesh/fe_interface.h"
// #include "libmesh/numeric_vector.h"
// #include "libmesh/periodic_boundaries.h"
// #include "libmesh/periodic_boundary.h"
// #include "libmesh/quadrature.h"
// #include "libmesh/quadrature_gauss.h"
// #include "libmesh/remote_elem.h"
// #include "libmesh/tensor_value.h"
// #include "libmesh/threads.h"
// #include "libmesh/enum_elem_type.h"
// #include "libmesh/enum_to_string.h"
// 
// #ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
// #include "libmesh/inf_fe.h"
// #include "libmesh/fe_interface_macros.h"
// #endif
#include "MGFE_L0_abstract.h"

#include<cassert>


MGFE0::MGFE0(const unsigned int d, const FEType & fet) :
  // _fe_map(nullptr ),
  _dim(d),
  _calculations_started(false), _calculate_nothing(false),  _calculate_map(false),
  _calculate_phi(false),  _calculate_dphi(false),  _calculate_d2phi(false),
  _calculate_curl_phi(false),  _calculate_div_phi(false),  _calculate_dphiref(false),
  _fe_type(fet),
  _elem_type(INVALID_ELEM),
  _elem_p_level(0), _p_level(0),
  _qrule(0),  // Quadrature -> int
  _shapes_on_quadrature(false),
  _n_total_qp(0),
  _add_p_level_in_reinit(false)
{
}

// =============================================================================
MGFE0::~MGFE0() = default;


std::unique_ptr<MGFE0> MGFE0::build(
  const unsigned int dim, const FEType & fet
){
  switch (dim) {
      // 0D
    case 0:
      {
        switch (fet.family)   {
          // case CLOUGH:
          //   return std::make_unique<FE<0,CLOUGH>>(fet);
          // case HERMITE:
          //   return std::make_unique<FE<0,HERMITE>>(fet);
          // case LAGRANGE:
          //   return std::make_unique<FE<0,LAGRANGE>>(fet);
          // case LAGRANGE_VEC:
          //   return std::make_unique<FE<0,LAGRANGE_VEC>>(fet);
          // case L2_LAGRANGE:
          //   return std::make_unique<FE<0,L2_LAGRANGE>>(fet);
          // case L2_LAGRANGE_VEC:
          //   return std::make_unique<FE<0,L2_LAGRANGE_VEC>>(fet);
          // case HIERARCHIC_VEC:
          //   return std::make_unique<FE<0,HIERARCHIC_VEC>>(fet);
          // case HIERARCHIC:
          //   return std::make_unique<FE<0,HIERARCHIC>>(fet);
          // case L2_HIERARCHIC:
          //   return std::make_unique<FE<0,L2_HIERARCHIC>>(fet);
          // case L2_HIERARCHIC_VEC:
          //   return std::make_unique<FE<0,L2_HIERARCHIC_VEC>>(fet);
          // case SIDE_HIERARCHIC:
          //   return std::make_unique<FE<0,SIDE_HIERARCHIC>>(fet);
          // case MONOMIAL:
          //   return std::make_unique<FE<0,MONOMIAL>>(fet);
          // case MONOMIAL_VEC:
          //   return std::make_unique<FE<0,MONOMIAL_VEC>>(fet);
          // case XYZ:
          //   return std::make_unique<FEXYZ<0>>(fet);
          // case SCALAR:
          //   return std::make_unique<FEScalar<0>>(fet);
          default:
            std::cerr<<"ERROR: Bad FEType.family= " << fet.family;
          }
      }
      // 1D
    case 1:
      {
        switch (fet.family)          {
          // case CLOUGH:
          //   return std::make_unique<FE<1,CLOUGH>>(fet);
          // case HERMITE:
          //   return std::make_unique<FE<1,HERMITE>>(fet);
          // case LAGRANGE:
          //   return std::make_unique<FE<1,LAGRANGE>>(fet);
          // case LAGRANGE_VEC:
          //   return std::make_unique<FE<1,LAGRANGE_VEC>>(fet);
          // case L2_LAGRANGE:
          //   return std::make_unique<FE<1,L2_LAGRANGE>>(fet);
          // case L2_LAGRANGE_VEC:
          //   return std::make_unique<FE<1,L2_LAGRANGE_VEC>>(fet);
          // case HIERARCHIC_VEC:
          //   return std::make_unique<FE<1,HIERARCHIC_VEC>>(fet);
          // case HIERARCHIC:
          //   return std::make_unique<FE<1,HIERARCHIC>>(fet);
          // case L2_HIERARCHIC:
          //   return std::make_unique<FE<1,L2_HIERARCHIC>>(fet);
          // case L2_HIERARCHIC_VEC:
          //   return std::make_unique<FE<1,L2_HIERARCHIC_VEC>>(fet);
          // case SIDE_HIERARCHIC:
          //   return std::make_unique<FE<1,SIDE_HIERARCHIC>>(fet);
          // case MONOMIAL:
          //   return std::make_unique<FE<1,MONOMIAL>>(fet);
          // case MONOMIAL_VEC:
          //   return std::make_unique<FE<1,MONOMIAL_VEC>>(fet);
          // case XYZ:
          //   return std::make_unique<FEXYZ<1>>(fet);
          // case SCALAR:
          //   return std::make_unique<FEScalar<1>>(fet);
          default:
            std::cerr<<"ERROR: Bad FEType.family= " << fet.family;
          }
      }
      // 2D ============================================================================================
    case 2:
      {
        switch (fet.family)
          {
          // case CLOUGH:
          //   return std::make_unique<FE<2,CLOUGH>>(fet);
          // case HERMITE:
          //   return std::make_unique<FE<2,HERMITE>>(fet);
          // case LAGRANGE:
          //   return std::make_unique<FE<2,LAGRANGE>>(fet);
          // case LAGRANGE_VEC:
          //   return std::make_unique<FE<2,LAGRANGE_VEC>>(fet);
          // case L2_LAGRANGE:
          //   return std::make_unique<FE<2,L2_LAGRANGE>>(fet);
          // case L2_LAGRANGE_VEC:
          //   return std::make_unique<FE<2,L2_LAGRANGE_VEC>>(fet);
          // case HIERARCHIC_VEC:
          //   return std::make_unique<FE<2,HIERARCHIC_VEC>>(fet);
          // case HIERARCHIC:
          //   return std::make_unique<FE<2,HIERARCHIC>>(fet);
          // case L2_HIERARCHIC:
          //   return std::make_unique<FE<2,L2_HIERARCHIC>>(fet);
          // case L2_HIERARCHIC_VEC:
          //   return std::make_unique<FE<2,L2_HIERARCHIC_VEC>>(fet);
          // case SIDE_HIERARCHIC:
          //   return std::make_unique<FE<2,SIDE_HIERARCHIC>>(fet);
          // case MONOMIAL:
          //   return std::make_unique<FE<2,MONOMIAL>>(fet);
          // case MONOMIAL_VEC:
          //   return std::make_unique<FE<2,MONOMIAL_VEC>>(fet);
          // case XYZ:
          //   return std::make_unique<FEXYZ<2>>(fet);
          // case SCALAR:
          //   return std::make_unique<FEScalar<2>>(fet);
          // case NEDELEC_ONE:
          //   return std::make_unique<FENedelecOne<2>>(fet);
          // case RAVIART_THOMAS:
          //   return std::make_unique<FERaviartThomas<2>>(fet);
          // case L2_RAVIART_THOMAS:
          //   return std::make_unique<FEL2RaviartThomas<2>>(fet);
          // case SUBDIVISION:
          //   return std::make_unique<FESubdivision>(fet);
          default:
            std::cerr<<"ERROR: Bad FEType.family= " << fet.family;
          }
      }
    case 3: // 3D-------------------------------------------------------------------------------
      {
        switch (fet.family)
          {
          // case CLOUGH:
          //   std::cerr<<"ERROR: Clough-Tocher elements currently only support 1D and 2D";
          // case HERMITE:
          //   return std::make_unique<FE<3,HERMITE>>(fet);
          // case LAGRANGE:
          //   return std::make_unique<FE<3,LAGRANGE>>(fet);
          // case LAGRANGE_VEC:
          //   return std::make_unique<FE<3,LAGRANGE_VEC>>(fet);
          // case L2_LAGRANGE:
          //   return std::make_unique<FE<3,L2_LAGRANGE>>(fet);
          // case L2_LAGRANGE_VEC:
          //   return std::make_unique<FE<3,L2_LAGRANGE_VEC>>(fet);
          // case HIERARCHIC_VEC:
          //   return std::make_unique<FE<3,HIERARCHIC_VEC>>(fet);
          // case HIERARCHIC:
          //   return std::make_unique<FE<3,HIERARCHIC>>(fet);
          // case L2_HIERARCHIC:
          //   return std::make_unique<FE<3,L2_HIERARCHIC>>(fet);
          // case L2_HIERARCHIC_VEC:
          //   return std::make_unique<FE<3,L2_HIERARCHIC_VEC>>(fet);
          // case SIDE_HIERARCHIC:
          //   return std::make_unique<FE<3,SIDE_HIERARCHIC>>(fet);
          // case MONOMIAL:
          //   return std::make_unique<FE<3,MONOMIAL>>(fet);
          // case MONOMIAL_VEC:
          //   return std::make_unique<FE<3,MONOMIAL_VEC>>(fet);
          // case XYZ:
          //   return std::make_unique<FEXYZ<3>>(fet);
          // case SCALAR:
          //   return std::make_unique<FEScalar<3>>(fet);
          // case NEDELEC_ONE:
          //   return std::make_unique<FENedelecOne<3>>(fet);
          // case RAVIART_THOMAS:
          //   return std::make_unique<FERaviartThomas<3>>(fet);
          // case L2_RAVIART_THOMAS:
          //   return std::make_unique<FEL2RaviartThomas<3>>(fet);
          default:
            std::cerr<<"ERROR: Bad FEType.family= " << fet.family;
          }
      }
    default:
      std::cerr<<"Invalid dimension dim = " << dim;
    }
}

//==============================================================================
int MGFE0::get_refspace_n_nodes(const ElemType itemType){
   // nodes.resize(Elem::type_to_n_nodes_map[itemType]);
  switch(itemType)  {
    case NODEELEM: return 1;    
    case EDGE3:    return 3;  
    case EDGE2:    return 2;
    case EDGE4:    return 4;
    case TRI7: return 7;
    case TRI6:  return 6;   
    case TRI3:  return 3;
    case TRISHELL3:   return 3;  
    case QUAD9:  return 9;
    case QUADSHELL9: return 9;
    case QUAD8:  return 8;
    case QUADSHELL8: return 8;
    case QUAD4: return 4;
    case QUADSHELL4:  return 4;
    case TET14: return 14;
    case TET10: return 10; 
    case TET4:  return 4;
    case HEX27: return 27;
    case HEX20: return 20;
    case HEX8: return 8;
    case PRISM21: return 21;
    case PRISM20:  return 20;
    case PRISM18:return 18;
    case PRISM15: return 15;
    case PRISM6:return 6;
    case PYRAMID18:return 18;
    case PYRAMID14: return 14;
    case PYRAMID13:return 13;
    case PYRAMID5: return 5;
    default:  std::cerr<<"ERROR: Unknown element type " << itemType; return -1;
    }
}
ElemType MGFE0::get_ref_elemtype(const ElemType elemType){
  
   switch(elemType)    {
    case NODEELEM: {int dim=get_dim(); switch(dim){
                          case 0: return NODEELEM;
                          case 1: return EDGE3;
                          case 2:  return QUAD9;
                          case 3:  return HEX27;
                          default:  abort();return NODEELEM;
                          }
                    }
    
    case EDGE3:   case EDGE2:  case EDGE4: return EDGE3; 
    case TRI7:  case TRI6:     case TRI3:   case TRISHELL3: return TRI6;
    case QUAD9:  case QUADSHELL9:   case QUAD8:    case QUADSHELL8: case QUAD4: case QUADSHELL4: return QUAD9;
    case TET14:  case TET10:  case TET4: return TET10;
    case HEX27:  case HEX20:   case HEX8: return HEX27;
      
    case PRISM21:  case PRISM20:   case PRISM18:    case PRISM15:  case PRISM6:  return PRISM15;
    case PYRAMID18: case PYRAMID14:   case PYRAMID13:  case PYRAMID5:  return PYRAMID14;
    default:   std::cerr<<"ERROR: Unknown element type " << elemType;
    }
    return NODEELEM;
}
// ====================================================
  
 ElemType MGFE0::get_ref_face_elemtype(const ElemType elemType){
  
   switch(elemType)    {
    case NODEELEM: {int dim=get_dim(); switch(dim){
      case 0: return NODEELEM;
                          case 1: return NODEELEM;
                          case 2: return EDGE3;
                          case 3:  return QUAD9;
                          default:  abort();return NODEELEM;
                          }
                    }
    
    case EDGE3:   case EDGE2:  case EDGE4: return NODEELEM; 
    case TRI7:  case TRI6:    return EDGE3;
    case TRI3:   case TRISHELL3: return EDGE2;
    case QUAD9:  case QUADSHELL9:  case QUADSHELL8: case QUAD8:    return EDGE3;
        case QUAD4: case QUADSHELL4: return  EDGE2;
    case TET14:  case TET10:  return TRI6; 
    case TET4: return TRI3;
    case HEX27:  case HEX20: return QUAD9;  
    case HEX8: return QUAD4;
      
    case PRISM21:  case PRISM20:   case PRISM18:    case PRISM15:  return QUAD9;
    case PRISM6:  return QUAD4;
    case PYRAMID18: case PYRAMID14:   case PYRAMID13: return QUAD9; 
    case PYRAMID5:  return QUAD4;
    default:   std::cerr<<"ERROR: Unknown element type " << elemType;
    }
    return NODEELEM;
} 
  
  
//==============================================================================
void MGFE0::get_refspace_nodes(const ElemType itemType,double nodes[]){ // _dim point
   // nodes.resize(Elem::type_to_n_nodes_map[itemType]);
  switch(itemType)    {
    case NODEELEM: {int dim=get_dim(); for(int i=0;i<dim;i++) nodes[i] = 0.;   return;    }
    
    case EDGE3:    { nodes[2] = 0.;}
    case EDGE2:    { nodes[1] = 1.;   //  nodes[1] = Point (1.,0.,0.);
                     nodes[0] =-1;;  //  nodes[0] = Point (-1.,0.,0.); 
        return;  }
    case EDGE4:  {   // not nested with EDGE3
        nodes[3] =1./3.;// nodes[3] - Point (1./3.,0.,0.);
        nodes[2] =-1/3.;// nodes[2] = Point (-1./3.,0.,0.);
        nodes[1] =1.;// nodes[1] = Point (1.,0.,0.);                        
        nodes[0] =-1.;//nodes[0] =  Point (-1.,0.,0.); 
        return;}
        
    case TRI7:  {       nodes[12] =1./3.;nodes[13] =1./3.;}// nodes[6] = Point (1./3.,1./3.,0.);   
    case TRI6:  {       nodes[10] =0.;   nodes[11] =.5;  }//nodes[5] = Point (0.,.5,0.);
                        nodes[8] =.5;    nodes[9]  =.5;//nodes[4] = Point (.5,.5,0.); 
                        nodes[6] =.5;    nodes[7]  =0.;// nodes[3] = Point (.5,0.,0.); 
    case TRI3:
    case TRISHELL3: {  nodes[0] =0.;   nodes[1] =0.;//nodes[0] = Point (0.,0.,0.); 
                       nodes[2] =1.;   nodes[3] =0.;//  nodes[1] = Point (1.,0.,0.); 
                       nodes[4] =0.;   nodes[5] =1.;//nodes[2] = Point (0.,1.,0.);  
          return;  }
    
    case QUAD9:
    case QUADSHELL9:  {  
       nodes[16] =0.;   nodes[17] =0.;}//nodes[8] = Point (0.,0.,0.); 
    case QUAD8:
    case QUADSHELL8: {  
       nodes[8] =0.;   nodes[9] =-1.;//nodes[4] = Point (0.,-1.,0.);
       nodes[10] =1.;   nodes[11] =0.;//  nodes[5] = Point (1.,0.,0.);
       nodes[12] =0.;   nodes[13] =1.;// nodes[6] = Point (0.,1.,0.);  
       nodes[14] =-1.;  nodes[15] =0.;}//nodes[7] = Point (-1.,0.,0.); 
    case QUAD4:
    case QUADSHELL4:  {   
       nodes[0] =-1.;  nodes[1] =-1.;// nodes[0] = Point (-1.,-1.,0.); 
       nodes[2] =1.;   nodes[3] =-1.;// nodes[1] = Point (1.,-1.,0.);
       nodes[4] =1.;   nodes[5] =1.;//  nodes[2] = Point (1.,1.,0.); 
       nodes[6] =-1.;  nodes[7] =1.;// nodes[3] = Point (-1.,1.,0.);
       return;  }
       
    case TET14: {
         nodes[30] =1./3.;nodes[31] =1./3.;nodes[32] =0.;//nodes[10] = Point (1/Real(3),1/Real(3),0.);  
         nodes[33] =1./3.;nodes[34] =0;    nodes[35] =1./3.;//nodes[11] = Point (1/Real(3),0.,1/Real(3));
         nodes[36] =1./3.;nodes[37] =1./3.;nodes[38] =1./3.;//nodes[12] = Point (1/Real(3),1/Real(3),1/Real(3));  
         nodes[39] =0.;   nodes[40] =1./3.;nodes[41] =1./3.;}//nodes[13] = Point (0.,1/Real(3),1/Real(3));
    case TET10:  {
         nodes[12] =.5;  nodes[13] =0.;nodes[14] =0.;//nodes[4] = Point (.5,0.,0.);
         nodes[15] =.5;  nodes[16] =.5;nodes[17] =0.;//nodes[5] = Point (.5,.5,0.);
         nodes[18] =0.;  nodes[19] =.5;nodes[20] =0.;//nodes[6] = Point (0.,.5,0.);
         nodes[21] =0.;  nodes[22] =0.;nodes[23] =0.5;//nodes[7] = Point (0.,0.,.5); 
         nodes[24] =.5;  nodes[25] =0.;nodes[26] =0.5;//nodes[8] = Point (.5,0.,.5);  
         nodes[27] =0.;  nodes[28] =.5;nodes[29] =0.5;}// nodes[9] = Point (0.,.5,.5);
    case TET4:      {
         nodes[0] =0.;  nodes[1] =0.; nodes[2] =0.;// nodes[0] = Point (0.,0.,0.); 
         nodes[3] =1.;  nodes[4] =0.; nodes[5] =0.;//nodes[1] = Point (1.,0.,0.);
         nodes[6] =0.;  nodes[7] =1.; nodes[8] =0.;//nodes[2] = Point (0.,1.,0.); 
         nodes[9] =0.;  nodes[10]=0.; nodes[11]=1.;//nodes[3] = Point (0.,0.,1.);
        return;
      }
      
    case HEX27: {
        nodes[60] =0.;  nodes[61]=0.;  nodes[62] =-1.;//nodes[20] = Point (0.,0.,-1.); 
        nodes[63] =0.;  nodes[64]=-1.; nodes[65] =0.;//nodes[21] = Point (0.,-1.,0.);    
        nodes[66] =1.;  nodes[67]=0.;  nodes[68] =0.;//nodes[22] = Point (1.,0.,0.);
        nodes[69] =0.;  nodes[70]=1.;  nodes[71] =0.;//nodes[23] = Point (0.,1.,0.);    
        nodes[72] =-1.; nodes[73]=0.;  nodes[74] =0.;//nodes[24] = Point (-1.,0.,0.);   
        nodes[75] =0.;  nodes[76]=0.;  nodes[77] =1.;//nodes[25] = Point (0.,0.,1.);   
        nodes[78] =0.;  nodes[79]=0.;  nodes[80] =0.;}//nodes[26] = Point (0.,0.,0.);
    case HEX20:   {
        nodes[24] =0.;  nodes[25] =-1.;nodes[26] =-1.;//nodes[8] = Point (0.,-1.,-1.);  
        nodes[27] =1.;  nodes[28] =0.; nodes[29] =-1.;//nodes[9] = Point (1.,0.,-1.); 
        nodes[30] =0.;  nodes[31] =1.; nodes[32] =-1.;//nodes[10] = Point (0.,1.,-1.);
        nodes[33] =-1.; nodes[34] =0.; nodes[35] =-1.;//nodes[11] = Point (-1.,0.,-1.); 
        nodes[36] =-1.; nodes[37] =-1.;nodes[38] =0.;//nodes[12] = Point (-1.,-1.,0.);   
        nodes[39] =1.;  nodes[40] =-1.;nodes[41] =0.;//nodes[13] = Point (1.,-1.,0.);
        nodes[42] =1.;  nodes[43] =1.; nodes[44] =0.;//nodes[14] = Point (1.,1.,0.);    
        nodes[45] =-1.; nodes[46] =1.; nodes[47] =0.;//nodes[15] = Point (-1.,1.,0.);   
        nodes[48] =0.;  nodes[49] =-1.;nodes[50] =1.;//nodes[16] = Point (0.,-1.,1.);
        nodes[51] =1.;  nodes[52] =0.; nodes[53] =1.;// nodes[17] = Point (1.,0.,1.);  
        nodes[54] =0.;  nodes[55] =1.; nodes[56] =1.;//nodes[18] = Point (0.,1.,1.);  
        nodes[57] =-1.; nodes[58] =0.; nodes[59] =1.;}//nodes[19] = Point (-1.,0.,1.);
    case HEX8:      {
      nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//  nodes[0] = Point (-1.,-1.,-1.);  
      nodes[3] =1.;  nodes[4] =-1.;nodes[5] =-1.;//nodes[1] = Point (1.,-1.,-1.);
      nodes[6] =1.;  nodes[7] =1.;nodes[8] =-1.;//nodes[2] = Point (1.,1.,-1.);
       nodes[9] =-1.;  nodes[10] =1.;nodes[11] =-1.;// nodes[3] = Point (-1.,1.,-1.);
      nodes[12] =-1.;  nodes[13] =-1.;nodes[14] =1.;// nodes[4] = Point (-1.,-1.,1.); 
      nodes[15] =1.;  nodes[16] =-1.;nodes[17] =1.;//nodes[5] = Point (1.,-1.,1.);
      nodes[18] =1.;  nodes[19] =1.;nodes[20] =1.;//  nodes[6] = Point (1.,1.,1.);  
      nodes[21] =-1.;  nodes[22] =1.;nodes[23] =1.;//nodes[7] = Point (-1.,1.,1.);
        return;}
      
    case PRISM21:   {
         nodes[60] =1./3.;  nodes[61] =1./3.;nodes[62] =0.;}// nodes[20] = Point (1/Real(3),1/Real(3),0);
    case PRISM20:   {
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[18] = Point (1/Real(3),1/Real(3),-1); 
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;}//nodes[19] = Point (1/Real(3),1/Real(3),1);
    case PRISM18:      {
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[15] = Point (.5,0.,0.);  
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[16] = Point (.5,.5,0.);    
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;}//nodes[17] = Point (0.,.5,0.);}
    case PRISM15:  {
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[6] = Point (.5,0.,-1.);  
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[7] = Point (.5,.5,-1.); 
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[8] = Point (0.,.5,-1.);
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[9] = Point (0.,0.,0.);   
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[10] = Point (1.,0.,0.);   
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[11] = Point (0.,1.,0.);
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[12] = Point (.5,0.,1.);   
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;//nodes[13] = Point (.5,.5,1.);   
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =-1.;}//nodes[14] = Point (0.,.5,1.);   
    case PRISM6:    {
        nodes[0] =0.;  nodes[1] =0.;nodes[2] =-1.;//nodes[0] = Point (0.,0.,-1.); 
        nodes[3] =1.;  nodes[4] =0.;nodes[5] =-1.;//nodes[1] = Point (1.,0.,-1.);    
        nodes[6] =0.;  nodes[7] =1.;nodes[8] =-1.;//nodes[2] = Point (0.,1.,-1.);     
        nodes[9] =0.;  nodes[10] =0.;nodes[11] =1.;//nodes[3] = Point (0.,0.,1.);
        nodes[12] =1.;  nodes[13] =0.;nodes[14] =1.;//nodes[4] = Point (1.,0.,1.);   
        nodes[15] =0.;  nodes[16] =1.;nodes[17] =1.;//nodes[5] = Point (0.,1.,1.);
        return;}
        
    case PYRAMID18:{ // triangle centers
        nodes[42] =-2./3.;  nodes[43] =0.;nodes[44] =1./3.;//nodes[14] = Point (-2/Real(3),0.,1/Real(3)); 
        nodes[45] =0.;  nodes[46] =2./3.;nodes[47] =1./3.;//nodes[15] = Point (0.,2/Real(3),1/Real(3));   
        nodes[48] =2./3.;  nodes[49] =0.;nodes[50] =1./3.;//nodes[16] = Point (2/Real(3),0.,1/Real(3));
        nodes[51] =0.;  nodes[52] =-2/3.;nodes[53] =1./3.;//nodes[17] = Point (0.,-2/Real(3),1/Real(3));
      }
    case PYRAMID14:      {    // base center
        nodes[39] =0.;  nodes[40] =0.;nodes[41] =0.;//nodes[13] = Point (0.,0.,0.);
      }
    case PYRAMID13:    {  // base midedge
        nodes[15] =0.;  nodes[16] =-1.;nodes[17] =0.;//nodes[5] = Point (0.,-1.,0.);    
        nodes[18] =1.;  nodes[19] =0.;nodes[20] =0.;//nodes[6] = Point (1.,0.,0.);   
        nodes[21] =0.;  nodes[22] =1.;nodes[23] =0.;//nodes[7] = Point (0.,1.,0.);  
        nodes[24] =-1.;  nodes[25] =0.;nodes[26] =0.;//nodes[8] = Point (-1,0.,0.);
        // lateral midedge
        nodes[27] =-.5;  nodes[28] =-.5;nodes[29] =.5;//nodes[9] = Point (-.5,-.5,.5); 
        nodes[30] =.5;  nodes[31] =-.5;nodes[32] =.5;//nodes[10] = Point (.5,-.5,.5);  
        nodes[33] =.5;  nodes[34] =.5;nodes[35] =.5;//nodes[11] = Point (.5,.5,.5);  
        nodes[36] =-.5;  nodes[37] =.5;nodes[38] =.5;//nodes[12] = Point (-.5,.5,.5);
        }
    case PYRAMID5:  {
        // base corners
        nodes[0] =-1.;  nodes[1] =-1.;nodes[2] =0.;//nodes[0] = Point (-1.,-1.,0.); 
        nodes[3] =1.;  nodes[4] =-1.;nodes[5] =0.;//nodes[1] = Point (1.,-1.,0.);  
        nodes[6] =1.;  nodes[7] =1.;nodes[8] =0.;//nodes[2] = Point (1.,1.,0.);   
        nodes[9] =-1.;  nodes[10] =1.;nodes[11] =0.;//nodes[3] = Point (-1.,1.,0.);
        // apex
        nodes[12] =0.;  nodes[13] =0.;nodes[14] =1.;//nodes[4] = Point (0.,0.,1.);
        return;
      }
  
    default:   std::cerr<<"ERROR: Unknown element type " << itemType;
    }
}
// ===============================================================================
bool MGFE0::on_reference_element(
    const double p[], 
    const ElemType t, 
    const double eps)
{//=================================================================================
  assert (eps> 0.);
const double xi   = p[0];
const double eta  = p[1];
const double zeta = p[2];

  switch (t)    {
    case NODEELEM:      {        return (!xi && !eta && !zeta);      }
    case EDGE2:    case EDGE3:    case EDGE4:  {
        // The reference 1D element is [-1,1].
        if ((xi >= -1.-eps) && (xi <=  1.+eps))  return true;
        return false;}
    case TRI3:    case TRISHELL3:    case TRI6:    case TRI7:  {
        // The reference triangle is isosceles
        // and is bound by xi=0, eta=0, and xi+eta=1.
        if ((xi  >= 0.-eps) && (eta >= 0.-eps) && ((xi + eta) <= 1.+eps)) return true;
        return false;}
    case QUAD4:   case QUADSHELL4:    case QUAD8:    case QUADSHELL8:    case QUAD9:   case QUADSHELL9: {
        // The reference quadrilateral element is [-1,1]^2.
        if ((xi  >= -1.-eps) && (xi  <=  1.+eps) &&  (eta >= -1.-eps) && (eta <=  1.+eps)) return true;
        return false;
      }
    case TET4:    case TET10:    case TET14: {
        // The reference tetrahedral is isosceles and is bound by xi=0, eta=0, zeta=0, and xi+eta+zeta=1.
        if ((xi   >= 0.-eps) && (eta  >= 0.-eps) && (zeta >= 0.-eps) && ((xi + eta + zeta) <= 1.+eps))  return true;
        return false;
      }
    case HEX8:    case HEX20:    case HEX27:{
        /* if ((xi   >= -1.) &&  (xi   <=  1.) &&  (eta  >= -1.) && (eta  <=  1.) && (zeta >= -1.) && (zeta <=  1.))  return true; */

        // The reference hexahedral element is [-1,1]^3.
        if ((xi   >= -1.-eps) && (xi   <=  1.+eps) && (eta  >= -1.-eps) && (eta  <=  1.+eps) &&  (zeta >= -1.-eps) &&   (zeta <=  1.+eps))
          {
            //    libMesh::out << "Strange Point:\n";   //    p.print();
            return true;
          }
        return false;
      }
    case PRISM6:    case PRISM15:    case PRISM18:    case PRISM20:    case PRISM21:    {
        // Figure this one out...
        // inside the reference triangle with zeta in [-1,1]
        if ((xi   >=  0.-eps) &&
            (eta  >=  0.-eps) &&
            (zeta >= -1.-eps) &&
            (zeta <=  1.+eps) &&
            ((xi + eta) <= 1.+eps))
          return true;

        return false;
      }
    case PYRAMID5:    case PYRAMID13:    case PYRAMID14:    case PYRAMID18:      {
        // Check that the point is on the same side of all the faces by testing whether:
        // n_i.(x - x_i) <= 0 for each i, where:
        //   n_i is the outward normal of face i,
        //   x_i is a point on face i.
        if ((-eta - 1. + zeta <= 0.+eps) && (  xi - 1. + zeta <= 0.+eps) &&  ( eta - 1. + zeta <= 0.+eps) &&  ( -xi - 1. + zeta <= 0.+eps) &&
            (            zeta >= 0.-eps))
          return true;
        return false;
      }
    default:
      std::cerr<<"ERROR: Unknown element type " << t;
    }

  // If we get here then the point is _not_ in the
  // reference element.   Better return false.

  return false;
}


//=================================================
// void MGFE0::print_JxW(std::ostream & os) const{  std::abort();}//this->_fe_map->print_JxW(os);}
//=================================================
// void MGFE0::print_xyz(std::ostream & os) const{   std::abort();}//this->_fe_map->print_xyz(os);}

// //=================================================
// void MGFE0::print_info(std::ostream & os) const{
//   os << "phi[i][j]: Shape function i at quadrature pt. j" << std::endl;
//   this->print_phi(os);
// 
//   os << "dphi[i][j]: Shape function i's gradient at quadrature pt. j" << std::endl;
//   this->print_dphi(os);
// 
//   os << "XYZ locations of the quadrature pts." << std::endl;
//   this->print_xyz(os);
// 
//   os << "Values of JxW at the quadrature pts." << std::endl;
//   this->print_JxW(os);
// }

//=================================================
std::ostream & operator << (std::ostream & os, const MGFE0 & fe){  fe.print_info(os);  return os;}

//===========================================================
unsigned int MGFE0::n_quadrature_points () const
 {
   // if (this->shapes_on_quadrature)
//     {
//       // libmesh_assert(this->qrule);
//       // libmesh_assert_equal_to(this->qrule->n_points(),
//                               // this->_n_total_qp);
//     }
   return this->_n_total_qp;
}

