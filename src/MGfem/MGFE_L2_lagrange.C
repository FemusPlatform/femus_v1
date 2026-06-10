// // The libMesh Finite Element Library.
// // Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
// 
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// 
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// // Lesser General Public License for more details.
// 
// // You should have received a copy of the GNU Lesser General Public
// // License along with this library; if not, write to the Free Software
// // Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// 
// 
// // Local includes
// #include "libmesh/dof_map.h"
// #include "libmesh/elem.h"
// #include "libmesh/enum_to_string.h"
#include "MGFE_L0_II_enum_order.h"
#include "MGFE_L2_FE_template.h"
// #include "libmesh/fe_interface.h"
// #include "libmesh/fe_macro.h"
// #include "libmesh/remote_elem.h"
// #include "libmesh/threads.h"
// 
#include<iostream>
// 
// #ifdef LIBMESH_ENABLE_AMR
// #ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
// // to have the macro 'inf_fe_switch' available.
// #include "libmesh/fe_interface_macros.h"
// #include "libmesh/inf_fe.h"
// #endif
// #endif
// 

// //==============================================================================
// int MGFE0::get_refspace_n_nodes(const ElemType itemType){
//    // nodes.resize(Elem::type_to_n_nodes_map[itemType]);
//   switch(itemType)  {
//     case NODEELEM: return 1;    
//     case EDGE3:    return 3;  
//     case EDGE2:    return 2;
//     case EDGE4:    return 4;
//     case TRI7: return 7;
//     case TRI6:  return 6;   
//     case TRI3:  return 3;
//     case TRISHELL3:   return 3;  
//     case QUAD9:  return 9;
//     case QUADSHELL9: return 9;
//     case QUAD8:  return 8;
//     case QUADSHELL8: return 8;
//     case QUAD4: return 4;
//     case QUADSHELL4:  return 4;
//     case TET14: return 14;
//     case TET10: return 10; 
//     case TET4:  return 4;
//     case HEX27: return 27;
//     case HEX20: return 20;
//     case HEX8: return 8;
//     case PRISM21: return 21;
//     case PRISM20:  return 20;
//     case PRISM18:return 18;
//     case PRISM15: return 15;
//     case PRISM6:return 6;
//     case PYRAMID18:return 18;
//     case PYRAMID14: return 14;
//     case PYRAMID13:return 13;
//     case PYRAMID5: return 5;
//     default:  std::cerr<<"ERROR: Unknown element type " << itemType; return -1;
//     }
// }
// namespace libMesh
// {
// 
// // global helper function
// void lagrange_nodal_soln(const Elem * elem,
//                          const Order order,
//                          const std::vector<Number> & elem_soln,
//                          std::vector<Number> &       nodal_soln,
//                          const bool add_p_level)
// {
//   const unsigned int n_nodes = elem->n_nodes();
//   const ElemType type        = elem->type();
// 
//   const Order totalorder = order + add_p_level*elem->p_level();
// 
//   nodal_soln.resize(n_nodes);
// 
// 
// 
//   switch (totalorder)
//     {
//       // linear Lagrange shape functions
//     case FIRST:
//       {
//         switch (type)
//           {
//           case EDGE3:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 2);
//               libmesh_assert_equal_to (nodal_soln.size(), 3);
// 
//               nodal_soln[0] = elem_soln[0];
//               nodal_soln[1] = elem_soln[1];
//               nodal_soln[2] = .5*(elem_soln[0] + elem_soln[1]);
// 
//               return;
//             }
// 
//           case EDGE4:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 2);
//               libmesh_assert_equal_to (nodal_soln.size(), 4);
// 
//               nodal_soln[0] = elem_soln[0];
//               nodal_soln[1] = elem_soln[1];
//               nodal_soln[2] = (2.*elem_soln[0] + elem_soln[1])/3.;
//               nodal_soln[3] = (elem_soln[0] + 2.*elem_soln[1])/3.;
// 
//               return;
//             }
// 
// 
//           case TRI7:
//             libmesh_assert_equal_to (nodal_soln.size(), 7);
//             nodal_soln[6] = (elem_soln[0] + elem_soln[1] + elem_soln[2])/3.;
//             libmesh_fallthrough();
//           case TRI6:
//             {
//               libmesh_assert (type == TRI7 || nodal_soln.size() == 6);
//               libmesh_assert_equal_to (elem_soln.size(), 3);
// 
//               nodal_soln[0] = elem_soln[0];
//               nodal_soln[1] = elem_soln[1];
//               nodal_soln[2] = elem_soln[2];
//               nodal_soln[3] = .5*(elem_soln[0] + elem_soln[1]);
//               nodal_soln[4] = .5*(elem_soln[1] + elem_soln[2]);
//               nodal_soln[5] = .5*(elem_soln[2] + elem_soln[0]);
// 
//               return;
//             }
// 
// 
//           case QUAD8:
//           case QUAD9:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 4);
// 
//               if (type == QUAD8)
//                 libmesh_assert_equal_to (nodal_soln.size(), 8);
//               else
//                 libmesh_assert_equal_to (nodal_soln.size(), 9);
// 
// 
//               nodal_soln[0] = elem_soln[0];
//               nodal_soln[1] = elem_soln[1];
//               nodal_soln[2] = elem_soln[2];
//               nodal_soln[3] = elem_soln[3];
//               nodal_soln[4] = .5*(elem_soln[0] + elem_soln[1]);
//               nodal_soln[5] = .5*(elem_soln[1] + elem_soln[2]);
//               nodal_soln[6] = .5*(elem_soln[2] + elem_soln[3]);
//               nodal_soln[7] = .5*(elem_soln[3] + elem_soln[0]);
// 
//               if (type == QUAD9)
//                 nodal_soln[8] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3]);
// 
//               return;
//             }
// 
// 
//           case TET14:
//             libmesh_assert_equal_to (nodal_soln.size(), 14);
//             nodal_soln[10] = (elem_soln[0] + elem_soln[1] + elem_soln[2])/3.;
//             nodal_soln[11] = (elem_soln[0] + elem_soln[1] + elem_soln[3])/3.;
//             nodal_soln[12] = (elem_soln[1] + elem_soln[2] + elem_soln[3])/3.;
//             nodal_soln[13] = (elem_soln[0] + elem_soln[2] + elem_soln[3])/3.;
//             libmesh_fallthrough();
//           case TET10:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 4);
//               libmesh_assert (type == TET14 || nodal_soln.size() == 10);
// 
//               nodal_soln[0] = elem_soln[0];
//               nodal_soln[1] = elem_soln[1];
//               nodal_soln[2] = elem_soln[2];
//               nodal_soln[3] = elem_soln[3];
//               nodal_soln[4] = .5*(elem_soln[0] + elem_soln[1]);
//               nodal_soln[5] = .5*(elem_soln[1] + elem_soln[2]);
//               nodal_soln[6] = .5*(elem_soln[2] + elem_soln[0]);
//               nodal_soln[7] = .5*(elem_soln[3] + elem_soln[0]);
//               nodal_soln[8] = .5*(elem_soln[3] + elem_soln[1]);
//               nodal_soln[9] = .5*(elem_soln[3] + elem_soln[2]);
// 
//               return;
//             }
// 
// 
//           case HEX20:
//           case HEX27:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 8);
// 
//               if (type == HEX20)
//                 libmesh_assert_equal_to (nodal_soln.size(), 20);
//               else
//                 libmesh_assert_equal_to (nodal_soln.size(), 27);
// 
//               nodal_soln[0]  = elem_soln[0];
//               nodal_soln[1]  = elem_soln[1];
//               nodal_soln[2]  = elem_soln[2];
//               nodal_soln[3]  = elem_soln[3];
//               nodal_soln[4]  = elem_soln[4];
//               nodal_soln[5]  = elem_soln[5];
//               nodal_soln[6]  = elem_soln[6];
//               nodal_soln[7]  = elem_soln[7];
//               nodal_soln[8]  = .5*(elem_soln[0] + elem_soln[1]);
//               nodal_soln[9]  = .5*(elem_soln[1] + elem_soln[2]);
//               nodal_soln[10] = .5*(elem_soln[2] + elem_soln[3]);
//               nodal_soln[11] = .5*(elem_soln[3] + elem_soln[0]);
//               nodal_soln[12] = .5*(elem_soln[0] + elem_soln[4]);
//               nodal_soln[13] = .5*(elem_soln[1] + elem_soln[5]);
//               nodal_soln[14] = .5*(elem_soln[2] + elem_soln[6]);
//               nodal_soln[15] = .5*(elem_soln[3] + elem_soln[7]);
//               nodal_soln[16] = .5*(elem_soln[4] + elem_soln[5]);
//               nodal_soln[17] = .5*(elem_soln[5] + elem_soln[6]);
//               nodal_soln[18] = .5*(elem_soln[6] + elem_soln[7]);
//               nodal_soln[19] = .5*(elem_soln[4] + elem_soln[7]);
// 
//               if (type == HEX27)
//                 {
//                   nodal_soln[20] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3]);
//                   nodal_soln[21] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[4] + elem_soln[5]);
//                   nodal_soln[22] = .25*(elem_soln[1] + elem_soln[2] + elem_soln[5] + elem_soln[6]);
//                   nodal_soln[23] = .25*(elem_soln[2] + elem_soln[3] + elem_soln[6] + elem_soln[7]);
//                   nodal_soln[24] = .25*(elem_soln[3] + elem_soln[0] + elem_soln[7] + elem_soln[4]);
//                   nodal_soln[25] = .25*(elem_soln[4] + elem_soln[5] + elem_soln[6] + elem_soln[7]);
// 
//                   nodal_soln[26] = .125*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3] +
//                                          elem_soln[4] + elem_soln[5] + elem_soln[6] + elem_soln[7]);
//                 }
// 
//               return;
//             }
// 
// 
//           case PRISM21:
//             nodal_soln[20]  = (elem_soln[9] + elem_soln[10] + elem_soln[11])/Real(3);
//             libmesh_fallthrough();
//           case PRISM20:
//             if (type == PRISM20)
//               libmesh_assert_equal_to (nodal_soln.size(), 20);
//             nodal_soln[18]  = (elem_soln[0] + elem_soln[1] + elem_soln[2])/Real(3);
//             nodal_soln[19]  = (elem_soln[3] + elem_soln[4] + elem_soln[5])/Real(3);
//             libmesh_fallthrough();
//           case PRISM18:
//             if (type == PRISM18)
//               libmesh_assert_equal_to (nodal_soln.size(), 18);
//             nodal_soln[15] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[4] + elem_soln[3]);
//             nodal_soln[16] = .25*(elem_soln[1] + elem_soln[2] + elem_soln[5] + elem_soln[4]);
//             nodal_soln[17] = .25*(elem_soln[2] + elem_soln[0] + elem_soln[3] + elem_soln[5]);
//             libmesh_fallthrough();
//           case PRISM15:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 6);
// 
//               if (type == PRISM15)
//                 libmesh_assert_equal_to (nodal_soln.size(), 15);
// 
//               nodal_soln[0]  = elem_soln[0];
//               nodal_soln[1]  = elem_soln[1];
//               nodal_soln[2]  = elem_soln[2];
//               nodal_soln[3]  = elem_soln[3];
//               nodal_soln[4]  = elem_soln[4];
//               nodal_soln[5]  = elem_soln[5];
//               nodal_soln[6]  = .5*(elem_soln[0] + elem_soln[1]);
//               nodal_soln[7]  = .5*(elem_soln[1] + elem_soln[2]);
//               nodal_soln[8]  = .5*(elem_soln[0] + elem_soln[2]);
//               nodal_soln[9]  = .5*(elem_soln[0] + elem_soln[3]);
//               nodal_soln[10] = .5*(elem_soln[1] + elem_soln[4]);
//               nodal_soln[11] = .5*(elem_soln[2] + elem_soln[5]);
//               nodal_soln[12] = .5*(elem_soln[3] + elem_soln[4]);
//               nodal_soln[13] = .5*(elem_soln[4] + elem_soln[5]);
//               nodal_soln[14] = .5*(elem_soln[3] + elem_soln[5]);
// 
//               return;
//             }
// 
//           case PYRAMID18:
//             {
//               libmesh_assert_equal_to (nodal_soln.size(), 18);
// 
//               nodal_soln[14] = (elem_soln[0] + elem_soln[1] + elem_soln[4])/Real(3);
//               nodal_soln[15] = (elem_soln[1] + elem_soln[2] + elem_soln[4])/Real(3);
//               nodal_soln[16] = (elem_soln[2] + elem_soln[3] + elem_soln[4])/Real(3);
//               nodal_soln[17] = (elem_soln[0] + elem_soln[3] + elem_soln[4])/Real(3);
// 
//               libmesh_fallthrough();
//             }
// 
//           case PYRAMID14:
//             {
//               if (type == PYRAMID14)
//                 libmesh_assert_equal_to (nodal_soln.size(), 14);
// 
//               nodal_soln[13] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3]);
// 
//               libmesh_fallthrough();
//             }
// 
//           case PYRAMID13:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 5);
// 
//               if (type == PYRAMID13)
//                 libmesh_assert_equal_to (nodal_soln.size(), 13);
// 
//               nodal_soln[0]  = elem_soln[0];
//               nodal_soln[1]  = elem_soln[1];
//               nodal_soln[2]  = elem_soln[2];
//               nodal_soln[3]  = elem_soln[3];
//               nodal_soln[4]  = elem_soln[4];
//               nodal_soln[5]  = .5*(elem_soln[0] + elem_soln[1]);
//               nodal_soln[6]  = .5*(elem_soln[1] + elem_soln[2]);
//               nodal_soln[7]  = .5*(elem_soln[2] + elem_soln[3]);
//               nodal_soln[8]  = .5*(elem_soln[3] + elem_soln[0]);
//               nodal_soln[9]  = .5*(elem_soln[0] + elem_soln[4]);
//               nodal_soln[10] = .5*(elem_soln[1] + elem_soln[4]);
//               nodal_soln[11] = .5*(elem_soln[2] + elem_soln[4]);
//               nodal_soln[12] = .5*(elem_soln[3] + elem_soln[4]);
// 
//               return;
//             }
//           default:
//             {
//               // By default the element solution _is_ nodal,
//               // so just copy it.
//               nodal_soln = elem_soln;
// 
//               return;
//             }
//           }
//       }
// 
//     case SECOND:
//       {
//         switch (type)
//           {
//           case EDGE4:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 3);
//               libmesh_assert_equal_to (nodal_soln.size(), 4);
// 
//               // Project quadratic solution onto cubic element nodes
//               nodal_soln[0] = elem_soln[0];
//               nodal_soln[1] = elem_soln[1];
//               nodal_soln[2] = (2.*elem_soln[0] - elem_soln[1] +
//                                8.*elem_soln[2])/9.;
//               nodal_soln[3] = (-elem_soln[0] + 2.*elem_soln[1] +
//                                8.*elem_soln[2])/9.;
//               return;
//             }
// 
//           case TRI7:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 6);
//               libmesh_assert_equal_to (nodal_soln.size(), 7);
// 
//               for (int i=0; i != 6; ++i)
//                 nodal_soln[i] = elem_soln[i];
// 
//               nodal_soln[6] = -1./9. * (elem_soln[0] + elem_soln[1] + elem_soln[2])
//                               +4./9. * (elem_soln[3] + elem_soln[4] + elem_soln[5]);
// 
//               return;
//             }
// 
//           case TET14:
//             {
//               libmesh_assert_equal_to (elem_soln.size(), 10);
//               libmesh_assert_equal_to (nodal_soln.size(), 14);
// 
//               for (int i=0; i != 10; ++i)
//                 nodal_soln[i] = elem_soln[i];
// 
//               nodal_soln[10] = -1./9. * (elem_soln[0] + elem_soln[1] + elem_soln[2])
//                                +4./9. * (elem_soln[4] + elem_soln[5] + elem_soln[6]);
//               nodal_soln[11] = -1./9. * (elem_soln[0] + elem_soln[1] + elem_soln[3])
//                                +4./9. * (elem_soln[4] + elem_soln[7] + elem_soln[8]);
//               nodal_soln[12] = -1./9. * (elem_soln[1] + elem_soln[2] + elem_soln[3])
//                                +4./9. * (elem_soln[5] + elem_soln[8] + elem_soln[9]);
//               nodal_soln[13] = -1./9. * (elem_soln[0] + elem_soln[2] + elem_soln[3])
//                                +4./9. * (elem_soln[6] + elem_soln[7] + elem_soln[9]);
// 
//               return;
//             }
// 
//           case PRISM21:
//             {
//               nodal_soln[20]  = (elem_soln[9] + elem_soln[10] + elem_soln[11])/Real(3);
//               libmesh_fallthrough();
//             }
//           case PRISM20:
//             {
//               if (type == PRISM20)
//                 libmesh_assert_equal_to (nodal_soln.size(), 20);
// 
//               for (int i=0; i != 18; ++i)
//                 nodal_soln[i] = elem_soln[i];
// 
//               nodal_soln[18]  = (elem_soln[0] + elem_soln[1] + elem_soln[2])/Real(3);
//               nodal_soln[19]  = (elem_soln[3] + elem_soln[4] + elem_soln[5])/Real(3);
//               return;
//             }
// 
//           case PYRAMID18:
//             {
//               libmesh_assert_equal_to (nodal_soln.size(), 18);
// 
//               for (int i=0; i != 14; ++i)
//                 nodal_soln[i] = elem_soln[i];
// 
//               nodal_soln[14] = (elem_soln[0] + elem_soln[1] + elem_soln[4])/Real(3);
//               nodal_soln[15] = (elem_soln[1] + elem_soln[2] + elem_soln[4])/Real(3);
//               nodal_soln[16] = (elem_soln[2] + elem_soln[3] + elem_soln[4])/Real(3);
//               nodal_soln[17] = (elem_soln[0] + elem_soln[3] + elem_soln[4])/Real(3);
// 
//               libmesh_fallthrough();
//             }
// 
//           default:
//             {
//               // By default the element solution _is_ nodal, so just
//               // copy the portion relevant to the nodal solution.
//               // (this is the whole nodal solution for true Lagrange
//               // elements, but a smaller part for "L2 Lagrange"
//               libmesh_assert_less_equal(nodal_soln.size(), elem_soln.size());
//               for (const auto i : index_range(nodal_soln))
//                 nodal_soln[i] = elem_soln[i];
// 
//               return;
//             }
//           }
//       }
// 
// 
// 
// 
//     default:
//       {
//         // By default the element solution _is_ nodal, so just copy
//         // the portion relevant to the nodal solution.  (this is the
//         // whole nodal solution for true Lagrange elements, but a
//         // smaller part for "L2 Lagrange"
//         libmesh_assert_less_equal(nodal_soln.size(), elem_soln.size());
//         for (const auto i : index_range(nodal_soln))
//           nodal_soln[i] = elem_soln[i];
// 
//         return;
//       }
//     }
// }
// 

//==============================================================================
int lagrange_refspace_n_nodes(const ElemType itemType,int geo_mesh){
   // nodes.resize(Elem::type_to_n_nodes_map[itemType]);
  switch( geo_mesh){
    case 0: //nodes in each Lagrange element
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
    case 1: //node in the  mesh geom element (example QUAD4->QUAD9(mesh)
     switch(itemType)  {
      case NODEELEM: return 1;    
      case EDGE3:   case EDGE2:  case EDGE4:    return 3;
      case TRI7:  case TRI6:  case TRI3:    case TRISHELL3:   return 6;  
      case QUAD9:   case QUADSHELL9:  case QUAD8:  case QUADSHELL8: case QUAD4:   case QUADSHELL4:  return 9;
      case TET14:  case TET10:  case TET4:  return 10;
    case HEX27: case HEX20:  case HEX8: return 27;
    case PRISM21:  case PRISM20:  case PRISM18:  case PRISM15:   case PRISM6:return 15;
    case PYRAMID18:  case PYRAMID14:   case PYRAMID13:  case PYRAMID5: return 14;
    default:  std::cerr<<"ERROR: Unknown element type " << itemType; return -1;
    }
    default:  std::cerr<<"ERROR: geom =0(elem) geom=1 (mesh) " << itemType; return -1;
  }
  return -1;
}


// ===============================================================
 unsigned int lagrange_n_dofs(
   const ElemType t, 
   const Order o
){// ===============================================================
  switch (o)  {
      // lagrange can only be constant on a single node
    case CONSTANT:  {   switch (t)  {
          case NODEELEM:        return 1;
          default:  std::cerr<<"ERROR: Bad ElemType = " << t << " for " << o << " order approximation!";
          }
      }
      // linear Lagrange shape functions
    case FIRST:   {   switch (t)     {
          case NODEELEM:  return 1;
          case EDGE2:  case EDGE3:   case EDGE4: return 2;
          case TRI3:  case TRISHELL3: case TRI3SUBDIVISION: case TRI6:  case TRI7:  return 3;
          case QUAD4: case QUADSHELL4:  case QUAD8:  case QUADSHELL8:  case QUAD9:  case QUADSHELL9:   return 4;
          case TET4:   case TET10:  case TET14:  return 4;
          case HEX8:  case HEX20:   case HEX27:  return 8;
          case PRISM6:  case PRISM15: case PRISM18:    case PRISM20: case PRISM21:     return 6;
          case PYRAMID5:  case PYRAMID13:   case PYRAMID14:   case PYRAMID18:  return 5;
          case INVALID_ELEM:return 0;
          default:   std::cerr<<"ERROR: Bad ElemType = " << t << " for " << o << " order approximation!";
          }
      }
      // quadratic Lagrange shape functions
    case SECOND: {     switch (t)   {
          case NODEELEM:  return 1;
          case EDGE3:  return 3;
          case TRI6:  case TRI7:   return 6;
          case QUAD8:  case QUADSHELL8: return 8;
          case QUAD9:  case QUADSHELL9:   return 9;
          case TET10:  case TET14:   return 10;
          case HEX20:  return 20;
          case HEX27:  return 27;
          case PRISM15:  return 15;
          case PRISM18:  case PRISM20: case PRISM21: return 18;    case PYRAMID13: return 13;
          case PYRAMID14:  case PYRAMID18: return 14;
          case INVALID_ELEM:   return 0;

          default:   std::cerr<<"ERROR: Bad ElemType = " << t << " for " << o << " order approximation!"; return 0;
          }
      }
    case THIRD:  {    switch (t)    {
          case NODEELEM:  return 1;
          case EDGE4: return 4;
          case PRISM20:  return 20;
          case PRISM21: return 21;
          case PYRAMID18: return 18;
          case TRI7: return 7;
          case TET14: return 14;
          case INVALID_ELEM:  return 0;
          default:  std::cerr<<"ERROR: Bad ElemType = " << t << " for " << o << " order approximation!";
          }
      }
    default: std::cerr<<"ERROR: Invalid Order " << o << " selected for LAGRANGE FE family!"; return 0;
    }
}
// 
// 
// 
// ==================================================================
unsigned int lagrange_n_dofs_at_node(
  const ElemType t,
  const Order o,
  const unsigned int n
){// =========================================================
  switch (o) { // lagrange can only be constant on a single node
    case CONSTANT:  {    switch (t)  {
          case NODEELEM:  return 1;
          default:  std::cerr<<"ERROR: Bad ElemType = " << t << " for " << o << " order approximation!";
          }
      }
      // linear Lagrange shape functions
    case FIRST:  {  switch (t)  {
          case NODEELEM: return 1;
          case EDGE2: case EDGE3: case EDGE4:   {   switch (n)  {
                case 0:    case 1:   return 1;
                default:  return 0;
                }
            }
          case TRI3:  case TRISHELL3: case TRI3SUBDIVISION: case TRI6:    case TRI7:   {   switch (n)   {
                case 0:    case 1:    case 2:   return 1;
                default:  return 0;
                }
            }
          case QUAD4: case QUADSHELL4:  case QUAD8:  case QUADSHELL8: case QUAD9:   case QUADSHELL9:   {   switch (n)    {
                case 0:  case 1:   case 2:    case 3:   return 1;
                default:  return 0;
                }
            }
          case TET4: case TET10:  case TET14:    {     switch (n)     {
                case 0:  case 1:  case 2: case 3:                  return 1;
                default:  return 0;
                }
            }
          case HEX8: case HEX20:  case HEX27: {      switch (n)      {
          case 0:  case 1:   case 2:   case 3:   case 4:    case 5:  case 6:   case 7:   return 1;
                default:  return 0;
                }
            }
          case PRISM6:  case PRISM15: case PRISM18: case PRISM20:  case PRISM21:     {     switch (n)       {
                case 0: case 1:   case 2:     case 3:   case 4:     case 5:        return 1;
                default:  return 0;
                }
            }
          case PYRAMID5: case PYRAMID13:  case PYRAMID14:   case PYRAMID18:      {      switch (n)           {
                case 0:    case 1:      case 2:       case 3:      case 4:    return 1;
                default: return 0;
                }
            }
          case INVALID_ELEM:  return 0;
          default:  std::cerr<<"ERROR: Bad ElemType = " << t << " for " << o << " order approximation!";
          }
      }
      // quadratic Lagrange shape functions
    case SECOND:      {        switch (t)          {
            // quadratic lagrange has one dof at each node
          case NODEELEM:    case EDGE3:   case TRI6:   case QUAD8:   case QUADSHELL8:    case QUAD9: case QUADSHELL9:          case TET10:
          case HEX20:     case HEX27:    case PRISM15:    case PRISM18:     case PYRAMID13:       case PYRAMID14:       return 1;
          case PRISM20:    case PRISM21:    return (n < 18);
          case PYRAMID18:   return (n < 14);
          case TRI7:      return (n < 6);
          case TET14:   return (n < 10);
          case INVALID_ELEM:  return 0;
          default: std::cerr<<"ERROR: Bad ElemType = " << t << " for " << o << " order approximation!";
          }
      }
    case THIRD:  {     switch (t)      {
          case NODEELEM:   case EDGE4: case PRISM20:   case PRISM21:  case PYRAMID18:  case TRI7:  case TET14:   return 1;
          case INVALID_ELEM:  return 0;
          default:   std::cerr<<"ERROR: Bad ElemType = " << t << " for " << o << " order approximation!";
          }
      }

    default:
      std::cerr<<"Unsupported order: "  << o ; return 0;
    }
}

// 
// } // anonymous namespace
// 
// 
// // Instantiate (side_) nodal_soln() function for every dimension
// LIBMESH_FE_NODAL_SOLN(LAGRANGE, lagrange_nodal_soln)
// LIBMESH_FE_SIDE_NODAL_SOLN(LAGRANGE)
//  template <> int MGFE<0,LAGRANGE>::get_refspace_n_nodes(const ElemType itemType, const int o){ return lagrange_refspace_n_nodes(itemType,o);}
//  template <> int MGFE<1,LAGRANGE>::get_refspace_n_nodes(const ElemType itemType, const int o){  return lagrange_refspace_n_nodes(itemType,o);}
//  template <> int MGFE<2,LAGRANGE>::get_refspace_n_nodes(const ElemType itemType, const int o){  return lagrange_refspace_n_nodes(itemType,o);}
//  template <> int MGFE<3,LAGRANGE>::get_refspace_n_nodes(const ElemType itemType, const int o){  return lagrange_refspace_n_nodes(itemType,o);}
// // 
// // Do full-specialization for every dimension, instead
// // of explicit instantiation at the end of this function.
// // This could be macro-ified.
 template <> unsigned int MGFE<0,LAGRANGE>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }
template <> unsigned int MGFE<1,LAGRANGE>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }
template <> unsigned int MGFE<2,LAGRANGE>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }
template <> unsigned int MGFE<3,LAGRANGE>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
template <> unsigned int MGFE<0,LAGRANGE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }
template <> unsigned int MGFE<1,LAGRANGE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }
template <> unsigned int MGFE<2,LAGRANGE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }
template <> unsigned int MGFE<3,LAGRANGE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }


// Lagrange elements have no dofs per element
// (just at the nodes)
template <> unsigned int MGFE<0,LAGRANGE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int MGFE<1,LAGRANGE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int MGFE<2,LAGRANGE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int MGFE<3,LAGRANGE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }

// Lagrange FEMs are always C^0 continuous
template <> FEContinuity MGFE<0,LAGRANGE>::get_continuity() const { return C_ZERO; }
template <> FEContinuity MGFE<1,LAGRANGE>::get_continuity() const { return C_ZERO; }
template <> FEContinuity MGFE<2,LAGRANGE>::get_continuity() const { return C_ZERO; }
template <> FEContinuity MGFE<3,LAGRANGE>::get_continuity() const { return C_ZERO; }

// Lagrange FEMs are not hierarchic
template <> bool MGFE<0,LAGRANGE>::is_hierarchic() const { return false; }
template <> bool MGFE<1,LAGRANGE>::is_hierarchic() const { return false; }
template <> bool MGFE<2,LAGRANGE>::is_hierarchic() const { return false; }
template <> bool MGFE<3,LAGRANGE>::is_hierarchic() const { return false; }

// Lagrange FEM shapes do not need reinit (is this always true?)
// template <> bool MGFE<0,LAGRANGE>::shapes_need_reinit() const { return false; }
// template <> bool MGFE<1,LAGRANGE>::shapes_need_reinit() const { return false; }
// template <> bool MGFE<2,LAGRANGE>::shapes_need_reinit() const { return false; }
// template <> bool MGFE<3,LAGRANGE>::shapes_need_reinit() const { return false; }

// // Methods for computing Lagrange constraints.  Note: we pass the
// // dimension as the last argument to the anonymous helper function.
// // Also note: we only need instantiations of this function for
// // Dim==2 and 3.
// #ifdef LIBMESH_ENABLE_AMR
// template <>
// void MGFE<2,LAGRANGE>::compute_constraints (DofConstraints & constraints,
//                                           DofMap & dof_map,
//                                           const unsigned int variable_number,
//                                           const Elem * elem)
// { lagrange_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/2); }
// 
// template <>
// void MGFE<3,LAGRANGE>::compute_constraints (DofConstraints & constraints,
//                                           DofMap & dof_map,
//                                           const unsigned int variable_number,
//                                           const Elem * elem)
// { lagrange_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/3); }
// #endif // LIBMESH_ENABLE_AMR
// 
// } // namespace libMesh
