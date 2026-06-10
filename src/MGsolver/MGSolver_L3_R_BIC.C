
// ===============================================================
// --------------   Vector system with rotation R ------------------
// ===============================================================
//                             MGSolverR_BIC.C
// BCuniT  ->  B1) MGSolR::set_unit_vectors_along_stream_
//             B2) MGSolR::set_test_fun_bound_normals_R
//             B3) MGSolR::set_unit_vectors_on_boundary_R


// ===============================================================
// control files 
// #include "Printinfo_conf.h"  // Print options

// // local alg lib ----------------------------------------------->
// #include "sparse_matrixM.h"  // algebra sparse matrices
// #include "numeric_vectorM.h" // algebra numerical vectors
// #include "linear_solverM.h"  // algebra solvers
//Femus lib include -------------------------------------->
#include "MGSolver_L3_R.h"       // Navier-Stokes class header file
// // #include "MGFE_conf.h"        // FEM approximation
#include "MGGeom_L1_El.h"     
#include "MGMesh_L1_Extended.h"  
// // #include "MGFE.h"          // Mesh class
// #include "EquationSystemsExtendedM.h"  // Equation map class
// 






// ============================================================================================
void  MGSolR::set_unit_vectors_along_stream_R(
    const int i_step,
    const double vel[], ///< unit velocity
    double stream[],    ///< unit stream vector
    double normal[],    ///< unit normal
    double tg2[]        ///< second unit vector on tangent plane (3D)
) {
// if( _mgmesh._allow_Oxyz==1){
//     double modulus_vel=0.;
//     stream[2]=0.;
//     for(int ivar = 0; ivar < _ndim_B; ivar++) {
//         stream[ivar]=vel[ivar];
//         modulus_vel  += vel[ivar] *  vel[ivar];
//     }
//     modulus_vel = sqrt(modulus_vel);
//     if(modulus_vel < 1.e-20) {
//         modulus_vel = 1.;
//         stream[0] =  1.;
//         stream[1] = 0. ;
//     }
// 
//     tg2[2] = 0.;  // for 2D case attention tg -> norm !!!!!!!!!!!!!!!!
// // normal -----------------------------------------------------------------------------
// //     double norm_norm = 0.;  // always 3D
//     for(int ivar = 0; ivar < _ndim_B; ivar++) {
//         stream[ivar] =  stream[ivar]/modulus_vel;
//         tg2[ivar] =  stream[ivar]* stream[ivar];  // modulus normal squared cmps
//         //     norm_norm += tg_tmp[ivar];
//     }
// 
// // max normal x (i,j,k)
//     int i_max1 = 2;
//     double sum_norm_max=0.;// tg_tmp[0] + tg_tmp[1];
// //         double sum_norm[MGFE<2,LAGRANGE>::MDIM_FEM]; sum_norm_max[0] = tg_tmp[0] + tg_tmp[1];  // k-dir 2D case=sum_norm_k i_max1=2
//     for(int ij=0; ij<3; ij++) {
//         double  sum_norm= tg2[ij]+tg2[(ij+1) %_ndim_B];
//         if(sum_norm_max<sum_norm) {
//             i_max1= (ij+2) %3;
//             sum_norm_max=sum_norm;
//         }
//     }
// //         double sum_norm_j = tg_tmp[2] + tg_tmp[0];    // j-dir =sum_norm_j i_max1=1
// //         double sum_norm_i = tg_tmp[1] + tg_tmp[2];  // i-dir =sum_norm_i i_max1=0
// //         if(sum_norm_max < sum_norm_j) {i_max1=1; sum_norm_max=sum_norm_j;}
// //         if(sum_norm_max < sum_norm_i) {i_max1=0; sum_norm_max=sum_norm_i;}
//     sum_norm_max = sqrt(sum_norm_max);
// // tg1 -------------------------------------------------------------------------------
//     double sum_n=0.;
//     normal[(i_max1) % 3] = 0.;
//     for(int ij=1; ij<3; ij++) {
//         normal[(i_max1+ij) %3] = (-1* (2-ij)+1* (ij-1)) *stream[(i_max1+ij*2) %3]/sum_norm_max;
//         sum_n +=  normal[(i_max1+ij) %3]*  normal[(i_max1+ij) %3];
//     }
// //         normal[(i_max1+1)%3] = -stream[(i_max1+2)%3]/sum_norm_max;
// //         normal[(i_max1+2)%3] =  stream[(i_max1+1)%3]/sum_norm_max;
// // tg2 -------------------------------------------------------------------------------
// 
//     sum_n=sqrt(sum_n);
//     if(fabs(sum_n) <1.e-10) {
//         std::cout << "abort  point= ";
//         abort();
//     }
//     for(int ivar = 0; ivar < _ndim_B; ivar++) {
//         normal[ivar] /=sum_n;
//     }
// // tg_tmp=vel_tmp x  normal_tmp -> tg_tmp[k]=vel_tmp[k+1]*normal_tmp[k+2]-vel_tmp[k+2]*normal_tmp[k+1];
//     for(int ij=0; ij<3; ij++) {
//         tg2[ij] =  stream[(ij+1) %3]*normal[(ij+2) %3]- stream[(ij+2) %3]*normal[(ij+1) %3];
//     }
// //   if(i_step==0) {
     stream[0]= 1.;
      stream[1] = 0.;
      stream[2] = 0.;
      normal[0]= 0.;
      normal[1] = 1.;
      normal[2] = 0.;
      tg2[0]= 0.;
      tg2[1] = 0.;
      tg2[2]= 1.;
// //   }
// }
    return;

} //  for(int in = 0;

// int tab3d_pt_g[27]={};
// int tab2d_pt_g[9]={0,6,8,2,3,7,5,1,4}


// ================================================================================================
// This function computes boundary normals from each point test function
void MGSolR::set_test_fun_bound_normals_R(
    double normi_tmp[],    ///<  boundary normals
    int Coupled            ///< coupled or segregate for storage
) { // ============================================================================================
if( _mgmesh._allow_Oxyz==1){
    int corner[27]; // counter for mulitple points on boundary
    double xm[3];
    const int Level = _NoLevels - 1;
    const int el_ngauss = _fe_D[2]->get_NoGauss();  // elem gauss points
    const int n_nodes=_mgmesh._NoNodes[_NoLevels - 1];
//  _mgmesh._Oxyz_pts.clear();
    _mgmesh._Oxyz_pts.resize(n_nodes*_ndim_B*_ndim_B);

    for(int idim=0; idim<_top_offset_D*_ndim_B; idim++) normi_tmp[idim ]=0.;  // point normals


    int ndof_lev = 0;  // for multilevel multi proc point ordering
    for(int ipr = 0; ipr < _mgmesh._n_subdom; ipr++) {   // -----------------------------------------------------
        ndof_lev += _mgmesh._off_el[0][ipr * _NoLevels + Level + 1] - _mgmesh._off_el[0][ipr * _NoLevels + Level];
        const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * ipr + 1];  // start element
        const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * ipr];      // stop element

        for(int iel = 0; iel < (nel_e - nel_b); iel++) { // -----------------------------------
            _mgmesh.get_el_nod_conn(0, Level, iel, _el_conn_D, _xx_qnds_R,ipr);
            _mgmesh.get_el_neighbor(_el_sides_D, 0, Level, iel, _el_neigh_D,ipr);
            int compute=0;
            for(int i = 0; i < _el_dof_D[0]; i++) {
                int in1=_el_conn_D[i];
                for(int tn = 0; tn < _ndim_B; tn++) {
                    for(int idim = 0; idim < _ndim_B; idim++) {
                        _mgmesh._Oxyz_pts[ n_nodes*_ndim_B*tn+in1*_ndim_B+idim]=0.; // set up
//             _DAdata_eq_D[2].mg_eqs[_DAdata_eq_D[2].tab_eqs[NS_F]+ (1-Coupled) *idim ]->set_x_aux(tn,_node_dof[_NoLevels-1][in1+ Coupled*idim*_top_offset_D],0.);
                    }
                    _mgmesh._Oxyz_pts[ n_nodes*_ndim_B*tn+in1*_ndim_B+tn]=1.;// setup only diagonal matrix
//           _DAdata_eq_D[2].mg_eqs[_DAdata_eq_D[2].tab_eqs[NS_F]+ (1-Coupled)*tn ]->set_x_aux(tn,_node_dof[_NoLevels-1][in1+ Coupled*tn*_top_offset_D],1.);
                }
                corner[i]=0;
            }
            // set corner and compute flag values -------------------------------------------------
            for(int  iside = 0; iside < _el_sides_D; iside++) if(_el_neigh_D[iside] == -1) {
                    for(int idim=0; idim<_ndim_B; idim++) xm[idim] =0.;
                    for(int i = 0; i < _face_dof_D[0]; i++)    {
                        const int lnode = _mgmesh._GeomEl._surf_top(i + _face_dof_D[0]*iside);
                        for(int idim=0; idim<_ndim_B; idim++) xm[idim] +=_xx_qnds_R[lnode+idim*_el_dof_D[0]];
                        corner[lnode]++;
                    }
                    compute=1;
                } // ----------------------------------------------------------------------------------

            // compute normal  on boundary elements -----------------------------------------------
            // n_x=int dphi/dx dx; n_y=int dphi/dy dx; n_z=int dphi/dz dx;
            if(compute==1) {
                for(int qp = 0; qp < el_ngauss; qp++) { // start element integration
                    const double det2 = _fe_D[2]->Jac(qp, _xx_qnds_R, _InvJac2_R);         // quadratic Jacobian _InvJac2_R
                    const double JxW_g2 = det2 * _fe_D[2]->get_weight(qp);  // quadratic weight _fe_D[2]->_weight1
                    _fe_D[2]->get_dphi_gl_g(qp,_InvJac2_R,_dphi_g_D[2]);     // global coord deriv  ->  _dphi_g_D[2] der test fun (quad)
                    for(int i = 0; i < _el_dof_D[0]; i++)  if(corner[i] >0) for(int ivar = 0; ivar < _ndim_B; ivar++)
                                normi_tmp[ivar*_top_offset_D+_el_conn_D[i]] += JxW_g2*_dphi_g_D[2][i+ivar*_el_dof_D[0]]; // normal
                }// end elment integration
            }  // end of compute==1 ---------------------------------------------------------------

            // setting end point multiplicity (m) on bc[][]  (bc[][]=10000*m+bc[][]%10000) -------------------------------------------
            for(int i = 0; i < _el_dof_D[0]; i++) {
                if(corner[i]>1)   {
                    _bc[0][_el_conn_D[i]]=_bc[0][_el_conn_D[i]]%10000+corner[i]*10000;
                }
            } // ---------------------------------------------------------------------------------------------------------------------

        }//int iel ------------------------------------------------------------------------
    }//int ipr  --------------------------------------------------------------------------------------------------
//     for(int i = 0; i < _top_offset_D; i++)  {
// //         std::cout << i << "  ";
//         for(int ivar = 0; ivar < _ndim_B; ivar++) {
//             std::cout <<   normi_tmp[ivar*_top_offset_D+i] <<std::endl; // normal
//         }// end elment integration
// //         std::cout << std::endl;
//     }
//     std::cout << "ABORTTTTTTTTTTTTTTTTTTTTTTTTTT" <<std::endl;

//     abort();
}

    return;
}

//=======================================================================================================
void MGSolR::set_unit_vectors_on_boundary_R(
    double  normi_tmp[],
    const  int option,
    const  int itime_flag,
    int Coupled
) { //=======================================================================================================
if( _mgmesh._allow_Oxyz==1){
    const int Level = _NoLevels - 1;               //top level
    std::vector<int> el_dof_indices(_el_mat_nrows_R);// element dof vector
    double oxyz_tmp[3][3];                         // Coordinate system
    int flag_chtype[27];                      // flag inlet outlet
    double norm_gauss[3];         // normal from test function
    double norm_tmp0[3];                   // tmp side normal to the boundary
    double tg_tmp0[3];                     // tmp side normal to the boundary
    double tg_tmp1[3];
    double vel[3];
    double xf2[3];
    double val1[3];
    vel[2]=0.;
    double normal_surface[6*3];
    int    elb_conn[9];          ///< boundary element connectivity
    int bc_bd [27*(3+1)];
    int bc_vol[27*(3+1)];
    const int n_nodes=_mgmesh._NoNodes[_NoLevels - 1];

    int ndof_lev = 0;  // for multilevel multi proc point ordering -----------------------------------
    for(int ipr = 0; ipr < _mgmesh._n_subdom; ipr++) {
        ndof_lev += _mgmesh._off_el[0][ipr * _NoLevels + Level + 1] - _mgmesh._off_el[0][ipr * _NoLevels + Level];

        const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * ipr + 1];  // start element
        const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * ipr];      // stop element

        // *************************************************************************************
        // Volume  geometry and element  fields +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // get element Connectivity (_el_conn)  and coordinates (xx_qnds)
        // _el_conn  -> connectivity       // _xx_qnds_R -> coordinates
        // el_sides -> number of _n_sides  // _el_neigh -> neighbour element on a side
        // el_dof_indices -> local to global map
        for(int iel = 0; iel < (nel_e - nel_b); iel++) {

            _mgmesh.get_el_nod_conn(0, Level, iel, _el_conn_D, _xx_qnds_R,ipr);
            _mgmesh.get_el_neighbor(_el_sides_D, 0, Level, iel, _el_neigh_D,ipr);
            // set element-nodes variables  bc (bc_q_dofs)
//       get_el_dof_bc(Level,iel+ndof_lev,_el_dof_D,_el_conn_D,_top_offset_D,el_dof_indices,bc_vol,bc_bd);
            
        for(int i=0; i < 3; i++){
            for(int j=0; j < 3; j++) oxyz_tmp[j][i] = 0.;
        }

            // init _xyz_D[], flag_chtype[] val1[idim]
            for(int idim = 0; idim < _ndim_B; idim++) {    // quad loop entities (vector)
                _xyz_D[idim] = 0.;val1[idim]=0.;
                for(int d = 0; d < _el_dof_D[0]; d++) {
                    _xyz_D[idim] += _xx_qnds_R[idim * _el_dof_D[0]+ d];
                    flag_chtype[d]=0;
                }
                for(int kdim=0; kdim<_ndim_B; kdim++) {oxyz_tmp[kdim][idim] = 0.;    }
                _xyz_D[idim] /= _el_dof_D[0];
                oxyz_tmp[idim][idim] = 1.;
            }


            // change of type bc points -----------------------------------------------------------------------
            for(int  iside1=0; iside1 < _el_sides_D; iside1++) {
                if(_el_neigh_D[iside1] == -1) {                // if -1 the side is a boundary

                    int face_bc_flag=100000 ;  // face boundary condition flag = min flag
                    for(int  lbnode1 = 0; lbnode1 < _face_dof_D[0]; lbnode1++) {    // quad quantities
                        int lnode1 = _mgmesh._GeomEl._surf_top(lbnode1 + _face_dof_D[0] * iside1); // local nodes
                        _sur_toply_D[lbnode1] = lnode1;        // lbnode -> lnode
                        const  int bc_tmp=_bc[0][_el_conn_D[lnode1]];
                        if(face_bc_flag> bc_tmp) face_bc_flag=bc_tmp;  // face  bc (min bc_flag)
                    } // face bc minimum found (face_bc_flag)

                    for(int  lbnode2 = 0; lbnode2 < _face_dof_D[0]; lbnode2++) {
                        const int in=_el_conn_D[_sur_toply_D[lbnode2]];
                        const int bc_tmp1=_bc[0][in]%100;
                        if(bc_tmp1 != face_bc_flag%100)  {// bc change point
                            flag_chtype[_sur_toply_D[lbnode2]]=1;
                            std::cout << " bc type change  " << in;
                        }
                    }

                }
            }



            // change of type bc points -----------------------------------------------------------------------
            for(int  iside1=0; iside1 < _el_sides_D; iside1++) {
                if(_el_neigh_D[iside1] == -1) {

                    int face_bc_flag=100000 ;  // face boundary condition flag = min flag
                    for(int  lbnode1 = 0; lbnode1 < _face_dof_D[0]; lbnode1++) {    // quad quantities
                        int lnode1 = _mgmesh._GeomEl._surf_top(lbnode1 + _face_dof_D[0] * iside1); // local nodes
                        _sur_toply_D[lbnode1] = lnode1;        // lbnode -> lnode


                        for(int idim = 0; idim < _ndim_B; idim++)  {
                            _xxb_qnds_R[idim * _face_dof_D[0] + lbnode1] =_xx_qnds_R[idim * _el_dof_D[0] + lnode1];
                        }


                        const  int bc_tmp=_bc[0][_el_conn_D[lnode1]];
                        if(face_bc_flag> bc_tmp) face_bc_flag=bc_tmp;  // face  bc (min bc_flag)
                    } // face bc minimum found (face_bc_flag)



                    // adjustment inoutput ---------------------------------------------------------
                    const int face_bc_flag_normal=abs((face_bc_flag%100)/10);
                    if(face_bc_flag_normal> 3  && face_bc_flag_normal < 9) {}
                    else {

                        _fe_D[2]->Oxy_face(_xxb_qnds_R, _xyz_D, norm_tmp0,tg_tmp0,tg_tmp1);
                        for(int  lbnode2 = 0; lbnode2 < _face_dof_D[0]; lbnode2++)   if(flag_chtype[_sur_toply_D[lbnode2]]==1) {
                                int lnode3 = _mgmesh._GeomEl._surf_top(lbnode2 + _face_dof_D[0] * iside1); // local nodes
                                _sur_toply_D[lbnode2] = lnode3;        // lbnode -> lnode
                                elb_conn[lbnode2] = _el_conn_D[lnode3]; // connctivity _el_conn->elb_conn
                                const int in=_el_conn_D[_sur_toply_D[lbnode2]];
                                double sp=0.;
                                for(int isp=0; isp<_ndim_B; isp++)   sp +=norm_tmp0[isp]*normi_tmp[isp*_top_offset_D+elb_conn[lbnode2]];
                                double sum=0.;
                                for(int k=0; k<_ndim_B; k++) {
                                    const double val1a= norm_tmp0[k];// geometrical normal
                                    const double vala=normi_tmp[k*_top_offset_D+elb_conn[lbnode2]]-sp*val1a;
                                    normi_tmp[k*_top_offset_D+elb_conn[lbnode2]]=vala;
                                    sum +=vala*vala;
                                }
                                sum=sqrt(sum);
                                if(sum < 1.e-10) {
                                    std::cout<< "error: correction normal do not work " << sum ;/* abort();*/
//                                    for(int k=0; k<_ndim_B; k++)  normi_tmp[k*_top_offset_D+elb_conn[lbnode2]] =0.; normi_tmp[0*_top_offset_D+elb_conn[lbnode2]] =1.;
//                                      sum =1.;

                                }

                                for(int k=0; k<_ndim_B; k++) {
                                    normi_tmp[k*_top_offset_D+elb_conn[lbnode2]] /=sum;
                                }
                            }

                    }
                    // -------------------------------------------------------------------------------------


                }
            }








            // boundary surfaces --------------------------------------------------------------------------
            for(int  iside = 0; iside < _el_sides_D; iside++) {
                if(_el_neigh_D[iside] == -1) {

                    int face_bc_flag=100000;  // face boundary condition flag = min flag
                    for(int idim = 0; idim < _ndim_B; idim++)    xf2[idim] =0.;
                    // setup boundary element  ----------------------------------------------------------------
                    for(int  lbnode = 0; lbnode < _face_dof_D[0]; lbnode++) {    // quad quantities

                        int lnode = _mgmesh._GeomEl._surf_top(lbnode + _face_dof_D[0] * iside); // local nodes
                        _sur_toply_D[lbnode] = lnode;        // lbnode -> lnode
                        elb_conn[lbnode] = _el_conn_D[lnode]; // connctivity _el_conn->elb_conn

                        // norm_gauss <- normal from test functions ------------------------------------
                        double norm=0.;
                        for(int idim = 0; idim < _ndim_B; idim++)  {

                            _xxb_qnds_R[idim * _face_dof_D[0] + lbnode] =_xx_qnds_R[idim * _el_dof_D[0] + lnode];
                            const double tmp=normi_tmp[idim*_top_offset_D+elb_conn[lbnode]];
                            norm_gauss[idim+_ndim_B*lbnode]=tmp;
//                             xf2[idim] += _xxb_qnds_R[idim * _face_dof_D[0] + lbnode];
                            norm +=tmp*tmp;
                        }
                        norm=sqrt(norm);
                        if(norm<1.e-15) {
                            for(int idim = 0; idim < _ndim_B; idim++)  norm_gauss[idim+_ndim_B*lbnode]=0.;
                            norm_gauss[0+_ndim_B*lbnode]=1.;
                            norm=1.;//   std::cout << " no gaussiam normal";abort();
                        }
                        for(int idim = 0; idim < _ndim_B; idim++)  {
                            norm_gauss[idim+_ndim_B*lbnode] /=norm;
                        }
                        // ---------------------------------------------------------------------------
                        // face_bc_flag=face  bc (min bc_flag)
                        const  int bc_tmp=_bc[0][elb_conn[lbnode]];
                        if(face_bc_flag> bc_tmp) face_bc_flag=bc_tmp;  // face  bc (min bc_flag)
                        // ---------------------------------------------------------------------------
                    }
                    // face boundary conditions   face_bc_flag=min(bc_flag)
                    const int face_bc_flag_normal=abs((face_bc_flag%100)/10);

                    //  point tangent-stream direction
                    for(int  lbnode = 0; lbnode < _face_dof_D[0]; lbnode++) {
                        //             int compute_oxy=0;// 0-> not compute Oxyz  1-> compute Oxyz
                        int bc_shift=0; // 0-> normal=stream 1->stream=tg
                        const int in     =elb_conn[lbnode];
                        int flag_bc      =_bc[0][in]%1000;
                        const int iflag_normal=abs((flag_bc%100)/10);

                        // ========================== no flag_chtype points on  bc ====================
                        
                        if(flag_chtype[_sur_toply_D[lbnode]] ==1)
                        {
                            int a=2;
                            double b = 3;
                        }

                        if(flag_chtype[_sur_toply_D[lbnode]] ==0 || face_bc_flag_normal==8) {
                            for(int k=0; k<_ndim_B; k++) oxyz_tmp[2][k] = 0.; // zero tg2
                            // no inlet-outlet  (bc_shift=1; normal->n(1); Osnt -> s(0)n(1)t(2))
                            // inlet + outlet    (bc_shift=0; normal->s(0); Osnt -> s(0)n(1)t(2))
                            if(face_bc_flag_normal> 3  && face_bc_flag_normal < 9)   bc_shift=1;
                            else   bc_shift=0;   // inlet + outlet    (bc_shift=0; normal->s(0); Osnt -> s(0)n(1)t(2))
                            // setting the direction of the test function normals
//                             if(flag_chtype[_sur_toply_D[lbnode]] !=0) { // change bc point ----------------------
// //                                 _fe_D[2]->Oxy_face(_xxb_qnds_R, _x_m, norm_tmp0,tg_tmp0,tg_tmp1);
//                                 double sum=0.;
//                             for(int k=0; k<_ndim_B; k++)   {
//                                 norm_tmp0[k]= _xxb_qnds_R[k* _face_dof_D[0] +lbnode]-xf2[k];
//                                 sum +=   norm_tmp0[k]*norm_tmp0[k];
//                             }
//                             sum =sqrt(sum);double sp=0.;
//                             for(int isp=0;isp<_ndim_B;isp++) {
//                                  norm_tmp0[isp] /=sum;
//                                 sp +=norm_tmp0[isp]*norm_gauss[isp+_ndim_B*lbnode];
//                             }
//
//                                 for(int k=0; k<_ndim_B; k++) {
//                                    const double val1= norm_tmp0[k];// geometrical normal
//                                    const double val=norm_gauss[k+_ndim_B*lbnode]-0.*sp*val1;   oxyz_tmp[k][2]= val*val;
//                                    oxyz_tmp[k][bc_shift]=val;
//                                 }
//                             }// -------------------------------------------------------------------------------
//                             else
                            {
                                // standard point ----------------------------------------------------------
                                for(int k=0; k<_ndim_B; k++) {
                                    const  double val=norm_gauss[k+_ndim_B*lbnode]; // test function normal
                                    oxyz_tmp[k][bc_shift]=val;
                                    oxyz_tmp[k][2]= val*val;
                                }
                            }// ---------------------------------------------------------------------------------
                            //  bc_shift=1; normal -> 1; bc_shift=0; normal -> 0; Oxyz=s(0)n(1)t(2)
                            // max normal x (i,j,k); k-dir 2D case=sum_norm_k i_max1=2
                            int i_max1 = 2;
                            double sum_norm_max=0.*(oxyz_tmp[0][2]+oxyz_tmp[1][2]);
                            for(int ij=0; ij<=3; ij++) {
                                double  sum_norm= oxyz_tmp[ij][2]+oxyz_tmp[(ij+1)%_ndim_B][2];
                                if(sum_norm_max<sum_norm) {
                                    i_max1=(ij+2)%3;
                                    sum_norm_max=sum_norm;
                                }
                            }
                            sum_norm_max=sqrt(sum_norm_max);
                            // tg1 -------------------------------------------------------------------------------
                            double sum_n=0.;
                            oxyz_tmp[(i_max1)%3][(1+bc_shift)%_ndim_B]=0.;
                            val1[i_max1%3]=0.; //
                            for(int ij=1; ij<3; ij++) {
                                const double val=(-1*(2-ij)+1*(ij-1)) *oxyz_tmp[(i_max1+ij*2)%3][0+bc_shift]/sum_norm_max;
                                 val1[(i_max1+ij)%3]=val;
//                                 oxyz_tmp[(i_max1+ij)%3][(1+bc_shift)%_ndim_B]=val;
                                sum_n += val*val;//
                            }
                            sum_n=sqrt(sum_n);
                            for(int ivar=0; ivar<_ndim_B; ivar++) oxyz_tmp[ivar][(1+bc_shift)%_ndim_B] =val1[ivar]/sum_n;  //oxyz_tmp[ivar][1] ???
                            // tg2 ------------------------------------------------------------------------------
                            if(_ndim_B==3) { // only 3d --------------------------------------------- 3d
                                for(int ij=0; ij<3; ij++) { // oxyz_tmp[2] =oxyz_tmp[0] cross oxyz_tmp[[1]
                                    oxyz_tmp[ij][(2+bc_shift)%_ndim_B] = (oxyz_tmp[(ij+1)%3][(0+bc_shift)%_ndim_B]*oxyz_tmp[(ij+2)%3][(1+bc_shift)%_ndim_B]-oxyz_tmp[(ij+1)%3][(1+bc_shift)%_ndim_B]*oxyz_tmp[(ij+2)%3][(0+bc_shift)%_ndim_B]);
                                }
                            } // only 3d  ----------------------------------------------------------- 3d
//                            if(_AxiSym_D==1) if(_xx_qnds_R[0 * _el_dof_D[0] + _sur_toply_D[lbnode]]<0.001)
//                            {
//                                 for(int ivar = 0; ivar <  _ndim_B; ivar++)  {for(int tn = 0; tn <  _ndim_B; tn++) oxyz_tmp[ivar][tn]=0;oxyz_tmp[ivar][ivar]=1;}
//
//                             }

                                
                                // for(int ivar = 0; ivar < _ndim_B; ivar++)  for(int tn = 0; tn <  _ndim_B; tn++)
                                //     _mgmesh._Oxyz_pts[n_nodes*_ndim_B*tn+in*_ndim_B+ivar]=oxyz_tmp[ivar][tn];
  for(int ivar = 0; ivar < _ndim_B; ivar++)  for(int tn = 0; tn <  _ndim_B; tn++){
                                    _mgmesh._Oxyz_pts[n_nodes*_ndim_B*tn+in*_ndim_B+ivar]=0.;
                              if(     tn==ivar)  _mgmesh._Oxyz_pts[n_nodes*_ndim_B*tn+in*_ndim_B+ivar]=1.;
  }
 
//                   _DAdata_eq_D[2].mg_eqs[_DAdata_eq_D[2].tab_eqs[NS_F]+ (1-Coupled) *ivar ]->set_x_aux(tn,_node_dof[_NoLevels-1][in+ Coupled*ivar*_top_offset_D],oxyz_tmp[ivar][tn]);
                        }
                    }//int  lbnode = 0;
                } // iside -1
            }  // -----------------------------  End Boundary -------------------------------------

        }  //  =============== End of element loop =============================================
    }
    el_dof_indices.clear();
//     delete []normi_tmp; // temporary point test function normals
    DEBUG_1(CHECK_R1());

}
    return;
}
// //=======================================================================================================
// void MGSolR::set_unit_vectors_on_boundary(
//   const double  normi_tmp[],
//   const  int option,
//   const  int itime_flag
// ) { //=======================================================================================================
//
// if(_allow_Oxyz==1){
//   const int Level = _NoLevels - 1;               //top level
//   std::vector<int> el_dof_indices(_el_mat_nrows);// element dof vector
//   double oxyz_tmp[3][3];                         // Coordinate system
//   int flag_inout[NDOF_FEM];                      // flag inlet outlet
//   double norm_gauss[NDOF_FEM*MGFE<2,LAGRANGE>::MDIM_FEM];         // normal from test function
//   double norm_tmp0[MGFE<2,LAGRANGE>::MDIM_FEM];                   // tmp side normal to the boundary
//   double tg_tmp0[MGFE<2,LAGRANGE>::MDIM_FEM];                     // tmp side normal to the boundary
//   double tg_tmp1[MGFE<2,LAGRANGE>::MDIM_FEM]; double vel[3];vel[2]=0.;
//
//
//   int ndof_lev = 0;  // for multilevel multi proc point ordering -----------------------------------
//   for(int ipr = 0; ipr < _mgmesh._n_subdom; ipr++) {
//     ndof_lev += _mgmesh._off_el[0][ipr * _NoLevels + Level + 1] - _mgmesh._off_el[0][ipr * _NoLevels + Level];
//
//     const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * ipr + 1];  // start element
//     const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * ipr];      // stop element
//
//     // *************************************************************************************
//     // Volume  geometry and element  fields +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     // get element Connectivity (_el_conn)  and coordinates (xx_qnds)
//     // _el_conn  -> connectivity       // _xx_qnds_R -> coordinates
//     // el_sides -> number of _n_sides  // _el_neigh -> neighbour element on a side
//     // el_dof_indices -> local to global map
//     for(int iel = 0; iel < (nel_e - nel_b); iel++) {
//
//       _mgmesh.get_el_nod_conn(0, Level, iel, _el_conn, _xx_qnds_R,ipr);
//       _mgmesh.get_el_neighbor(_el_sides_D, 0, Level, iel, _el_neigh,ipr);
//       // set element-nodes variables  bc (bc_q_dofs)
//       get_el_dof_bc(Level,iel+ndof_lev,_el_dof,_el_conn,_top_offset,el_dof_indices,_bc_vol,bc_bd);
//       get_el_vel(_el_ndof,_el_conn,_u_1txy);
//
//       // init _xyz_D[], flag_inout[] vel[idim]
//       for(int idim = 0; idim < _ndim_B; idim++) {    // quad loop entities (vector)
//         vel[idim]=0.;  _xyz_D[idim] = 0.;  for(int d = 0; d < _el_dof_D[0]; d++) {
//           _xyz_D[idim] += _xx_qnds_R[idim * _el_dof_D[0]+ d];
//           flag_inout[d]=0; vel[idim] += _u_1txy[idim*_el_dof_D[0]+d];
//         }
//         _xyz_D[idim] /= _el_dof_D[0]; vel[idim]/= _el_dof_D[0];
//       }
//
//
//       // boundary surfaces --------------------------------------------------------------------------
//       for(int  iside = 0; iside < _el_sides_D; iside++) {
//         if(_el_neigh[iside] == -1) {
//
//           int face_bc_flag=1000;  // face boundary condition flag = min flag
//           // setup boundary element  ----------------------------------------------------------------
//           for(int  lbnode = 0; lbnode < _face_dof_D[0]; lbnode++) {    // quad quantities
//
//             int lnode = _mgmesh._GeomEl._surf_top[lbnode + _face_dof_D[0] * iside]; // local nodes
//             _sur_toply_D[lbnode] = lnode;        // lbnode -> lnode
//             elb_conn[lbnode] = _el_conn_D[lnode]; // connctivity _el_conn->elb_conn
//
//             // norm_gauss <- normal from test functions ------------------------------------
//             double norm=0.;
//             for(int idim = 0; idim < _ndim_B; idim++)  {
//               _xxb_qnds_R[idim * _face_dof_D[0] + lbnode] =_xx_qnds_R[idim * _el_dof_D[0] + lnode];
//               const double tmp=normi_tmp[idim*_top_offset+elb_conn[lbnode]];
//               norm_gauss[idim+_ndim_B*lbnode]=tmp;
//               norm +=tmp*tmp;
//             }
//             norm=sqrt(norm); if(norm<1.e-15) {std:cout << " no gaussiam normal"; abort();}
//             for(int idim = 0; idim < _ndim_B; idim++)  norm_gauss[idim+_ndim_B*lbnode] /=norm;
//             // ---------------------------------------------------------------------------
//             // face  bc (min bc_flag)
//             const  int bc_tmp=_bc[0][elb_conn[lbnode]];
//             if(face_bc_flag> bc_tmp) face_bc_flag=bc_tmp;  // face  bc (min bc_flag)
//               // ---------------------------------------------------------------------------
//           }
//           // face boundary conditions   face_bc_flag=min(bc_flag)
//           const int face_bc_flag_normal=abs((face_bc_flag%100)/10);
//
//           // data for in-out flow surfaces   (flag_inout=1)
//           for(int  lbnode = 0; lbnode < _face_dof_D[0]; lbnode++) {
//             const  int bc_flag=_bc[0][elb_conn[lbnode]];
//             if(bc_flag != face_bc_flag) {
//               if(face_bc_flag_normal < 4) flag_inout[_sur_toply_D[lbnode] ]=1;
//               if(bc_flag/10000>1) flag_inout[_sur_toply_D[lbnode] ]=0;
// //               if(face_bc_flag_normal==9)    {
// //                   flag_inout[_sur_toply_D[lbnode] ]=-1;
// //               }
//             }
//           }
//
//           _fe_D[2]->Oxy_face(_xxb_qnds_R, _x_m, norm_tmp0,tg_tmp0,tg_tmp1);
//
//
//           //  point tangent-stream direction
//           for(int  lbnode = 0; lbnode < _face_dof_D[0]; lbnode++) {
//
//             int compute_oxy=0;// 0-> not compute Oxyz  1-> compute Oxyz
//             int bc_shift=0; // 0-> normal=stream 1->stream=tg
//             const int in     =elb_conn[lbnode];
//             int flag_bc      =_bc[0][in]%1000;
//             const int iflag_normal=abs((flag_bc%100)/10);
//             int flag_io_surf = flag_inout[_sur_toply_D[lbnode] ];
//
//
//
//             // ========================== internal points on  boundary surfaces (flag_io_surf==0)====================
//
//             if(flag_bc==face_bc_flag && flag_io_surf==0) {  // bc (bc_flag) only over its own surface (min bc_flag)
//               compute_oxy=1; // normal -> 1 -> s(0)n(1)t(2)
//               int bc_sign=1-2*(int)(face_bc_flag_normal/9); // -1 only if  face_bc_flag_normal==9 (inlet)
//               for(int k=0; k<_ndim_B; k++) oxyz_tmp[2][k] = 0.;
//
//               // no inlet-outlet  (bc_shift=1; normal->n(1); Osnt -> s(0)n(1)t(2))
//               // inlet + outlet    (bc_shift=0; normal->s(0); Osnt -> s(0)n(1)t(2))
//               if(face_bc_flag_normal> 3  && face_bc_flag_normal < 9) { bc_shift=1;   }
//               else {  bc_shift=0;}
//
//               // setting the direction of the test function normals
//               for(int k=0; k<_ndim_B; k++) {
//                  double val=bc_sign*norm_gauss[k+_ndim_B*lbnode];
//                 if(_bc[1][_node_dof[_NoLevels-1][in]]/10000>1) {
//                     val= norm_tmp0[k];
//                 }
//                 oxyz_tmp[k][bc_shift]=val; oxyz_tmp[k][2]= val*val;
//               }
//
//
//               //  bc_shift=1; normal -> 1; bc_shift=0; normal -> 0; Oxyz=s(0)n(1)t(2)
//               // max normal x (i,j,k); k-dir 2D case=sum_norm_k i_max1=2
//               int i_max1 = 2; double sum_norm_max=0.;
//               for(int ij=0; ij<3; ij++) {
//                 double  sum_norm= oxyz_tmp[ij][2]+oxyz_tmp[(ij+1)%_ndim_B][2];
//                 if(sum_norm_max<sum_norm) {i_max1=(ij+2)%3; sum_norm_max=sum_norm;}
//               }
//               sum_norm_max=sqrt(sum_norm_max);
//               // tg1 -------------------------------------------------------------------------------
//               double sum_n=0.; oxyz_tmp[(i_max1)%3][(1+bc_shift)%_ndim_B]=0.;
//               for(int ij=1; ij<3; ij++) {
//                 const double val=(-1*(2-ij)+1*(ij-1)) *oxyz_tmp[(i_max1+ij*2)%3][0+bc_shift]/sum_norm_max;
//                 oxyz_tmp[(i_max1+ij)%3][(1+bc_shift)%_ndim_B]=val;
//
//                 sum_n += val*val;
//               }
// //               double scpr=0.;
// //               int ssign=1; if( face_bc_flag_normal>7 && scpr<-1.e-05 ) ssign=-1;
//               sum_n=sqrt(sum_n);  for(int ivar=0; ivar<_ndim_B; ivar++) {
//                   oxyz_tmp[ivar][(1+bc_shift)%_ndim_B] /=sum_n;  //oxyz_tmp[ivar][1] ???
//               }
//
//               // tg2 -------------------------------------------------------------------------------
//           double scpr0=0.;  for(int ivar=0; ivar<_ndim_B; ivar++) scpr0  +=oxyz_tmp[ivar][0]*vel[ivar];
//
//
//
//               if(_ndim_B==3) { // only 3d ------------------------------
//                 for(int ij=0; ij<3; ij++) { // oxyz_tmp[2] =oxyz_tmp[0] cross oxyz_tmp[[1]
//                   oxyz_tmp[ij][(2+bc_shift)%_ndim_B] = (oxyz_tmp[(ij+1)%3][(0+bc_shift)%_ndim_B]*oxyz_tmp[(ij+2)%3][(1+bc_shift)%_ndim_B]-oxyz_tmp[(ij+1)%3][(1+bc_shift)%_ndim_B]*oxyz_tmp[(ij+2)%3][(0+bc_shift)%_ndim_B]);
//                 }
//
//
//               //check always positive for tg direction on boundary -> -------------------------------
//
//             double scpr2=0.;
//                for(int ivar=0; ivar<_ndim_B; ivar++) scpr2  +=oxyz_tmp[ivar][2]*vel[ivar];
//
//                if(fabs(scpr2)>fabs(scpr0)) { for(int ivar=0; ivar<_ndim_B; ivar++){
//                  double tmp=oxyz_tmp[ivar][0];oxyz_tmp[ivar][0]=oxyz_tmp[ivar][2];oxyz_tmp[ivar][2]=tmp;
//                  }
//                  scpr0=scpr2;
//                }
//             } // only 3d  ------------------------------------------
//
//                 if(scpr0<-1.e-5) for(int ivar=0; ivar<_ndim_B; ivar++) {
//                   oxyz_tmp[ivar][0] *=-1;  //oxyz_tmp[ivar][1] ???
//               }
//
//               for(int ivar = 0; ivar < _ndim_B; ivar++)  for(int tn = 0; tn <  _ndim_B; tn++)
//                   _DAdata_eq_D[2].mg_eqs[_DAdata_eq_D[2].tab_eqs[NS_F]+ (1-_Coupled) *ivar ]->set_x_aux(tn,_node_dof[_NoLevels-1][in+ _Coupled*ivar*_top_offset],oxyz_tmp[ivar][tn]);
//             }
//             // ========================== end internal point on boundary surfaces ============================================
//             // =======================================================================================
//             // ==========================  outlet contour  ============================================
//
//             if(flag_bc !=face_bc_flag && flag_io_surf == 1) {  //  outlet contour (dim-2)
//
//               oxyz_tmp[2][0]=0.; oxyz_tmp[2][2] = 0.;
// //              double sum1=0.;
//               double sum0=0.;
//               for(int k=0; k<_ndim_B; k++) {
//                 oxyz_tmp[k][1]=flag_io_surf*norm_gauss[k+_ndim_B*lbnode]; // gaussian normal -> normal
//                 oxyz_tmp[k][0]=norm_tmp0[k]; // surface normal -> stream (s)
//                 sum0 +=oxyz_tmp[k][0]*oxyz_tmp[k][0];
// //                sum1 +=oxyz_tmp[k][1]*oxyz_tmp[k][1];
//               }
//               // correction to guassian normal   n=n'-(s.n')s (orthogonal to stream (s))
//               double prd=0.; double tmp0=0.; sum0=sqrt(sum0);
//               for(int k=0; k<_ndim_B; k++) {
//                 oxyz_tmp[k][0] /=sum0;
//                 prd += oxyz_tmp[k][0]* oxyz_tmp[k][1];
//               }
//               for(int k=0; k<_ndim_B; k++) {
//                 oxyz_tmp[k][1] -= prd* oxyz_tmp[k][0];
//                 tmp0 += oxyz_tmp[k][1]* oxyz_tmp[k][1];
//               }
//               tmp0=sqrt(tmp0);  for(int k=0; k<_ndim_B; k++) oxyz_tmp[k][1] /=tmp0;
//
//               if(_ndim_B==3) {
//                 for(int ij=0; ij<3; ij++) oxyz_tmp[ij][2] = (oxyz_tmp[(ij+1)%3][0]*oxyz_tmp[(ij+2)%3][1]-oxyz_tmp[(ij+1)%3][1]*oxyz_tmp[(ij+2)%3][0]);
//               }
//               //check always positive for tg direction on boundary -> -------------------------------
//               for(int ivar = 0; ivar < _ndim_B; ivar++)  for(int tn = 0; tn <  _ndim_B; tn++)
//                   _DAdata_eq_D[2].mg_eqs[_DAdata_eq_D[2].tab_eqs[NS_F]+ (1-_Coupled) *ivar ]->set_x_aux(tn,_node_dof[_NoLevels-1][in+ _Coupled*ivar*_top_offset],oxyz_tmp[ivar][tn]);
//             }
//             // ============================  end  in-out surfaces  ===================================
//
//             // =======================================================================================
//             // ===================== Inlet Boundary countour flag_io_surf==-1 ========================
//             if(flag_bc!=face_bc_flag && flag_io_surf==-1) {
//               int sign_s=1;int flag_tg=0;
//               double prd0=0.; double prd1=0.;  for(int k=0; k<_ndim_B; k++) {
//                 prd0 +=vel[k]*tg_tmp0[k]; prd1 +=vel[k]*tg_tmp1[k];
//               }
//                if(prd0<0.) sign_s=-1;
//                if(1.e-5+fabs(prd0)<fabs(prd1)) {sign_s=1;flag_tg=1;  if(prd1<0.) sign_s=-1;        }
//               oxyz_tmp[2][0]=0.; oxyz_tmp[2][2] = 0.;
//               double sum0=0.;
//               for(int k=0; k<_ndim_B; k++) {
//                 oxyz_tmp[k][1]=flag_io_surf*norm_gauss[k+_ndim_B*lbnode];
//                 oxyz_tmp[k][0]=(1-flag_tg)*tg_tmp0[k]+flag_tg*tg_tmp1[k]; // also 2d
//                 sum0 +=oxyz_tmp[k][0]*oxyz_tmp[k][0];
//               }
//
//               double prd=0.; double tmp0=0.; sum0=sign_s*sqrt(sum0);
//               for(int k=0; k<_ndim_B; k++) {
//                 oxyz_tmp[k][0] /=sum0;
//                 prd += oxyz_tmp[k][0]* oxyz_tmp[k][1];
//               }
//               for(int k=0; k<_ndim_B; k++) {
//                 oxyz_tmp[k][1] -= prd* oxyz_tmp[k][0];
//                 tmp0 += oxyz_tmp[k][1]* oxyz_tmp[k][1];
//               }
//               tmp0=sqrt(tmp0);
//               for(int k=0; k<_ndim_B; k++) oxyz_tmp[k][1] /=tmp0;
//               if(_ndim_B==3) {
//                 for(int ij=0; ij<3; ij++) oxyz_tmp[ij][2] = (oxyz_tmp[(ij+1)%3][0]*oxyz_tmp[(ij+2)%3][1]-oxyz_tmp[(ij+1)%3][1]*oxyz_tmp[(ij+2)%3][0]);
//               }
//               //check always positive for tg direction on boundary -> -------------------------------
//               for(int ivar = 0; ivar < _ndim_B; ivar++)  for(int tn = 0; tn <  _ndim_B; tn++)
//                   _DAdata_eq_D[2].mg_eqs[_DAdata_eq_D[2].tab_eqs[NS_F]+ (1-_Coupled) *ivar ]->set_x_aux(tn,_node_dof[_NoLevels-1][in+ _Coupled*ivar*_top_offset],oxyz_tmp[ivar][tn]);
//             }
//             // ========================Inlet Boundary countour flag_io_surf==-1======================
//
//
//           }//int  lbnode = 0;
//         } // iside -1
//
//
//       }  // -----------------------------  End Boundary -------------------------------------
//
//
//
//     }  //  =============== End of element loop =============================================
//   }
//   el_dof_indices.clear();
// //     delete []normi_tmp; // temporary point test function normals
//
// }
//   return;
// }
// #endif  //ENDIF NS_EQUATIONS
// #endif  // NS_equation is personal


void MGSolR::init_ntgbg() {
    double *normi_tmp=new double [_top_offset_D*_ndim_B];
    set_test_fun_bound_normals_R(normi_tmp,8);
    // set_unit_vectors_on_boundary_R(normi_tmp,0,0,8);
    delete[] normi_tmp;
      _mgmesh.active_Oxyz(0);
    return;
}












