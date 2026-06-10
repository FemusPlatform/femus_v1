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


// libMesh includes
#include "MGFE_L0_II_enum_order.h"
#include "MGquadrature_L1_gauss.h"
#include "MGenum_quadrature_type.h"

// namespace libMesh
// {

// See the files:
// quadrature_gauss_1D.C
// quadrature_gauss_2D.C
// quadrature_gauss_3D.C
// for implementation of specific element types.


QuadratureType QGauss::type() const
{
  return QGAUSS;
}

std::unique_ptr<QBase> QGauss::clone() const
{
  return std::make_unique<QGauss>(*this);
}

void QGauss::keast_rule(
  const double rule_data[][4],
  const unsigned int n_pts
){
  // Like the Dunavant rule, the input data should have 4 columns.  These columns
  // have the following format and implied permutations (w=weight).
  // {a, 0, 0, w} = 1-permutation  (a,a,a)
  // {a, b, 0, w} = 4-permutation  (a,b,b), (b,a,b), (b,b,a), (b,b,b)
  // {a, 0, b, w} = 6-permutation  (a,a,b), (a,b,b), (b,b,a), (b,a,b), (b,a,a), (a,b,a)
  // {a, b, c, w} = 12-permutation (a,a,b), (a,a,c), (b,a,a), (c,a,a), (a,b,a), (a,c,a)
  //                               (a,b,c), (a,c,b), (b,a,c), (b,c,a), (c,a,b), (c,b,a)

  // Always insert into the points & weights vector relative to the offset
  unsigned int offset=0;


  for (unsigned int p=0; p<n_pts; ++p)  {
      // There must always be a non-zero entry to start the row
      assert(rule_data[p][0] != static_cast<double>(0.0));
      // A zero weight may imply you did not set up the raw data correctly
      assert(rule_data[p][3]!= static_cast<double>(0.0));
      // What kind of point is this?
      // One non-zero entry in first 3 cols   ? 1-perm (centroid) point = 1
      // Two non-zero entries in first 3 cols ? 3-perm point            = 3
      // Three non-zero entries               ? 6-perm point            = 6
      unsigned int pointtype=1;
      if (rule_data[p][1] != static_cast<double>(0.0))   {
          if (rule_data[p][2] != static_cast<double>(0.0))   pointtype = 12;
          else  pointtype = 4;
        }
      else   { // The second entry is zero.  What about the third?
          if (rule_data[p][2] != static_cast<double>(0.0))   pointtype = 6;
        }


      switch (pointtype) {
        case 1:  {// Be sure we have enough space to insert this point
            assert(offset + 0< _n_points);
            const double a = rule_data[p][0];
            // The point has only a single permutation (the centroid!)
            for(int kdim=0;kdim<_dim;kdim++) _points[(offset  + 0)*_dim+kdim] = a;
            // The weight is always the last entry in the row.
            _weights[offset + 0] = rule_data[p][3];
            offset += pointtype;
            break;
          }

        case 4:  {
            // Be sure we have enough space to insert these points
            assert(offset + 3< _n_points);

            const double a  = rule_data[p][0];
            const double b  = rule_data[p][1];
            const double wt = rule_data[p][3];

            // Here it's understood the second entry is to be used twice, and
            // thus there are three possible permutations.
            _points[ (offset + 0)*_dim+0] = a;
            _points[ (offset + 1)*_dim+0] = b;
            _points[ (offset + 2)*_dim+0] = b;
            _points[ (offset + 3)*_dim+0] = b;
                                  
            _points[ (offset + 0)*_dim+1] = b;
            _points[ (offset + 1)*_dim+1] = a;
            _points[ (offset + 2)*_dim+1] = b;
            _points[ (offset + 3)*_dim+1] = b;
                                  
            _points[ (offset + 0)*_dim+2] = b;
            _points[ (offset + 1)*_dim+2] = b;
            _points[ (offset + 2)*_dim+2] = a;
            _points[ (offset + 3)*_dim+2] = b;

            for (unsigned int j=0; j<pointtype; ++j)  _weights[offset + j] = wt;
            offset += pointtype;
            break;
          }

        case 6:  {   // Be sure we have enough space to insert these points
            assert (offset + 5< _n_points);

            const double a  = rule_data[p][0];
            const double b  = rule_data[p][2];
            const double wt = rule_data[p][3];

            // Three individual entries with six permutations.
            _points[(offset + 0)*3+0] = a;
            _points[(offset + 1)*3+0] = a;
            _points[(offset + 2)*3+0] = b;
            _points[(offset + 3)*3+0] = b;
            _points[(offset + 4)*3+0] = b;
            _points[(offset + 5)*3+0] = a;
            
            _points[(offset + 0)*3+1] = a;
            _points[(offset + 1)*3+1] = b;
            _points[(offset + 2)*3+1] = b;
            _points[(offset + 3)*3+1] = a;
            _points[(offset + 4)*3+1] = a;
            _points[(offset + 5)*3+1] = b;
            
            _points[(offset + 0)*3+2] = b;
            _points[(offset + 1)*3+2] = b;
            _points[(offset + 2)*3+2] = a;
            _points[(offset + 3)*3+2] = b;
            _points[(offset + 4)*3+2] = a;
            _points[(offset + 5)*3+2] = a;
            
            for (unsigned int j=0; j<pointtype; ++j)  _weights[offset + j] = wt;
            offset += pointtype;
            break;
          }


        case 12:          {  // Be sure we have enough space to insert these points
            assert (offset + 11< _n_points);
            const double a  = rule_data[p][0];
            const double b  = rule_data[p][1];
            const double c  = rule_data[p][2];
            const double wt = rule_data[p][3];

            // Three individual entries with six permutations.
            _points[(offset + 0)*3+0] = a; _points[(offset + 6 )*3+0]  = a;
            _points[(offset + 1)*3+0] = a; _points[(offset + 7 )*3+0]  = a;
            _points[(offset + 2)*3+0] = b; _points[(offset + 8 )*3+0]  = b;
            _points[(offset + 3)*3+0] = c; _points[(offset + 9 )*3+0]  = b;
            _points[(offset + 4)*3+0] = a; _points[(offset + 10)*3+0] =  c;
            _points[(offset + 5)*3+0] = a; _points[(offset + 11)*3+0] =  c;
            
            
                _points[(offset + 0)*3+1] =a; _points[(offset + 6 )*3+1]  =b;
                _points[(offset + 1)*3+1] =a; _points[(offset + 7 )*3+1]  =c;
                _points[(offset + 2)*3+1] =a; _points[(offset + 8 )*3+1]  =a;
                _points[(offset + 3)*3+1] =a; _points[(offset + 9 )*3+1]  =c;
                _points[(offset + 4)*3+1] =b; _points[(offset + 10)*3+1] = a;
                _points[(offset + 5)*3+1] =c; _points[(offset + 11)*3+1] = b;
            
            
            
             _points[ (offset + 0)*3+2] =b; _points[(offset + 6 )*3+2]  = c;
             _points[ (offset + 1)*3+2] =c; _points[(offset + 7 )*3+2]  = b;
             _points[ (offset + 2)*3+2] =a; _points[(offset + 8 )*3+2]  = c;
             _points[ (offset + 3)*3+2] =a; _points[(offset + 9 )*3+2]  = a;
             _points[ (offset + 4)*3+2] =a; _points[(offset + 10)*3+2] =  b;
             _points[ (offset + 5)*3+2] =a; _points[(offset + 11)*3+2] =  a;
            

            for (unsigned int j=0; j<pointtype; ++j)   _weights[offset + j] = wt;

            offset += pointtype;
            break;
          }

        default:  std::cout<<"Don't know what to do with this many permutation points!"; abort();
        }

    }

}


// A number of different rules for triangles can be described by
// permutations of the following types of points:
// I:   "1"-permutation, (1/3,1/3)  (single point only)
// II:   3-permutation, (a,a,1-2a)
// III:  6-permutation, (a,b,1-a-b)
// The weights for a given set of permutations are all the same.
void QGauss::dunavant_rule2(const double * wts,
                            const double * a,
                            const double * b,
                            const unsigned int * permutation_ids,
                            unsigned int n_wts)
{
  // Figure out how many total points by summing up the entries
  // in the permutation_ids array, and resize the _points and _weights
  // vectors appropriately.
  unsigned int total_pts = 0;
  for (unsigned int p=0; p<n_wts; ++p) total_pts += permutation_ids[p];

  // Resize point and weight vectors appropriately.
   _n_points= total_pts;                     _points=new double[_n_points*_dim];  _weights=new double[_n_points];
  // _points.resize(total_pts);  _weights.resize(total_pts);

  // Always insert into the points & weights vector relative to the offset
  unsigned int offset=0;

  for (unsigned int p=0; p<n_wts; ++p)    {
      switch (permutation_ids[p])  {
        case 1:  {
            // The point has only a single permutation (the centroid!)
            // So we don't even need to look in the a or b arrays.
            // _points [offset  + 0] = Point(double(1)/3, double(1)/3);
            _points [(offset  + 0)*2+0] = double(1.)/3.;
            _points [(offset  + 0)*2+1] =  double(1.)/3.;
            _weights[offset + 0] = wts[p];

            offset += 1;
            break;
          }


        case 3: {
            // For this type of rule, don't need to look in the b array.
            // _points[offset + 0] = Point(a[p],         a[p]);         // (a,a)
            // _points[offset + 1] = Point(a[p],         1-2*a[p]); // (a,1-2a)
            // _points[offset + 2] = Point(1-2*a[p], a[p]);         // (1-2a,a)
            
_points[(offset + 0)*2+0] = a[p]   ;      // (a,a)
_points[(offset + 1)*2+0] = a[p]   ; // (a,1-2a)
_points[(offset + 2)*2+0] = 1-2*a[p];  // (1-2a,a)

_points[(offset + 0)*2+1] =     a[p];         // (a,a)
_points[(offset + 1)*2] =       1-2*a[p]; // (a,1-2a)
_points[(offset + 2)*2+1] = a[p];         // (1-2a,a)
            

            for (unsigned int j=0; j<3; ++j)  _weights[offset + j] = wts[p];
            offset += 3;
            break;
          }

        case 6:    {  // This type of point uses all 3 arrays...
            // _points[offset + 0] = Point(a[p], b[p]);
            // _points[offset + 1] = Point(b[p], a[p]);
            // _points[offset + 2] = Point(a[p], 1-a[p]-b[p]);
            // _points[offset + 3] = Point(1-a[p]-b[p], a[p]);
            // _points[offset + 4] = Point(b[p], 1-a[p]-b[p]);
            // _points[offset + 5] = Point(1-a[p]-b[p], b[p]);
            
            
             _points[(offset + 0)*2+0] =a[p];
             _points[(offset + 1)*2+0] =b[p];
             _points[(offset + 2)*2+0] =a[p];
             _points[(offset + 3)*2+0] =1-a[p]-b[p];
             _points[(offset + 4)*2+0] =b[p];
             _points[(offset + 5)*2+0] =1-a[p]-b[p];

              _points[(offset + 0)*2+1] =  b[p];
              _points[(offset + 1)*2+1] =  a[p];
              _points[(offset + 2)*2+1] =  1-a[p]-b[p];
              _points[(offset + 3)*2+1] =  a[p];
              _points[(offset + 4)*2+1] =  1-a[p]-b[p];
              _points[(offset + 5)*2+1] =  b[p];
             
             
             
             
             
            for (unsigned int j=0; j<6; ++j)  _weights[offset + j] = wts[p];

            offset += 6;
            break;
          }

        default:   std::cout<<"Unknown permutation id: " << permutation_ids[p] << "!";
        }
    }

}


void QGauss::dunavant_rule(const double rule_data[][4],
                           const unsigned int n_pts)
{
  // The input data array has 4 columns.  The first 3 are the permutation points.
  // The last column is the weights for a given set of permutation points.  A zero
  // in two of the first 3 columns implies the point is a 1-permutation (centroid).
  // A zero in one of the first 3 columns implies the point is a 3-permutation.
  // Otherwise each point is assumed to be a 6-permutation.

  // Always insert into the points & weights vector relative to the offset
  unsigned int offset=0;


  for (unsigned int p=0; p<n_pts; ++p)  {

      // There must always be a non-zero entry to start the row
      assert( rule_data[p][0]!= static_cast<double>(0.0) );

      // A zero weight may imply you did not set up the raw data correctly
      assert( rule_data[p][3]!= static_cast<double>(0.0) );

      // What kind of point is this?
      // One non-zero entry in first 3 cols   ? 1-perm (centroid) point = 1
      // Two non-zero entries in first 3 cols ? 3-perm point            = 3
      // Three non-zero entries               ? 6-perm point            = 6
      unsigned int pointtype=1;

      if (rule_data[p][1] != static_cast<double>(0.0))  {
          if (rule_data[p][2] != static_cast<double>(0.0))  pointtype = 6;
          else            pointtype = 3;
        }

      switch (pointtype)  {
        case 1:  { // Be sure we have enough space to insert this point
            assert (offset + 0< _n_points);

            // The point has only a single permutation (the centroid!)
            // _points[offset  + 0] = Point(rule_data[p][0], rule_data[p][0]);
             _points[(offset  + 0)*2+0] = rule_data[p][0];
             _points[(offset  + 0)*2+1] =  rule_data[p][0];
            // The weight is always the last entry in the row.
             _weights[offset + 0] = rule_data[p][3];

            offset += 1;
            break;
          }

        case 3:     {  // Be sure we have enough space to insert these points
            assert (offset + 2< _n_points);

            // Here it's understood the second entry is to be used twice, and
            // thus there are three possible permutations.
//             _points[offset + 0] = Point(rule_data[p][0], rule_data[p][1]);
//             _points[offset + 1] = Point(rule_data[p][1], rule_data[p][0]);
//             _points[offset + 2] = Point(rule_data[p][1], rule_data[p][1]);
//             
            _points[(offset + 0)*2+0] = rule_data[p][0]; 
            _points[(offset + 1)*2+0] = rule_data[p][1]; 
            _points[(offset + 2)*2+0] = rule_data[p][1]; 
            _points[(offset + 0)*2+1] =  rule_data[p][1];
            _points[(offset + 1)*2+1] =  rule_data[p][0];
            _points[(offset + 2)*2+1] =  rule_data[p][1];
            
            
            
            
            

            for (unsigned int j=0; j<3; ++j)   _weights[offset + j] = rule_data[p][3];

            offset += 3;
            break;
          }

        case 6:
          {
            // Be sure we have enough space to insert these points
            assert (offset + 5< _n_points);

            // Three individual entries with six permutations.
            // _points[offset + 0] = Point(rule_data[p][0], rule_data[p][1]);
            // _points[offset + 1] = Point(rule_data[p][0], rule_data[p][2]);
            // _points[offset + 2] = Point(rule_data[p][1], rule_data[p][0]);
            // _points[offset + 3] = Point(rule_data[p][1], rule_data[p][2]);
            // _points[offset + 4] = Point(rule_data[p][2], rule_data[p][0]);
            // _points[offset + 5] = Point(rule_data[p][2], rule_data[p][1]);
            
            
             _points[(offset + 0)*2+0] =rule_data[p][0];
             _points[(offset + 1)*2+0] =rule_data[p][0];
             _points[(offset + 2)*2+0] =rule_data[p][1];
             _points[(offset + 3)*2+0] =rule_data[p][1];
             _points[(offset + 4)*2+0] =rule_data[p][2];
             _points[(offset + 5)*2+0] =rule_data[p][2];
             _points[(offset + 0)*2+1] =rule_data[p][1];
             _points[(offset + 1)*2+1] =rule_data[p][2];
             _points[(offset + 2)*2+1] =rule_data[p][0];
             _points[(offset + 3)*2+1] =rule_data[p][2];
             _points[(offset + 4)*2+1] =rule_data[p][0];
             _points[(offset + 5)*2+1] =rule_data[p][1];
            
            

            for (unsigned int j=0; j<6; ++j)  _weights[offset + j] = rule_data[p][3];

            offset += 6;
            break;
          }

        default:
          std::cout<<"Don't know what to do with this many permutation points!";
        }
    }
}

// } // namespace libMesh
