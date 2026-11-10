/** This file is part of Corsair simulation.
 *
 *  Copyright 2011, 2012 Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <cmath>
#include <limits>

template<typename REAL> void crossProduct(const REAL* const vec1,const REAL* const vec2,REAL* const res);
template<int SIZE,typename REAL> REAL dotProduct(const REAL* const vec1,const REAL* const vec2);
template<int SIZE,typename REAL> void unitVector(REAL* vec);
template<int SIZE,typename REAL> REAL vectorMagnitude(const REAL* const vec);
template<int SIZE,typename REAL> REAL vectorMagnitude2(const REAL* const vec);

template<typename INT> INT matrixIndex(INT i,INT j);
template<typename REAL> void invertMatrix3by3(const REAL* const mat,REAL* invmat);
template<typename REAL> void matrixMatrixMultiply3by3(const REAL* const A,const REAL* const B,REAL* result);
template<int SIZE,typename REAL> void matrixVectorMultiply(const REAL* const mat,const REAL* const vec,REAL* result);

/** Calculate cross product between vectors vec1 and vec2, and store result to res.
 * @param vec1 First vector.
 * @param vec2 Second vector.
 * @param res Vector in which the result is written.*/
template<typename REAL> inline
void crossProduct(const REAL* const vec1,const REAL* const vec2,REAL* const res) {
   res[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
   res[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
   res[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

/** Calculate the dot product of two vectors.
 * @param vec1 First vector.
 * @param vec2 Second vector.
 * @return Dot product of vec1 and vec2.*/
template<int SIZE,typename REAL> inline
REAL dotProduct(const REAL* const vec1,const REAL* const vec2) {
   REAL result = 0.0;
   for (int i=0; i<SIZE; ++i) {
      result += vec1[i]*vec2[i];
   }
   return result;
}

/** Normalize a vector. Template parameter SIZE gives the size of the vector.
 * @param vec Vector to be normalized.*/
template<int SIZE,typename REAL> inline
void unitVector(REAL* vec) {
   REAL mag = 0.0;
   for (int i=0; i<SIZE; ++i) {
      mag += vec[i]*vec[i];
   }
   mag = sqrt(mag) + std::numeric_limits<REAL>::min();

   for (int i=0; i<SIZE; ++i) {
      vec[i] /= mag;
   }
}

template<int SIZE,typename REAL> inline
REAL vectorMagnitude(const REAL* const vec) {
   REAL mag = 0.0;
   for (int i=0; i<SIZE; ++i) {
      mag += vec[i]*vec[i];
   }
   return sqrt(mag);
}

template<int SIZE,typename REAL> inline
REAL vectorMagnitude2(const REAL* const vec) {
   REAL mag2 = 0.0;
   for (int i=0; i<SIZE; ++i) {
      mag2 += vec[i]*vec[i];
   }
   return mag2;
}

template<typename INT> inline 
INT matrixIndex(INT i,INT j) {
   return i*3+j;
}

/** Invert a 3x3 matrix.
 * @param mat Matrix to be inverted.
 * @param invmat Inverted matrix.*/
template<typename REAL> inline
void invertMatrix3by3(const REAL* const mat,REAL* invmat) {
   REAL det = mat[0]*mat[4]*mat[8] + mat[1]*mat[5]*mat[6];
   det += mat[2]*mat[3]*mat[7] - mat[0]*mat[5]*mat[7];
   det -= mat[1]*mat[3]*mat[8];
   det -= mat[2]*mat[4]*mat[6];

   invmat[0] = mat[4]*mat[8] - mat[5]*mat[7];
   invmat[3] = mat[5]*mat[6] - mat[3]*mat[8];
   invmat[6] = mat[3]*mat[7] - mat[4]*mat[6];

   invmat[1] = mat[2]*mat[7] - mat[1]*mat[8];
   invmat[4] = mat[0]*mat[8] - mat[2]*mat[6];
   invmat[7] = mat[1]*mat[6] - mat[0]*mat[7];

   invmat[2] = mat[1]*mat[5] - mat[2]*mat[4];
   invmat[5] = mat[2]*mat[3] - mat[0]*mat[5];
   invmat[8] = mat[0]*mat[4] - mat[1]*mat[3];
   for (int i=0; i<9; ++i) invmat[i] = invmat[i] / det;
}

template<typename REAL> inline
void matrixMatrixMultiply3by3(const REAL* const A,const REAL* const B,REAL* result) {
   result[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
   result[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
   result[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];

   result[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
   result[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
   result[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];

   result[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
   result[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
   result[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}

template<int SIZE,typename REAL> inline
void matrixVectorMultiply(const REAL* const mat,const REAL* const vec,REAL* result) {
   for (int i=0; i<SIZE; ++i) result[i] = 0.0;
   for (int j=0; j<SIZE; ++j) for (int i=0; i<SIZE; ++i) result[j] += mat[j*SIZE+i]*vec[i];
}

#endif
