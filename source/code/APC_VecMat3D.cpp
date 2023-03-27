//------------------------------------------------------------------------------
//
// APC_VecMat3D.cpp
// 
// Implements three-dimensional vectors, 3*3 matrices and
//          operators for vector/matrix arithmetics
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vcl.h>
#include <windows.h>    // Header file for windows
#include <stdio.h>      // Header file for standard Input/Output


#include "APC_Const.h"
#include "APC_VecMat3D.h"

#ifdef __GNUC__  // GNU C++ adaptation
#include "GNU_iomanip.h"
#endif

using namespace std;


//
// Implementation of class Polar 
//
Polar::Polar()
 : phi(0.0),
   theta(0.0),
   r(0.0)
{
}


Polar::Polar(double Az, double Elev, double R)
 : phi(Az),
   theta(Elev),
   r(R)
{
}


//
// Implementation of class Vec3D 
//

//
// Calculate polar components
//
void Vec3D::CalcPolarAngles ()
{
  // Length of projection in x-y-plane:
  const double rhoSqr = m_Vec[0] * m_Vec[0] + m_Vec[1] * m_Vec[1]; 
  
  // Norm of vector
  m_r = sqrt ( rhoSqr + m_Vec[2] * m_Vec[2] );
    
  // Azimuth of vector
  if ( (m_Vec[0]==0.0) && (m_Vec[1]==0.0) ) 
    m_phi = 0.0;
  else
    m_phi = atan2 (m_Vec[1], m_Vec[0]);
  if ( m_phi < 0.0 )  m_phi += 2.0*pi;
    
  // Altitude of vector
  const double rho = sqrt ( rhoSqr );
  if ( (m_Vec[2]==0.0) && (rho==0.0) ) 
    m_theta = 0.0;
  else
    m_theta = atan2(m_Vec[2], rho);
}


//
// Default constructor (used for creation of temporary Vec3D-objects)
//
Vec3D::Vec3D ()
 : m_r(0.0),
   m_theta(0.0),
   m_phi(0.0),
   m_bPolarValid(false)
{
   m_Vec[0] = 0.0;
   m_Vec[1] = 0.0;
   m_Vec[2] = 0.0;
}


//
// Constructor for component-wise initialization
//
Vec3D::Vec3D (double X, double Y, double Z)
 : m_r(0.0),
   m_theta(0.0),
   m_phi(0.0),
   m_bPolarValid(false)
{
   m_Vec[0] = X;
   m_Vec[1] = Y;
   m_Vec[2] = Z;
}


//
// Conversion constructor from polar components
//
Vec3D::Vec3D (const Polar& polar)
 : m_bPolarValid(true),
   m_theta(polar.theta),
   m_phi(polar.phi),
   m_r(polar.r)
{
  const double cosEl = cos(m_theta);

  m_Vec[0] = polar.r * cos(m_phi) * cosEl;
  m_Vec[1] = polar.r * sin(m_phi) * cosEl;
  m_Vec[2] = polar.r              * sin(m_theta);
}


//
// Component access (index may be "x", "y", or "z")
//
double Vec3D::operator [] (index Index) const
{
  return m_Vec[Index];
}


//
// Retrieval of polar components (index may be "r", "theta", or "phi")
//
double Vec3D::operator [] (pol_index Index)
{
  if ( !m_bPolarValid ) {
    // Recalculate norm and polar angles of vector
    CalcPolarAngles();
    m_bPolarValid = true;
  }
  
  switch ( Index ) {
    case r:
      return m_r;
      break;
    case theta:
      return m_theta;
      break;
    default:
      return m_phi;
  }
}


//
// In-place addition of vector
//
void Vec3D::operator += (const Vec3D& Vec)
{
  m_Vec[0] += Vec.m_Vec[0];
  m_Vec[1] += Vec.m_Vec[1];
  m_Vec[2] += Vec.m_Vec[2];
    
  m_bPolarValid = false;
}


//
// In-place subtraction of vector
//
void Vec3D::operator -= (const Vec3D& Vec)
{
  m_Vec[0] -= Vec.m_Vec[0];
  m_Vec[1] -= Vec.m_Vec[1];
  m_Vec[2] -= Vec.m_Vec[2];
  
  m_bPolarValid = false;
}


//
// Dot product
//
double Dot (const Vec3D& left, const Vec3D& right)
{
  return left.m_Vec[0] * right.m_Vec[0] +
         left.m_Vec[1] * right.m_Vec[1] +
         left.m_Vec[2] * right.m_Vec[2];
}


//
// Norm of vector
//
double Norm (const Vec3D& Vec)
{
  return sqrt(Dot(Vec, Vec));
}


//
// Scalar multiplication
//
Vec3D operator * (double fScalar, const Vec3D& Vec)
{
  return Vec3D ( fScalar * Vec.m_Vec[0], 
                 fScalar * Vec.m_Vec[1],
                 fScalar * Vec.m_Vec[2] );
}


//
// Scalar multiplication
//
Vec3D operator * (const Vec3D& Vec, double fScalar)
{
  return fScalar * Vec;
}


//
// Scalar division
//
Vec3D operator / (const Vec3D& Vec, double fScalar)
{
  return Vec3D ( Vec.m_Vec[0] / fScalar, 
                 Vec.m_Vec[1] / fScalar,
                 Vec.m_Vec[2] / fScalar );
}


//
// Unary minus of vector
//
Vec3D operator - (const Vec3D& Vec)
{
  return Vec3D( -Vec.m_Vec[0], -Vec.m_Vec[1], -Vec.m_Vec[2] );
}


//
// Addition of vectors
//
Vec3D operator + (const Vec3D& left, const Vec3D& right)   
{
  return Vec3D ( left.m_Vec[0] + right.m_Vec[0], 
                 left.m_Vec[1] + right.m_Vec[1], 
                 left.m_Vec[2] + right.m_Vec[2]  );
}


//
// Subtraction of vectors
//
Vec3D operator - (const Vec3D& left, const Vec3D& right)    
{
  return Vec3D ( left.m_Vec[0] - right.m_Vec[0], 
                 left.m_Vec[1] - right.m_Vec[1], 
                 left.m_Vec[2] - right.m_Vec[2]  );
}


//
// Vector product
//
Vec3D Cross(const Vec3D& left, const Vec3D& right)
{
  Vec3D vResult;
    
  vResult.m_Vec[0] = left.m_Vec[1] * right.m_Vec[2] - 
                     left.m_Vec[2] * right.m_Vec[1];
  vResult.m_Vec[1] = left.m_Vec[2] * right.m_Vec[0] - 
                     left.m_Vec[0] * right.m_Vec[2];
  vResult.m_Vec[2] = left.m_Vec[0] * right.m_Vec[1] - 
                     left.m_Vec[1] * right.m_Vec[0];
    
  return vResult;
}


//
// Simple vector output
//
ostream& operator << (ostream& os, const Vec3D& Vec)
{
  int w = os.width();

  os << setw(0) << " (";
  os << setw(w) << Vec.m_Vec[0] << setw(0) << ",";
  os << setw(w) << Vec.m_Vec[1] << setw(0) << ",";
  os << setw(w) << Vec.m_Vec[2];
  os << setw(0) << ")" << endl;

  return os;
}


//
// Implementation of class Mat3D 
//

//
// Default constructor (sets null matrix)
//
Mat3D::Mat3D ()
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      m_Mat[i][j] = 0.0;
}


//
// Constructor for matrix from column vectors
//
Mat3D::Mat3D (const Vec3D& e_1, const Vec3D& e_2, const Vec3D& e_3 )
{
  for (int i = 0; i < 3; i++) {
    m_Mat[i][0]=e_1.m_Vec[i];
    m_Mat[i][1]=e_2.m_Vec[i];
    m_Mat[i][2]=e_3.m_Vec[i];
  }
}


//
// Component access: retrieve given column of Matrix Mat
//
Vec3D Col(const Mat3D& Mat, index Index)
{
  Vec3D Res;

  for (int i = 0; i < 3; i++)
    Res.m_Vec[i] = Mat.m_Mat[i][int(Index)];

  return Res;

}


//
// Component access: retrieve given row of Matrix Mat
//
Vec3D Row(const Mat3D& Mat, index Index)
{
  Vec3D Res;

  for (int j = 0; j < 3; j++)
    Res.m_Vec[j] = Mat.m_Mat[int(Index)][j];
  
  return Res;
}


//
// Create identity matrix
//
Mat3D Id3D()
{
  Mat3D Id;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Id.m_Mat[i][j] = (i==j) ? 1 : 0;
   
  return Id;
}


//
// Elementary rotation matrices
//
Mat3D R_x(double RotAngle)
{
  const double S = sin (RotAngle);
  const double C = cos (RotAngle);

  Mat3D U;
  
  U.m_Mat[0][0] = 1.0;  U.m_Mat[0][1] = 0.0;  U.m_Mat[0][2] = 0.0;
  U.m_Mat[1][0] = 0.0;  U.m_Mat[1][1] =  +C;  U.m_Mat[1][2] =  +S;
  U.m_Mat[2][0] = 0.0;  U.m_Mat[2][1] =  -S;  U.m_Mat[2][2] =  +C;

  return U;
}


Mat3D R_y(double RotAngle)
{
  const double S = sin (RotAngle);
  const double C = cos (RotAngle);

  Mat3D U;
  
  U.m_Mat[0][0] =  +C;  U.m_Mat[0][1] = 0.0;  U.m_Mat[0][2] =  -S;
  U.m_Mat[1][0] = 0.0;  U.m_Mat[1][1] = 1.0;  U.m_Mat[1][2] = 0.0;
  U.m_Mat[2][0] =  +S;  U.m_Mat[2][1] = 0.0;  U.m_Mat[2][2] =  +C;

  return U;
}


Mat3D R_z(double RotAngle)
{
  const double S = sin (RotAngle);
  const double C = cos (RotAngle);

  Mat3D U;

  U.m_Mat[0][0] =  +C;  U.m_Mat[0][1] =  +S;  U.m_Mat[0][2] = 0.0;
  U.m_Mat[1][0] =  -S;  U.m_Mat[1][1] =  +C;  U.m_Mat[1][2] = 0.0;
  U.m_Mat[2][0] = 0.0;  U.m_Mat[2][1] = 0.0;  U.m_Mat[2][2] = 1.0;

  return U;
}


//
// Transpose of matrix
//
Mat3D Transp(const Mat3D& Mat)
{
  Mat3D T;

  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      T.m_Mat[i][j] = Mat.m_Mat[j][i];
    
  return T;
}


//
// Scalar multiplication
//
Mat3D operator * (double fScalar, const Mat3D& Mat)
{
  Mat3D Result;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Result.m_Mat[i][j] = fScalar*Mat.m_Mat[i][j];

  return Result;
}


Mat3D operator * (const Mat3D& Mat, double fScalar)
{
  return fScalar * Mat;
}


//
// Scalar division
//
Mat3D operator / (const Mat3D& Mat, double fScalar)
{
  Mat3D Result;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Result.m_Mat[i][j] = Mat.m_Mat[i][j]/fScalar;
  
  return Result;
}


//
// Matrix-vector product
//
Vec3D operator * (const Mat3D& Mat, const Vec3D& Vec)
{
  Vec3D Result;

  for (int i = 0; i < 3; i++) {
    double Scalp = 0.0;
    for (int j = 0; j < 3; j++)
      Scalp += Mat.m_Mat[i][j] * Vec.m_Vec[j];
    Result.m_Vec[i] = Scalp;
  }

  return Result;
}


//
// Vector-matrix product
//
Vec3D operator * (const Vec3D& Vec, const Mat3D& Mat)
{
  Vec3D Result;

  for (int j = 0; j < 3; j++) {
    double Scalp = 0.0;
    for (int i = 0; i < 3; i++)
      Scalp += Vec.m_Vec[i] * Mat.m_Mat[i][j] ;
    Result.m_Vec[j] = Scalp;
  }

  return Result;
}


//
// Unary minus of matrix
//
Mat3D operator - (const Mat3D& Mat)
{
  Mat3D Result;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Result.m_Mat[i][j] = -Mat.m_Mat[i][j];

  return Result;
}


//
// Addition of matrices
//
Mat3D operator + (const Mat3D& left, const Mat3D& right)
{
  Mat3D Result;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Result.m_Mat[i][j] = left.m_Mat[i][j] + right.m_Mat[i][j];
    }
  }
    
  return Result;
}


//
// Subtraction of matrices
//
Mat3D operator - (const Mat3D& left, const Mat3D& right)
{
  Mat3D Result;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Result.m_Mat[i][j] = left.m_Mat[i][j] - right.m_Mat[i][j];
    }
  }
    
  return Result;
}


//
// Matrix product
//
Mat3D operator * (const Mat3D& left, const Mat3D& right)
{
  Mat3D Result;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double Scalp = 0.0;
      for (int k = 0; k < 3; k++)
        Scalp += left.m_Mat[i][k] * right.m_Mat[k][j];
      Result.m_Mat[i][j] = Scalp;
    }
  }
    
  return Result;
}


//
// Simple matrix output
//
ostream& operator << (ostream& os, const Mat3D& Mat)
{
  int w = os.width();

  for (int i=0; i<3; i++) {
    os << setw(0) << " (";
    os << setw(w) << Mat.m_Mat[i][0] << setw(0) << ",";
    os << setw(w) << Mat.m_Mat[i][1] << setw(0) << ",";
    os << setw(w) << Mat.m_Mat[i][2];
    os << setw(0) << ")" << endl;
  }

  return os;
}
