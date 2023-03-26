//------------------------------------------------------------------------------
//
// APC_VecMat3D.h
// 
// Classes for three-dimensional vectors and 3*3 matrices,
//          functions and operators for vector/matrix calculations
//
//------------------------------------------------------------------------------

#ifndef INC_APC_VECMAT3D_H
#define INC_APC_VECMAT3D_H

#include <iostream>


//
// Types
//

// identification for components of Vec3D objects
enum index { 
  x = 0,
  y = 1,
  z = 2
};


// declarations for polar angle representation of Vec3D objects
enum pol_index { 
  phi   = 0,  // azimuth of vector
  theta = 1,  // altitude of vector
  r     = 2   // norm of vector
};


//  Declaration of struct Polar
struct Polar {  

  // Constructors
  Polar();
  Polar(double Az, double Elev, double R = 1.0);

  // Members
  double phi;    // azimuth of vector
  double theta;  // altitude of vector
  double r;      // norm of vector
};


class Mat3D; // Forward declaration of class Mat3D used in class Vec3D


//
// Vec3D: three dimensional vectors
//
class Vec3D
{

  public:
    
    // Constructors
    Vec3D ();  // default constructor initializes to zero vector
    Vec3D (double X, double Y, double Z);
    Vec3D (const Polar& polar);
    
    // component access (read only)
    double operator [] (index Index) const;      
    
    // retrieves polar angles or norm of vector
    double operator [] (pol_index Index);         

    // component access
    friend Vec3D Col(const Mat3D& Mat, index Index);      
    friend Vec3D Row(const Mat3D& Mat, index Index);      

    // in-place addition of another vector
    void operator += (const Vec3D& Vec);
    
    // in-place subtraction of another vector
    void operator -= (const Vec3D& Vec);

    // dot product
    friend double Dot (const Vec3D& left, const Vec3D& right);

    // norm of vector
    friend double Norm (const Vec3D& Vec);

    // scalar multiplication
    friend Vec3D operator * (double fScalar, const Vec3D& Vec);
    friend Vec3D operator * (const Vec3D& Vec, double fScalar);
    
    // scalar division
    friend Vec3D operator / (const Vec3D& Vec, double fScalar);

    // unary minus of vector
    friend Vec3D operator - (const Vec3D& Vec);

    // addition of vectors
    friend Vec3D operator + (const Vec3D& left, const Vec3D& right);    
    
    // subtraction of vectors
    friend Vec3D operator - (const Vec3D& left, const Vec3D& right);    
    
    // vector product
    friend Vec3D Cross (const Vec3D& left, const Vec3D& right);

    // matrix-vector product
    friend Vec3D operator * (const Mat3D& Mat, const Vec3D& Vec);

    // vector-matrix product
    friend Vec3D operator * (const Vec3D& Vec, const Mat3D& Mat);

    // simple vector output
    friend std::ostream& operator << (std::ostream& os, const Vec3D& Vec);

    friend class Mat3D;

  private:
        
    // Members
    double m_Vec[3];        // components of vector
    double m_phi;           // polar angle (azimuth)
    double m_theta;         // polar angle (altitude)
    double m_r;             // norm of vector
    bool   m_bPolarValid;   // status flag for validity of polar coordinates
    
    // On-demand calculation of polar components
    void CalcPolarAngles ();  
};


//
// Mat3D: 3 dimensional transformation matrices
//
class Mat3D
{
  public:

    Mat3D (); // default constructor for null matrix
    
    // constructor for matrix from column vectors
    Mat3D ( const Vec3D& e_1, const Vec3D& e_2, const Vec3D& e_3 );

    // component access
    friend Vec3D Col(const Mat3D& Mat, index Index);      
    friend Vec3D Row(const Mat3D& Mat, index Index);      

    // identity matrix
    friend Mat3D Id3D();

    // elementary rotations
    friend Mat3D R_x(double RotAngle);
    friend Mat3D R_y(double RotAngle);
    friend Mat3D R_z(double RotAngle);

    // transposed matrix
    friend Mat3D Transp(const Mat3D& Mat);

    // scalar multiplication
    friend Mat3D operator * (double fScalar, const Mat3D& Mat);
    friend Mat3D operator * (const Mat3D& Mat, double fScalar);
    
    // scalar division
    friend Mat3D operator / (const Mat3D& Mat, double fScalar);

    // matrix-vector product
    friend Vec3D operator * (const Mat3D& Mat, const Vec3D& Vec);
    
    // vector-matrix product
    friend Vec3D operator * (const Vec3D& Vec, const Mat3D& Mat);

    // unary minus of matrix
    friend Mat3D operator - (const Mat3D& Mat);

    // addition of matrices
    friend Mat3D operator + (const Mat3D& left, const Mat3D& right);
    
    // subtraction of matrices
    friend Mat3D operator - (const Mat3D& left, const Mat3D& right);    

    // multiplication of matrices
    friend Mat3D operator * (const Mat3D& left, const Mat3D& right);
        
    // simple matrix output
    friend std::ostream& operator << (std::ostream& os, const Mat3D& Mat);
       
  private:

    double m_Mat[3][3];  // matrix elements
};


#endif  // include blocker
