//------------------------------------------------------------------------------
//
// APC_Cheb.cpp
//
//    Implements three dimensional Chebyshev approximation
//
//    The original Fortran/Pascal/C++ codes is described in
//
//  - Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
//    Freeman and Comp., San Francisco (1975)
//  - Montenbruck O., Pfleger T.: Astronomie Mit Dem Personal Computer.
//    Springer-Verlag Berlin Heidelberg (1989).
//  - Montenbruck O., Pfleger T.: Astronomy on the Personal Computer.
//    Springer-Verlag Berlin Heidelberg (2000).
//
//------------------------------------------------------------------------------

//#include <cmath>

#include <math.h>
#include <vcl.h>
#include <windows.h>    // Header file for windows
#include <stdio.h>      // Header file for standard Input/Output


#include "APC_Cheb.h"
#include "APC_Const.h"
#include "APC_VecMat3D.h"

//
// Constructor
//
Cheb3D::Cheb3D(C3Dfunct f, int n, double dt) :
    m_n(n), m_f(f), m_dt(dt), m_Valid(false)
{
    Cx = new double[n + 1];
    Cy = new double[n + 1];
    Cz = new double[n + 1];
};

//
// Destructor
//
Cheb3D::~Cheb3D()
{
    delete[] Cx;
    delete[] Cy;
    delete[] Cz;
};

//
// Fit: approximates function m_f over the interval [ta, tb]
//
void Cheb3D::Fit(double ta, double tb)
{
    //
    // Variables
    //
    int j, k;
    double tau, t, fac;
    double* T = new double[m_n + 1];
    Vec3D f;

    // Clear all coefficients
    for (j = 0; j <= m_n; j++)
        Cx[j] = Cy[j] = Cz[j] = 0.0;

    // Loop over roots of T^(n+1)
    for (k = 0; k <= m_n; k++) {
		tau = Cos((2 * k + 1) * pi / (2 * m_n + 2)); // Grid point tau_k
		t = ((tb - ta) / 2.0) * tau + ((tb + ta) / 2.0); // Time t
		f = m_f(t); // Function value at time t

		// Calculate coefficients C_j using recurrence relation
		for (j = 0; j <= m_n; j++) {
			switch (j) {
				case 0:
					T[j] = 1.0;
					break;
				case 1:
					T[j] = tau;
					break;
				default:
                    T[j] = 2.0 * tau * T[j - 1] - T[j - 2];
            };
            // Increment coefficient C_j by f(t)*T_j(tau)
            Cx[j] += f[x] * T[j];
            Cy[j] += f[y] * T[j];
            Cz[j] += f[z] * T[j];
        };
    };

    // Scaling
    fac = 2.0 / (m_n + 1);
    for (j = 0; j <= m_n; j++) {
        Cx[j] *= fac;
        Cy[j] *= fac;
        Cz[j] *= fac;
    }

    // Update some members
    m_ta = ta;
    m_tb = tb;

    m_Valid = true;

    delete[] T;
}

//
// Value: evaluation of approximation at argument t
//
Vec3D Cheb3D::Value(double t)
{
    //
    // Constants
    //
    const double eps = 0.01; // Fractional overlap

    //
    // Variables
    //
    Vec3D f1, f2, old_f1;
    double tau, k;

    // Generate new coefficients as required
    if (!m_Valid || (t < m_ta) || (m_tb < t)) {
        k = floor(t / m_dt);
        Fit((k - eps) * m_dt, (k + 1 + eps) * m_dt);
    }

    // Evaluate approximation
    tau = (2.0 * t - m_ta - m_tb) / (m_tb - m_ta);

    for (int i = m_n; i >= 1; i--) {
        old_f1 = f1;
        f1 = 2.0 * tau * f1 - f2 + Vec3D(Cx[i], Cy[i], Cz[i]);
        f2 = old_f1;
    };

    return tau * f1 - f2 + 0.5 * Vec3D(Cx[0], Cy[0], Cz[0]);
}

