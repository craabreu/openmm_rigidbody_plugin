#include "internal/eigenDecomposition.h"
#include "openmm/OpenMMException.h"
#include <math.h>
#include <limits>

/* References:
1) Oliver K. Smith
   Eigenvalues of a symmetric 3 × 3 matrix
   Communications of the ACM, 4 (1961), p. 168

2) Joachim Kopp
   Efficient numerical diagonalization of hermitian 3x3 matrices
   Int. J. Mod. Phys. C, 19 (2008), p. 523
*/

using namespace RigidBodyPlugin;
using namespace OpenMM;

// Constants
#define M_SQRT3    1.73205080756887729352744634151   // sqrt(3)
#define PI         3.14159265358979323846264338328
#define eps        std::numeric_limits<double>::epsilon()

// Macros
#define SQR(x)      ((x)*(x))                        // x^2

//--------------------------------------------------------------------------------------------------

void swap(double& a, double& b) {
    double c = a;
    a = b;
    b = c;
}

//--------------------------------------------------------------------------------------------------

void computeEigenvector( Vec3& q, Mat3 a, double n1tmp, double n2tmp, double thresh ) {
    double norm = q.dot(q);
    double n1 = n1tmp + a[0][0]*a[0][0];
    double n2 = n2tmp + a[1][1]*a[1][1];
    double error = n1*n2;

    // If the first column is zero, then (1, 0, 0) is an eigenvector
    if (n1 <= thresh)
        q = Vec3(1.0, 0.0, 0.0);

    // If the second column is zero, then (0, 1, 0) is an eigenvector
    else if (n2 <= thresh)
        q = Vec3(0.0, 1.0, 0.0);

    // If angle between A(*,1) and A(*,2) is too small, don't use
    //  cross product, but calculate v ~ (1, -A0/A1, 0)
    else if (norm < 4096.0*eps*eps*error) {
        double t = fabs(a[0][1]);
        double f = -a[0][0]/a[0][1];
        if (fabs(a[1][1]) > t) {
            t = fabs(a[1][1]);
            f = -a[0][1]/a[1][1];
        }
        if (fabs(a[1][2]) > t)
            f = -a[0][2]/a[1][2];
        norm = 1.0/sqrt(1.0 + f*f);
        q = Vec3(norm, f*norm, 0.0);
    }
    // This is the standard branch
    else
        q *= sqrt(1.0/norm);
}

//--------------------------------------------------------------------------------------------------

Vec3 eigenvalues(Mat3 A) {
    double p1 = SQR(A[0][1]) + SQR(A[0][2]) + SQR(A[1][2]);
    if (p1 < eps) {
        Vec3 w = A.diag();
        if (w[0] < w[1]) swap(w[0], w[1]);
        if (w[0] < w[2]) swap(w[0], w[2]);
        if (w[1] < w[2]) swap(w[1], w[2]);
        return w;
    }
    else {
        Vec3 d = A.diag();
        double TrA = d[0] + d[1] + d[2];
        double q = TrA/3.0;
        d -= Vec3(q, q, q);
        double p2 = d.dot(d) + 2.0*p1;
        double p = sqrt(p2/6.0);
        double r = (A - Diag3(q)).det()*(3.0/(p*p2));
        double phi;
        if (r <= -1.0)
            phi = PI/3.0;
        else if (r >= 1.0)
            phi = 0.0;
        else
            phi = acos(r)/3.0;
        double w0 = q + 2.0*p*cos(phi);
        double w2 = q + 2.0*p*cos(phi + 2.0*PI/3.0);
        return Vec3(w0, TrA - (w0 + w2), w2);
    }
}

//--------------------------------------------------------------------------------------------------

Mat3 eigenvectors(Mat3 A, Vec3& w) {
    double wmax8eps = 8.0*eps*fabs(w[0]);
    double thresh = wmax8eps*wmax8eps;

    // Prepare calculation of eigenvectors
    Vec3 q0, q1;
    Mat3 a = A.symmetric();
    double n1 = SQR(a[0][1]) + SQR(a[0][2]);
    double n2 = SQR(a[0][1]) + SQR(a[1][2]);
    q0[0] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
    q1[0] = q0[0];
    q0[1] = a[0][2]*a[0][1] - a[1][2]*a[0][0];
    q1[1] = q0[1];
    q1[2] = SQR(a[0][1]);

    // Calculate first eigenvector by the formula v(1) = (A - lambda(1)).e1 x (A - lambda(1)).e2
    a[0][0] -= w[0];
    a[1][1] -= w[0];
    q0 = Vec3(q1[0] + a[0][2]*w[0], q1[1] + a[1][2]*w[0], a[0][0]*a[1][1] - q1[2]);
    computeEigenvector( q0, a, n1, n2, thresh );

    // Prepare calculation of second eigenvector
    double  t = w[0] - w[1];

    // Is this eigenvalue degenerate?
    if (fabs(t) > wmax8eps) {
        // For non-degenerate eigenvalue, calculate second eigenvector by the formula
        //         v[1] = (A - lambda[1]).e1 x (A - lambda[1]).e2
        a[0][0] += t;
        a[1][1] += t;
        q1 = Vec3(q1[0] + a[0][2]*w[1], q1[1] + a[1][2]*w[1], a[0][0]*a[1][1] - q1[2]);
        computeEigenvector( q1, a, n1, n2, thresh );

    }
    else {
        // For degenerate eigenvalue, calculate second eigenvector according to
        //         v[1] = v(1) x (A - lambda[1]).e[i]

        // This would really get too complicated if we could not assume all of A to
        //       contain meaningful values.
        a[0][0] += w[0];
        a[1][1] += w[0];
        bool success = false;
        for (int i = 0; i < 3 && !success; i++) {
            a[i][i] -= w[1];
            Vec3 ai = a.col(i);
            n1 = ai.dot(ai);
            success = n1 > thresh;
            if (success) {
                q1 = q0.cross(ai);
                double norm = q1.dot(q1);
                success = norm > 65536.0*eps*eps*n1;
                if (success)
                    q1 *= sqrt(1.0/norm);
            }
        }

        // This means that any vector orthogonal to v(1) is an EV.
        if (!success) {
            int i = 0;
            while (q0[i] == 0.0)
                i++;
            int j = i % 3;
            double  norm = 1.0/sqrt(SQR(q0[i]) + SQR(q0[j]));
            q1[i] =  q0[j]*norm;
            q1[j] = -q0[i]*norm;
            q1[i+1 % 3] = 0.0;
        }
    }

    // Calculate third eigenvector according to v[2] = v(1) x v[1]
    return Mat3(q0, q1, q0.cross(q1));
}
