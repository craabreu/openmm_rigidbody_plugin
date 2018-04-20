/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "internal/MatVec.h"
#include <math.h>

using namespace RigidBodyPlugin;
using namespace OpenMM;

/*------------------------------------------------------------------------------
 FULL 3x3 MATRICES
------------------------------------------------------------------------------*/

// Create a Mat3 whose elements are all 0
Mat3::Mat3() {
    row.resize(3);
    row[0] = row[1] = row[2] = Vec3();
}

// Create a Mat3 with specified x, y, and z rows
Mat3::Mat3(Vec3 x, Vec3 y, Vec3 z) {
    row.resize(3);
    row[0] = x;
    row[1] = y;
    row[2] = z;
}

// Extract a row
Vec3 Mat3::operator[](int index) const {
    assert(index >= 0 && index < 3);
    return row[index];
}

// Extract the diagonal
Vec3 Mat3::diag() const {
    return Vec3(row[0][0], row[1][1], row[2][2]);
}

// Assign row from Vec3
Vec3& Mat3::operator[](int index) {
    assert(index >= 0 && index < 3);
    return row[index];
}

// Extract column
Vec3 Mat3::col(int index) const {
    assert(index >= 0 && index < 3);
    return Vec3(row[0][index], row[1][index], row[2][index]);
}

// Equality
bool Mat3::operator==(const Mat3& rhs) const {
    return (row[0] == rhs[0] && row[1] == rhs[1] && row[2] == rhs[2]);
}

// Inequality
bool Mat3::operator!=(const Mat3& rhs) const {
    return (row[0] != rhs[0] || row[1] != rhs[1] || row[2] != rhs[2]);
}

// Unary plus
Mat3 Mat3::operator+() const {
    return Mat3(*this);
}

// Plus
Mat3 Mat3::operator+(const Mat3& rhs) const {
    const Mat3& lhs = *this;
    return Mat3(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
}

// Plus equal
Mat3& Mat3::operator+=(const Mat3& rhs) {
    row[0] += rhs[0];
    row[1] += rhs[1];
    row[2] += rhs[2];
    return *this;
}

// Unary minus
Mat3 Mat3::operator-() const {
    const Mat3& lhs = *this;
    return Mat3(-lhs[0], -lhs[1], -lhs[2]);
}

// Minus
Mat3 Mat3::operator-(const Mat3& rhs) const {
    const Mat3& lhs = *this;
    return Mat3(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
}

Mat3 Mat3::operator-(const Diag3& rhs) const {
    Mat3 a = *this;
    Vec3 b = rhs.diag();
    a[0][0] -= b[0];
    a[1][1] -= b[1];
    a[2][2] -= b[2];
    return a;
}

// Minus equal
Mat3& Mat3::operator-=(const Mat3& rhs) {
    row[0] -= rhs[0];
    row[1] -= rhs[1];
    row[2] -= rhs[2];
    return *this;
}

// Product with scalar
Mat3 Mat3::operator*(double rhs) const {
    const Mat3& lhs = *this;
    return Mat3(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs);
}

// Times equal
Mat3& Mat3::operator*=(double rhs) {
    row[0] *= rhs;
    row[1] *= rhs;
    row[2] *= rhs;
    return *this;
}

// Division by scalar
Mat3 Mat3::operator/(double rhs) const {
    const Mat3& lhs = *this;
    double scale = 1.0/rhs;
    return Mat3(lhs[0]*scale, lhs[1]*scale, lhs[2]*scale);
}

// Equal divided
Mat3& Mat3::operator/=(double rhs) {
    double scale = 1.0/rhs;
    row[0] *= scale;
    row[1] *= scale;
    row[2] *= scale;
    return *this;
}

// Symmetrization
Mat3 Mat3::symmetric() const {
    Mat3 a = *this;
    a[1][0] = a[0][1];
    a[2][0] = a[0][2];
    a[2][1] = a[1][2];
    return a;
}

// Trace
double Mat3::trace() const {
    const Mat3& lhs = *this;
    return lhs[0][0] + lhs[1][1] + lhs[2][2];
}

// Determinant
double Mat3::det() const {
    const Mat3& lhs = *this;
    const Vec3& a0 = lhs.row[0];
    const Vec3& a1 = lhs.row[1];
    const Vec3& a2 = lhs.row[2];
    return a0[0]*(a1[1]*a2[2] - a2[1]*a1[2]) -
           a0[1]*(a1[0]*a2[2] - a2[0]*a1[2]) +
           a0[2]*(a1[0]*a2[1] - a2[0]*a1[1]);
}

// Transposition
Mat3 Mat3::t() const {
    const Mat3& lhs = *this;
    return Mat3(lhs.col(0), lhs.col(1), lhs.col(2));
}

// Product with vector
Vec3 Mat3::operator*(const Vec3& rhs) const {
    const Mat3& lhs = *this;
    return Vec3(lhs[0].dot(rhs), lhs[1].dot(rhs), lhs[2].dot(rhs));
}

// Product with another full matrix
Mat3 Mat3::operator*(const Mat3& rhs) const {
    const Mat3& lhs = *this;
    return Mat3(lhs*(rhs.col(0)), lhs*(rhs.col(1)), lhs*rhs.col(2)).t();
}

// Product with a diagonal matrix
Mat3 Mat3::operator*(const Diag3& rhs) const {
    const Mat3& lhs = *this;
    return Mat3(rhs*lhs[0], rhs*lhs[1], rhs*lhs[2]);
}

/*------------------------------------------------------------------------------
 DIAGONAL 3x3 MATRICES
------------------------------------------------------------------------------*/

// Create a Diag3 whose elements are all 0
Diag3::Diag3() {
    data = Vec3();
}

// Create a Diag3 with specified x, y, and z rows
Diag3::Diag3(Vec3 x) {
    data = x;
}

// Create a Diag3 with specified x in all diagonal entries
Diag3::Diag3(double x) {
    data = Vec3(x, x, x);
}

// Extract a row
Vec3 Diag3::operator[](int index) const {
    assert(index >= 0 && index < 3);
    Vec3 x;
    x[index] = data[index];
    return x;
}

// Extract the diagonal
Vec3 Diag3::diag() const {
    return data;
}

// Equality
bool Diag3::operator==(const Diag3& rhs) const {
    return (data == rhs.data);
}

// Inequality
bool Diag3::operator!=(const Diag3& rhs) const {
    return (data != rhs.data);
}

// Unary plus
Diag3 Diag3::operator+() const {
    return Diag3(*this);
}

// Plus
Diag3 Diag3::operator+(const Diag3& rhs) const {
    const Diag3& lhs = *this;
    return Diag3(lhs.data + rhs.data);
}

// Plus equal
Diag3& Diag3::operator+=(const Diag3& rhs) {
    data += rhs.data;
    return *this;
}

// Unary minus
Diag3 Diag3::operator-() const {
    const Diag3& lhs = *this;
    return Diag3(-lhs.data);
}

// Minus
Diag3 Diag3::operator-(const Diag3& rhs) const {
    const Diag3& lhs = *this;
    return Diag3(lhs.data - rhs.data);
}

// Minus equal
Diag3& Diag3::operator-=(const Diag3& rhs) {
    data -= rhs.data;
    return *this;
}

// Product with scalar
Diag3 Diag3::operator*(double rhs) const {
    const Diag3& lhs = *this;
    return Diag3(lhs.data*rhs);
}

// Times equal
Diag3& Diag3::operator*=(double rhs) {
    data *= rhs;
    return *this;
}

// Division by scalar
Diag3 Diag3::operator/(double rhs) const {
    const Diag3& lhs = *this;
    return Diag3(lhs.data/rhs);
}

// Equal divided
Diag3& Diag3::operator/=(double rhs) {
    data /= rhs;
    return *this;
}

// Transposition
Diag3 Diag3::t() const {
    return *this;
}

// Product with a vector
Vec3 Diag3::operator*(const Vec3& rhs) const {
    const Diag3& lhs = *this;
    return Vec3(lhs.data[0]*rhs[0], lhs.data[1]*rhs[1], lhs.data[2]*rhs[2]);
}

// Product with another diagonal matrix
Diag3 Diag3::operator*(const Diag3& rhs) const {
    const Diag3& lhs = *this;
    return Diag3(lhs*rhs.data);
}

// Product with a full matrix
Mat3 Diag3::operator*(const Mat3& rhs) const {
    const Diag3& lhs = *this;
    return Mat3(rhs[0]*lhs.data[0], rhs[1]*lhs.data[1], rhs[2]*lhs.data[2]);
}

// Inverse:
Diag3 Diag3::inv() const {
    return Diag3(Vec3(1.0/data[0], 1.0/data[1], 1.0/data[2]));
}

// Convertion to full matrix format
Mat3 Diag3::asMat3() const {
    return Mat3(Vec3(data[0], 0.0, 0.0),
                Vec3(0.0, data[1], 0.0),
                Vec3(0.0, 0.0, data[2]));
}

/*------------------------------------------------------------------------------
 4D VECTORS (QUATERNIONS)
------------------------------------------------------------------------------*/

// Create a Quat whose elements are all 0.
Quat::Quat() {
    data[0] = data[1] = data[2] = data[3] = 0.0;
}

// Create a Quat with specified x, y, and z components.
Quat::Quat(double x, double y, double z, double w) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = w;
}

// Create unit quaternion from rotation matrix (Shepperd, 1978)
Quat::Quat(const Mat3& A) {
    double a11 = A[0][0];
    double a22 = A[1][1];
    double a33 = A[2][2];
    Quat Q2(1.0+a11+a22+a33, 1.0+a11-a22-a33, 1.0-a11+a22-a33, 1.0-a11-a22+a33);
    int imax = Q2.maxloc();
    double Q2max = Q2[imax];
    double factor = 0.5/sqrt(Q2max);
    if (imax == 0) {
        data[1] = (A[1][2]-A[2][1])*factor;
        data[2] = (A[2][0]-A[0][2])*factor;
        data[3] = (A[0][1]-A[1][0])*factor;
    }
    else if (imax == 1) {
        data[0] = (A[1][2]-A[2][1])*factor;
        data[2] = (A[0][1]+A[1][0])*factor;
        data[3] = (A[0][2]+A[2][0])*factor;
    }
    else if (imax == 2) {
        data[0] = (A[2][0]-A[0][2])*factor;
        data[1] = (A[0][1]+A[1][0])*factor;
        data[3] = (A[1][2]+A[2][1])*factor;
    }
    else {
        data[0] = (A[0][1]-A[1][0])*factor;
        data[1] = (A[0][2]+A[2][0])*factor;
        data[2] = (A[1][2]+A[2][1])*factor;
    }
    data[imax] = Q2max*factor;
}

// Extract element
double Quat::operator[](int index) const {
    assert(index >= 0 && index < 4);
    return data[index];
}

// Assign element
double& Quat::operator[](int index) {
    assert(index >= 0 && index < 4);
    return data[index];
}

// Equality
bool Quat::operator==(const Quat& rhs) const {
    return (data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2] && data[3] == rhs[3]);
}

// Inequality
bool Quat::operator!=(const Quat& rhs) const {
    return (data[0] != rhs[0] || data[1] != rhs[1] || data[2] != rhs[2] || data[3] != rhs[3]);
}

// Unary plus
Quat Quat::operator+() const {
    return Quat(*this);
}

// Plus
Quat Quat::operator+(const Quat& rhs) const {
    const Quat& lhs = *this;
    return Quat(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3]);
}

// Plus equal
Quat& Quat::operator+=(const Quat& rhs) {
    data[0] += rhs[0];
    data[1] += rhs[1];
    data[2] += rhs[2];
    data[3] += rhs[3];
    return *this;
}

// Unary minus
Quat Quat::operator-() const {
    const Quat& lhs = *this;
    return Quat(-lhs[0], -lhs[1], -lhs[2], -lhs[3]);
}

// Minus
Quat Quat::operator-(const Quat& rhs) const {
    const Quat& lhs = *this;
    return Quat(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3]);
}

// Minus equal
Quat& Quat::operator-=(const Quat& rhs) {
    data[0] -= rhs[0];
    data[1] -= rhs[1];
    data[2] -= rhs[2];
    data[3] -= rhs[3];
    return *this;
}

// Product with scalar
Quat Quat::operator*(double rhs) const {
    const Quat& lhs = *this;
    return Quat(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs, lhs[3]*rhs);
}

// Times equal
Quat& Quat::operator*=(double rhs) {
    data[0] *= rhs;
    data[1] *= rhs;
    data[2] *= rhs;
    data[3] *= rhs;
    return *this;
}

// Division with scalar
Quat Quat::operator/(double rhs) const {
    const Quat& lhs = *this;
    double scale = 1.0/rhs;
    return Quat(lhs[0]*scale, lhs[1]*scale, lhs[2]*scale, lhs[3]*scale);
}

// Divided equal
Quat& Quat::operator/=(double rhs) {
    double scale = 1.0/rhs;
    data[0] *= scale;
    data[1] *= scale;
    data[2] *= scale;
    data[3] *= scale;
    return *this;
}

// Dot product
double Quat::dot(const Quat& rhs) const {
    const Quat& lhs = *this;
    return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2] + lhs[3]*rhs[3];
}

double Quat::norm() const {
    const Quat& lhs = *this;
    return sqrt(lhs.dot(lhs));
}

// Location of maximum component
int Quat::maxloc() const {
    const Quat& lhs = *this;
    int imax = 0;
    double vmax = lhs[0];
    for (int i = 1; i < 4; i++)
        if (data[i] > vmax) {
            vmax = lhs[i];
            imax = i;
        }
    return imax;
}

// Premultiplication of B(q) by a vector
Quat Quat::B(Vec3 x) const {
    return Quat(-data[1]*x[0] - data[2]*x[1] - data[3]*x[2],
                 data[0]*x[0] - data[3]*x[1] + data[2]*x[2],
                 data[3]*x[0] + data[0]*x[1] - data[1]*x[2],
                -data[2]*x[0] + data[1]*x[1] + data[0]*x[2]);
}

// Premultiplication of C(q) by a vector
Quat Quat::C(Vec3 x) const {
    return Quat(-data[1]*x[0] - data[2]*x[1] - data[3]*x[2],
                 data[0]*x[0] + data[3]*x[1] - data[2]*x[2],
                -data[3]*x[0] + data[0]*x[1] + data[1]*x[2],
                 data[2]*x[0] - data[1]*x[1] + data[0]*x[2]);
}

// Premultiplication of B^t(q) by a quaternion
Vec3 Quat::Bt(Quat y) const {
    return Vec3(-data[1]*y[0] + data[0]*y[1] + data[3]*y[2] - data[2]*y[3],
                -data[2]*y[0] - data[3]*y[1] + data[0]*y[2] + data[1]*y[3],
                -data[3]*y[0] + data[2]*y[1] - data[1]*y[2] + data[0]*y[3]);
}

// Premultiplication of C^t(q) by a quaternion
Vec3 Quat::Ct(Quat y) const {
    return Vec3(-data[1]*y[0] + data[0]*y[1] - data[3]*y[2] + data[2]*y[3],
                -data[2]*y[0] + data[3]*y[1] + data[0]*y[2] - data[1]*y[3],
                -data[3]*y[0] - data[2]*y[1] + data[1]*y[2] + data[0]*y[3]);
}

// Premultiplication by permutation matrices:

Quat Quat::B1() const {
    const Quat& q = *this;
    return Quat(-q[1],  q[0],  q[3], -q[2]);
}

Quat Quat::B2() const {
    const Quat& q = *this;
    return Quat(-q[2], -q[3],  q[0],  q[1]);
}

Quat Quat::B3() const {
    const Quat& q = *this;
    return Quat(-q[3],  q[2], -q[1],  q[0]);
}

Vec3 Quat::A(Vec3 x) const {
    const Quat& q = *this;
    return q.Bt(q.C(x));
}

Vec3 Quat::At(Vec3 x) const {
    const Quat& q = *this;
    return q.Ct(q.B(x));
}

/*------------------------------------------------------------------------------
 MISCELLANEOUS
------------------------------------------------------------------------------*/

Projection::Projection(Vec3 x) {
    row.resize(3);
    row[0] = -x*x[0];
    row[1] = -x*x[1];
    row[2] = -x*x[2];
    double xtx = x.dot(x);
    row[0][0] += xtx;
    row[1][1] += xtx;
    row[2][2] += xtx;
}
