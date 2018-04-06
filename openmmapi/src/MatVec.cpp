/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "internal/MatVec.h"

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
    return Vec3(row[1][1], row[2][2], row[3][3]);
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

// Transposition
Mat3 Mat3::t() const {
    return Mat3(this->col(0), this->col(1), this->col(2));
}

// Product with vector
Vec3 Mat3::operator*(const Vec3& rhs) const {
    const Mat3& lhs = *this;
    return Vec3(lhs[0].dot(rhs), lhs[1].dot(rhs), lhs[2].dot(rhs));
}

// Product with another full matrix
Mat3 Mat3::operator*(const Mat3& rhs) const {
    const Mat3& lhs = *this;
    return Mat3(lhs*rhs.col(0), lhs*rhs.col(1), lhs*rhs*col(2)).t();
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
 4D VECTORS
------------------------------------------------------------------------------*/

/**
 * Create a Vec4 whose elements are all 0.
 */
Vec4::Vec4() {
    data[0] = data[1] = data[2] = data[3] = 0.0;
}
/**
 * Create a Vec4 with specified x, y, and z components.
 */
Vec4::Vec4(double x, double y, double z, double w) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = w;
}
double Vec4::operator[](int index) const {
    assert(index >= 0 && index < 4);
    return data[index];
}
double& Vec4::operator[](int index) {
    assert(index >= 0 && index < 4);
    return data[index];
}
bool Vec4::operator==(const Vec4& rhs) const {
    return (data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2] && data[3] == rhs[3]);
}
bool Vec4::operator!=(const Vec4& rhs) const {
    return (data[0] != rhs[0] || data[1] != rhs[1] || data[2] != rhs[2] || data[3] != rhs[3]);
}
// Arithmetic operators
// unary plus
Vec4 Vec4::operator+() const {
    return Vec4(*this);
}
// plus
Vec4 Vec4::operator+(const Vec4& rhs) const {
    const Vec4& lhs = *this;
    return Vec4(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3]);
}

Vec4& Vec4::operator+=(const Vec4& rhs) {
    data[0] += rhs[0];
    data[1] += rhs[1];
    data[2] += rhs[2];
    data[3] += rhs[3];
    return *this;
}
// unary minus
Vec4 Vec4::operator-() const {
    const Vec4& lhs = *this;
    return Vec4(-lhs[0], -lhs[1], -lhs[2], -lhs[3]);
}

// minus
Vec4 Vec4::operator-(const Vec4& rhs) const {
    const Vec4& lhs = *this;
    return Vec4(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3]);
}
Vec4& Vec4::operator-=(const Vec4& rhs) {
    data[0] -= rhs[0];
    data[1] -= rhs[1];
    data[2] -= rhs[2];
    data[3] -= rhs[3];
    return *this;
}
// scalar product
Vec4 Vec4::operator*(double rhs) const {
    const Vec4& lhs = *this;
    return Vec4(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs, lhs[3]*rhs);
}
Vec4& Vec4::operator*=(double rhs) {
    data[0] *= rhs;
    data[1] *= rhs;
    data[2] *= rhs;
    data[3] *= rhs;
    return *this;
}
// scalar division
Vec4 Vec4::operator/(double rhs) const {
    const Vec4& lhs = *this;
    double scale = 1.0/rhs;
    return Vec4(lhs[0]*scale, lhs[1]*scale, lhs[2]*scale, lhs[3]*scale);
}
Vec4& Vec4::operator/=(double rhs) {
    double scale = 1.0/rhs;
    data[0] *= scale;
    data[1] *= scale;
    data[2] *= scale;
    data[3] *= scale;
    return *this;
}
// dot product
double Vec4::dot(const Vec4& rhs) const {
    const Vec4& lhs = *this;
    return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2] + lhs[3]*rhs[3];
}
// location of maximum component
int Vec4::maxloc() {
    const Vec4& lhs = *this;
    int imax = 0;
    double vmax = lhs[0];
    for (int i = 1; i < 4; i++)
        if (data[i] > vmax) {
            vmax = lhs[i];
            imax = i;
        }
    return imax;
}

/*------------------------------------------------------------------------------
 MISCELLANEOUS
------------------------------------------------------------------------------*/

Mat3 RankOne(Vec3 x, Vec3 y) {
    return Diag3(x)*Mat3(y,y,y);
}
