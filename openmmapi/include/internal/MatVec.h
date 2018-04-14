#ifndef OPENMM_MATVEC_H_
#define OPENMM_MATVEC_H_

/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include <cassert>
#include <iosfwd>
#include "openmm/Vec3.h"
#include <vector>

using namespace OpenMM;
using std::vector;

namespace RigidBodyPlugin {

// Forward declaration:
class Diag3;

/**
 * This class represents a 3x3 matrix
 */

class Mat3 {
public:
    Mat3();
    Mat3(Vec3 x, Vec3 y, Vec3 z);

    Vec3 operator[](int index) const;
    Vec3& operator[](int index);
    Vec3 col(int index) const;
    Vec3 diag() const;
    bool operator==(const Mat3& rhs) const;
    bool operator!=(const Mat3& rhs) const;
    Mat3 operator+() const;
    Mat3 operator+(const Mat3& rhs) const;
    Mat3& operator+=(const Mat3& rhs);
    Mat3 operator-() const;
    Mat3 operator-(const Mat3& rhs) const;
    Mat3 operator-(const Diag3& rhs) const;
    Mat3& operator-=(const Mat3& rhs);
    Mat3 operator*(double rhs) const;
    Mat3& operator*=(double rhs);
    Mat3 operator/(double rhs) const;
    Mat3& operator/=(double rhs);
    Mat3 symmetric() const;
    double trace() const;
    double det() const;
    Mat3 t() const;
    Vec3 operator*(const Vec3& rhs) const;
    Mat3 operator*(const Mat3& rhs) const;
    Mat3 operator*(const Diag3& rhs) const;
protected:
    vector<Vec3> row;
};

class Diag3 {
public:
    Diag3();
    Diag3(Vec3 x);
    Diag3(double x);

    Vec3 operator[](int index) const;
    Vec3 diag() const;
    bool operator==(const Diag3& rhs) const;
    bool operator!=(const Diag3& rhs) const;
    Diag3 operator+() const;
    Diag3 operator+(const Diag3& rhs) const;
    Diag3& operator+=(const Diag3& rhs);
    Diag3 operator-() const;
    Diag3 operator-(const Diag3& rhs) const;
    Diag3& operator-=(const Diag3& rhs);
    Diag3 operator*(double rhs) const;
    Diag3& operator*=(double rhs);
    Diag3 operator/(double rhs) const;
    Diag3& operator/=(double rhs);
    Diag3 t() const;
    Vec3 operator*(const Vec3& rhs) const;
    Diag3 operator*(const Diag3& rhs) const;
    Mat3 operator*(const Mat3& rhs) const;
    Diag3 inv() const;
    Mat3 asMat3() const;
private:
    Vec3 data;
};

/**
 * This class represents a four component vector.  It is used for storing quaternions.
 */

class Quat {
public:
    Quat();
    Quat(double x, double y, double z, double w);
    Quat(const Mat3& A);

    double operator[](int index) const;
    double& operator[](int index);
    bool operator==(const Quat& rhs) const;
    bool operator!=(const Quat& rhs) const;
    Quat operator+() const;
    Quat operator+(const Quat& rhs) const;
    Quat& operator+=(const Quat& rhs);
    Quat operator-() const;
    Quat operator-(const Quat& rhs) const;
    Quat& operator-=(const Quat& rhs);
    Quat operator*(double rhs) const;
    Quat& operator*=(double rhs);
    Quat operator/(double rhs) const;
    Quat& operator/=(double rhs);
    double dot(const Quat& rhs) const;
    double norm() const;
    int maxloc() const;
    Quat B(Vec3 x) const;
    Quat C(Vec3 x) const;
    Vec3 Bt(Quat y) const;
    Vec3 Ct(Quat y) const;
    Vec3 A(Vec3 x) const;
    Vec3 At(Vec3 x) const;
private:
    double data[4];
};

class Projection : public Mat3 {
public:
    Projection(Vec3 x);
};

template <class CHAR, class TRAITS>
std::basic_ostream<CHAR,TRAITS>& operator<<(std::basic_ostream<CHAR,TRAITS>& o, const Mat3& M) {
    o<<'['<<M[0]<<"; "<<M[1]<<"; "<<M[2]<<']';
    return o;
}

template <class CHAR, class TRAITS>
std::basic_ostream<CHAR,TRAITS>& operator<<(std::basic_ostream<CHAR,TRAITS>& o, const Diag3& D) {
    o<<'['<<D[0]<<"; "<<D[1]<<"; "<<D[2]<<']';
    return o;
}

template <class CHAR, class TRAITS>
std::basic_ostream<CHAR,TRAITS>& operator<<(std::basic_ostream<CHAR,TRAITS>& o, const Quat& q) {
    o<<'['<<q[0]<<", "<<q[1]<<", "<<q[2]<<", "<<q[3]<<']';
    return o;
}

} // namespace RigidBodyPlugin

#endif /*OPENMM_MATVEC_H_*/
