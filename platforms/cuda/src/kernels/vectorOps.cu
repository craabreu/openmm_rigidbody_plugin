/**
 * This file defines vector operations to simplify code elsewhere.
 */

// Versions of make_x() that take a single value and set all components to that.

inline __device__ float3 make_float3(float a) {
    return make_float3(a, a, a);
}

inline __device__ float4 make_float4(float a) {
    return make_float4(a, a, a, a);
}

inline __device__ double3 make_double3(double a) {
    return make_double3(a, a, a);
}

inline __device__ double4 make_double4(double a) {
    return make_double4(a, a, a, a);
}

// Negate a vector.

inline __device__ float3 operator-(float3 a) {
    return make_float3(-a.x, -a.y, -a.z);
}

inline __device__ float4 operator-(float4 a) {
    return make_float4(-a.x, -a.y, -a.z, -a.w);
}

inline __device__ double3 operator-(double3 a) {
    return make_double3(-a.x, -a.y, -a.z);
}

inline __device__ double4 operator-(double4 a) {
    return make_double4(-a.x, -a.y, -a.z, -a.w);
}

// Add two vectors.

inline __device__ float3 operator+(float3 a, float3 b) {
    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

inline __device__ float4 operator+(float4 a, float4 b) {
    return make_float4(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}

inline __device__ double3 operator+(double3 a, double3 b) {
    return make_double3(a.x+b.x, a.y+b.y, a.z+b.z);
}

inline __device__ double4 operator+(double4 a, double4 b) {
    return make_double4(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}

// Subtract two vectors.

inline __device__ float3 operator-(float3 a, float3 b) {
    return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
}

inline __device__ float4 operator-(float4 a, float4 b) {
    return make_float4(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}

inline __device__ double3 operator-(double3 a, double3 b) {
    return make_double3(a.x-b.x, a.y-b.y, a.z-b.z);
}

inline __device__ double4 operator-(double4 a, double4 b) {
    return make_double4(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}

// Multiply two vectors.

inline __device__ float3 operator*(float3 a, float3 b) {
    return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);
}

inline __device__ float4 operator*(float4 a, float4 b) {
    return make_float4(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}

inline __device__ double3 operator*(double3 a, double3 b) {
    return make_double3(a.x*b.x, a.y*b.y, a.z*b.z);
}

inline __device__ double4 operator*(double4 a, double4 b) {
    return make_double4(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}

// Divide two vectors.

inline __device__ float3 operator/(float3 a, float3 b) {
    return make_float3(a.x/b.x, a.y/b.y, a.z/b.z);
}

inline __device__ float4 operator/(float4 a, float4 b) {
    return make_float4(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}

inline __device__ double3 operator/(double3 a, double3 b) {
    return make_double3(a.x/b.x, a.y/b.y, a.z/b.z);
}

inline __device__ double4 operator/(double4 a, double4 b) {
    return make_double4(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}

// += operator

inline __device__ void operator+=(float3& a, float3 b) {
    a.x += b.x; a.y += b.y; a.z += b.z;
}

inline __device__ void operator+=(float4& a, float4 b) {
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}

inline __device__ void operator+=(double3& a, double3 b) {
    a.x += b.x; a.y += b.y; a.z += b.z;
}

inline __device__ void operator+=(double4& a, double4 b) {
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}

// -= operator

inline __device__ void operator-=(float3& a, float3 b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

inline __device__ void operator-=(float4& a, float4 b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}

inline __device__ void operator-=(double3& a, double3 b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

inline __device__ void operator-=(double4& a, double4 b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}

// *= operator

inline __device__ void operator*=(float3& a, float3 b) {
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}

inline __device__ void operator*=(float4& a, float4 b) {
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}

inline __device__ void operator*=(double3& a, double3 b) {
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
}

inline __device__ void operator*=(double4& a, double4 b) {
    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;
}

// /= operator

inline __device__ void operator/=(float3& a, float3 b) {
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
}

inline __device__ void operator/=(float4& a, float4 b) {
    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;
}

inline __device__ void operator/=(double3& a, double3 b) {
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
}

inline __device__ void operator/=(double4& a, double4 b) {
    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;
}

// Multiply a vector by a constant.

inline __device__ float3 operator*(float3 a, float b) {
    return make_float3(a.x*b, a.y*b, a.z*b);
}

inline __device__ float4 operator*(float4 a, float b) {
    return make_float4(a.x*b, a.y*b, a.z*b, a.w*b);
}

inline __device__ double3 operator*(double3 a, double b) {
    return make_double3(a.x*b, a.y*b, a.z*b);
}

inline __device__ double4 operator*(double4 a, double b) {
    return make_double4(a.x*b, a.y*b, a.z*b, a.w*b);
}

// Divide a vector by a constant.

inline __device__ float3 operator/(float3 a, float b) {
    float scale = 1.0f/b;
    return a*scale;
}

inline __device__ float4 operator/(float4 a, float b) {
    float scale = 1.0f/b;
    return a*scale;
}

inline __device__ double3 operator/(double3 a, double b) {
    double scale = 1.0/b;
    return a*scale;
}

inline __device__ double4 operator/(double4 a, double b) {
    double scale = 1.0/b;
    return a*scale;
}

// *= operator (multiply vector by constant)

inline __device__ void operator*=(float3& a, float b) {
    a.x *= b; a.y *= b; a.z *= b;
}

inline __device__ void operator*=(float4& a, float b) {
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

inline __device__ void operator*=(double3& a, double b) {
    a.x *= b; a.y *= b; a.z *= b;
}

inline __device__ void operator*=(double4& a, double b) {
    a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

// Dot product

inline __device__ float dot(float3 a, float3 b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

inline __device__ double dot(double3 a, double3 b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

inline __device__ float dot(float4 a, float4 b) {
    return a.x*b.x+a.y*b.y+a.z*b.z+a.w*b.w;
}

inline __device__ double dot(double4 a, double4 b) {
    return a.x*b.x+a.y*b.y+a.z*b.z+a.w*b.w;
}

// Cross product

inline __device__ float3 cross(float3 a, float3 b) {
    return make_float3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

inline __device__ double3 cross(double3 a, double3 b) {
    return make_double3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

// Normalize a vector

inline __device__ float3 normalize(float3 a) {
    return a*rsqrtf(a.x*a.x+a.y*a.y+a.z*a.z);
}

inline __device__ float4 normalize(float4 a) {
    return a*rsqrtf(a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w);
}

inline __device__ double3 normalize(double3 a) {
    return a*rsqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

inline __device__ double4 normalize(double4 a) {
    return a*rsqrt(a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w);
}

// Strip off the fourth component of a vector.

inline __device__ float3 trim(float4 v) {
    return make_float3(v.x, v.y, v.z);
}

inline __device__ double3 trim(double4 v) {
    return make_double3(v.x, v.y, v.z);
}

// Add fourth component to a vector.

inline __device__ float4 fuse(float3 v, float a) {
    return make_float4(v.x, v.y, v.z, a);
}

inline __device__ double4 fuse(double3 v, double a) {
    return make_double4(v.x, v.y, v.z, a);
}
