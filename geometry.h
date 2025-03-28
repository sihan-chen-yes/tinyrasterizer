//
// Created by csh on 2/26/25.
// Eigen like tiny library
//
#ifndef TINYRENDERER_GEOMETRY_H
#define TINYRENDERER_GEOMETRY_H
#include <cmath>
#include <vector>
#define EPSILON 1e-6

class Matrix;

template <class t> struct Vec2 {
    union {
        struct {t u, v;};
        struct {t x, y;};
        t raw[2];
    };
    Vec2() : u(0), v(0) {}
    Vec2(t _u, t _v) : u(_u),v(_v) {}
    template <class u> Vec2(const Vec2<u> &v);
    inline Vec2<t> operator +(const Vec2<t> &V) const { return Vec2<t>(u+V.u, v+V.v); }
    inline Vec2<t> operator -(const Vec2<t> &V) const { return Vec2<t>(u-V.u, v-V.v); }
    inline Vec2<t> operator *(float f)          const { return Vec2<t>(u*f, v*f); }
    t& operator[](const int i) { return raw[i]; }
    float norm () const { return std::sqrt(x*x+y*y); }
    Vec2<t> & normalize(t l=1) { *this = (*this)*(l/norm()); return *this; }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec2<t>& v);
};

template <class t> struct Vec3 {
    union {
        struct {t x, y, z;};
        struct { t ivert, iuv, inorm; };
        t raw[3];
    };
    Vec3() : x(0), y(0), z(0) {}
    Vec3(t _x, t _y, t _z) : x(_x),y(_y),z(_z) {}
    Vec3(Matrix m);  // Declare it properly
    template <class u> Vec3(const Vec3<u> &v);
    inline Vec3<t> operator ^(const Vec3<t> &v) const { return Vec3<t>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); }
    inline Vec3<t> operator +(const Vec3<t> &v) const { return Vec3<t>(x+v.x, y+v.y, z+v.z); }
    inline Vec3<t> operator -(const Vec3<t> &v) const { return Vec3<t>(x-v.x, y-v.y, z-v.z); }
    inline Vec3<t> operator *(float f)          const { return Vec3<t>(x*f, y*f, z*f); }
    inline Vec3<t> operator /(float f)          const { return Vec3<t>(x/f, y/f, z/f); }
    inline t       operator *(const Vec3<t> &v) const { return x*v.x + y*v.y + z*v.z; }
    t& operator[](const int i) { return raw[i]; }
    float norm () const { return std::sqrt(x*x+y*y+z*z); }
    Vec3<t> & normalize(t l=1) { *this = (*this)*(l/norm()); return *this; }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec3<t>& v);
};

template <class t> struct Vec4 {
    union {
        struct {t x, y, z, w;};
        t raw[4];
    };
    Vec4() : x(0), y(0), z(0), w(0) {}
    Vec4(t _x, t _y, t _z, t _w) : x(_x),y(_y),z(_z),w(_w) {}
    Vec4(Matrix m);
    template <class u> Vec4(const Vec4<u> &v);
    inline Vec4<t> operator +(const Vec4<t> &v) const { return Vec3<t>(x+v.x, y+v.y, z+v.z, w+v.w); }
    inline Vec4<t> operator -(const Vec4<t> &v) const { return Vec3<t>(x-v.x, y-v.y, z-v.z, w-v.w); }
    inline Vec4<t> operator *(float f)          const { return Vec3<t>(x*f, y*f, z*f, w*f); }
    inline Vec4<t> operator /(float f)          const { return Vec3<t>(x/f, y/f, z/f, w/f); }
    inline t       operator *(const Vec4<t> &v) const { return x*v.x + y*v.y + z*v.z + w*v.w; }
    t& operator[](const int i) { return raw[i]; }
    float norm () const { return std::sqrt(x*x+y*y+z*z+w*w); }
    Vec4<t> & normalize(t l=1) { *this = (*this)*(l/norm()); return *this; }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec4<t>& v);
};

typedef Vec2<float> Vec2f;
typedef Vec2<int>   Vec2i;
typedef Vec3<float> Vec3f;
typedef Vec3<int>   Vec3i;
typedef Vec4<float>   Vec4f;
typedef Vec4<int>   Vec4i;

template <> template <> Vec2<int>::Vec2(const Vec2<float> &v);
template <> template <> Vec2<float>::Vec2(const Vec2<int> &v);
template <> template <> Vec3<int>::Vec3(const Vec3<float> &v);
template <> template <> Vec3<float>::Vec3(const Vec3<int> &v);
template <> template <> Vec4<float>::Vec4(const Vec4<int> &v);
template <> template <> Vec4<int>::Vec4(const Vec4<float> &v);

template <typename T>
Vec3<T> to_cartesian(const Vec4<T> &v);

template <> Vec3<float> to_cartesian(const Vec4<float> &v);

template <typename T>
Vec4<T> to_homogeneous(const Vec3<T> &v);

template <> Vec4<float> to_homogeneous(const Vec3<float> &v);

template <class t> inline Vec2<t> operator*(float f, const Vec2<t>& v) {
    return Vec2<t>(v.x * f, v.y * f);
}

template <class t> inline Vec3<t> operator*(float f, const Vec3<t>& v) {
    return Vec3<t>(v.x * f, v.y * f, v.z * f);
}

template <class t> inline Vec4<t> operator*(float f, const Vec4<t>& v) {
    return Vec4<t>(v.x * f, v.y * f, v.z * f, v.w * f);
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec2<t>& v) {
    s << "(" << v.x << ", " << v.y << ")\n";
    return s;
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec3<t>& v) {
    s << "(" << v.x << ", " << v.y << ", " << v.z << ")\n";
    return s;
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec4<t>& v) {
    s << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")\n";
    return s;
}

const int DEFAULT_ALLOC=4;

class Matrix {
    std::vector<std::vector<float>> m;
    int rows, cols;
public:
    Matrix(int r=DEFAULT_ALLOC, int c=DEFAULT_ALLOC);
    template <class t> Matrix(Vec3<t> v);
    template <class t> Matrix(Vec4<t> v);
    inline int nrows();
    inline int ncols();

    static Matrix identity(int dimensions);
    std::vector<float>& operator[](const int i);
    Matrix operator*(const Matrix& a);
    Matrix transpose();
    Matrix inverse();

    void set_row(int i, Vec3f r);
    void set_col(int j, Vec3f c);

    friend std::ostream& operator<<(std::ostream& s, Matrix& m);
};

#endif //TINYRENDERER_GEOMETRY_H
