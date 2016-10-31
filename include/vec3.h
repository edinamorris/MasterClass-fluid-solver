#ifndef VEC3
#define VEC3

#include <iostream>
#include <cmath>

class vec3 {

public:

  vec3(float x = 0, float y = 0, float z = 0) { _element[0] = x; _element[1] = y; _element[2] = z; }

  vec3(const vec3& v) { *this = v; }
  vec3& operator=(float s) { _element[0] = _element[1] = _element[2] = s; return *this; }

  // Access methods
  operator const float*() const { return _element; }
  float& operator[](int i) { return _element[i]; }

  vec3& operator+=(const vec3& v)     { x += v.x; y += v.y; z += v.z; return *this; }
  vec3& operator-=(const vec3& v)     { x -= v.x; y -= v.y; z -= v.z; return *this; }
  vec3& operator/=(const vec3& v)     { x /= v.x; y /= v.y; z /= v.z; return *this; }
  vec3& operator/=(const float b)     { x /= b; y /= b; z /= b; return *this; }
  vec3& operator*=(const float b)     { x *= b; y *= b; z *= b;   return *this; }
  vec3 operator+(const vec3 b) const { return vec3(x+b.x,y+b.y,z+b.z); }
  vec3 operator-(const vec3 &b) const { return vec3(x-b.x,y-b.y,z-b.z); }
  vec3 operator*(const float b) const { return vec3(x*b,y*b,z*b); }
  vec3 operator/(const float b) const { return vec3(x/b,y/b,z/b); }

  // this computes the cross product
  vec3 operator^(const vec3 &v) const { return vec3(y*v.z-v.y*z,-x*v.z+v.x*z,x*v.y-v.x*y); }

  // *this one* does the dot product
  float dot(const vec3 &b) const      { return x*b.x+y*b.y+z*b.z; }

  // the data
  union {
     struct { float x,y,z; };
     struct { float r,g,b; };
     float _element[3];
  };
};

inline std::istream &operator>>(std::istream &in, vec3& v)
{ return in >> v[0] >> v[1] >> v[2]; }

inline std::ostream &operator<<(std::ostream &out, vec3& v)
{ return out << v[0] << " " << v[1] << " " << v[2]; }

inline vec3 operator*(const float a, const vec3& b)
{ return vec3(b.x * a, b.y * a, b.z * a); }

#endif // VEC3

