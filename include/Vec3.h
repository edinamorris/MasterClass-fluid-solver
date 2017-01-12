/// \file Vec3.h
/// \brief simple vector class for 3 floats, includes basic operations
/// \author Edina Morris

#ifndef VEC3_H_
#define VEC3_H_

#include <iostream>
#include <cmath>

class vec3 {

public:

    vec3(float _x = 0, float _y = 0, float _z = 0);

    vec3(const vec3& _value);
    vec3& operator=(float _s) { m_element[0] = m_element[1] = m_element[2] = _s; return *this; }

    // Access methods
    operator const float*() const { return m_element; }
    float& operator[](int _i) { return m_element[_i]; }

    vec3& operator+=(const vec3& _value);
    vec3& operator-=(const vec3& _value);
    vec3 operator+(const vec3 _value) const;
    vec3 operator-(const vec3 &_value) const;
    vec3 operator*(const float _value) const;
    vec3 operator/(const float _value) const;

    //dot product
    float dot(const vec3 &_value) const;

    // the data
    union
    {
       struct { float m_x,m_y,m_z; };
       struct { float m_r,m_g,m_b; };
       float m_element[3];
    };

};

inline vec3 operator*(const float a, const vec3& b)
{ return vec3(b.m_x * a, b.m_y * a, b.m_z * a); }

#endif // VEC3_H_

