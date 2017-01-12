/// \file Vec3.cpp
/// \brief simple vector class for 3 floats, includes basic operations
/// \author Edina Morris

#include "Vec3.h"

vec3::vec3(float _x, float _y, float _z)
{
    m_element[0] = _x;
    m_element[1] = _y;
    m_element[2] = _z;
}

vec3::vec3(const vec3 &_value)
{
    *this = _value;
}

float vec3::dot(const vec3 &_value) const
{
    return m_x*_value.m_x+m_y*_value.m_y+m_z*_value.m_z;
}

vec3& vec3::operator +=(const vec3& _value)
{
    m_x += _value.m_x;
    m_y += _value.m_y;
    m_z += _value.m_z;
    return *this;
}

vec3& vec3::operator -=(const vec3& _value)
{
    m_x -= _value.m_x;
    m_y -= _value.m_y;
    m_z -= _value.m_z;
    return *this;
}

vec3 vec3::operator+(const vec3 _value) const
{
    return vec3(m_x+_value.m_x,
                m_y+_value.m_y,
                m_z+_value.m_z);
}

vec3 vec3::operator-(const vec3 &_value) const
{
    return vec3(m_x-_value.m_x,
                m_y-_value.m_y,
                m_z-_value.m_z);
}

vec3 vec3::operator*(const float _value) const
{
    return vec3(m_x*_value,m_y*_value,m_z*_value);
}

vec3 vec3::operator/(const float _value) const
{
    return vec3(m_x/_value,m_y/_value,m_z/_value);
}
