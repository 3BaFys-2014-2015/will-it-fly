#include "vector_2d_c.hpp"


namespace wif_core
{


vector_2d_c::vector_2d_c(double x, double y) :
	x(x),
	y(y)
{
	//
}


vector_2d_c::vector_2d_c(const vector_2d_c & other) :
	x(other.x),
	y(other.y)
{
	//
}


vector_2d_c::vector_2d_c() :
	x(0.0),
	y(0.0)
{
	//
}


bool vector_2d_c::operator==(const vector_2d_c & other) const
{
	return (*this - other).get_length_sq() < 0.000001;
}


vector_2d_c & vector_2d_c::operator=(const vector_2d_c & other)
{
	this->x = other.x;
	this->y = other.y;
	return *this;
}


vector_2d_c vector_2d_c::operator+(const vector_2d_c & other) const
{
	return vector_2d_c(this->x + other.x, this->y + other.y);
}


vector_2d_c vector_2d_c::operator-() const
{
	return vector_2d_c(-this->x, -this->y);
}


vector_2d_c vector_2d_c::operator-(const vector_2d_c & other) const
{
	return vector_2d_c(this->x - other.x, this->y - other.y);
}


vector_2d_c vector_2d_c::operator*(double r) const
{
	return vector_2d_c(r * this->x, r * this->y);
}


vector_2d_c vector_2d_c::operator/(double r) const
{
	return vector_2d_c(this->x / r, this->y / r);
}


vector_2d_c & vector_2d_c::operator+=(const vector_2d_c & other)
{
	this->x += other.x;
	this->y += other.y;
	return * this;
}


vector_2d_c & vector_2d_c::operator-=(const vector_2d_c & other)
{
	this->x -= other.x;
	this->y -= other.y;
	return * this;
}


vector_2d_c & vector_2d_c::operator*=(double r)
{
	this->x *= r;
	this->y *= r;
	return * this;
}


vector_2d_c & vector_2d_c::operator/=(double r)
{
	this->x /= r;
	this->y /= r;
	return * this;
}


double vector_2d_c::get_length() const
{
	return std::sqrt(get_length_sq());
}


double vector_2d_c::get_length_sq() const
{
	return (this->x * this->x + this->y * this->y);
}


double vector_2d_c::get_angle() const
{
	return std::atan2(this->y, this->x);
}


double vector_2d_c::dot(const vector_2d_c & other) const
{
	return (this->x * other.x + this->y * other.y);
}


double vector_2d_c::cross(const vector_2d_c & other) const
{
	return (this->x * other.y - this->y * other.x);
}


vector_2d_c operator * (double n, const vector_2d_c & vec)
{
	return vector_2d_c(n * vec.x, n * vec.y);
}


/**
 * creates a new vector_2d_c object, but with the length and the angle with the x-axis as parameters
 */
vector_2d_c vector_2d_radian(double radian, double angle)
{
	return vector_2d_c(radian * std::cos(angle), radian * std::sin(angle));
}


std::ostream & operator <<(std::ostream & output, const vector_2d_c & vec)
{
	output << "(" << vec.x << ", " << vec.y << ")";
	return output;
}


} // namespace wif_core
