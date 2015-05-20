#ifndef __VECTOR_2D_C_HPP_INCLUDED__
#define __VECTOR_2D_C_HPP_INCLUDED__

#include <cmath>
#include <ostream>


namespace wif_core
{


class vector_2d_c
{
public:
	vector_2d_c();

	vector_2d_c(double x, double y);

	vector_2d_c(const vector_2d_c & other);

	bool operator==(const vector_2d_c & other) const;

	vector_2d_c & operator=(const vector_2d_c & other);

	vector_2d_c operator+(const vector_2d_c & other) const;

	vector_2d_c operator-() const;

	vector_2d_c operator-(const vector_2d_c & other) const;

	vector_2d_c operator*(double r) const;

	vector_2d_c operator/(double r) const;

	vector_2d_c & operator+=(const vector_2d_c & other);

	vector_2d_c & operator-=(const vector_2d_c & other);

	vector_2d_c & operator*=(double r);

	vector_2d_c & operator/=(double r);

	/**
	 * gives the norm of the vector
	 */
	double get_length() const;

	/**
	 * same as pow(get_length(),2)
	 */
	double get_length_sq() const;

	/**
	 * get angle with x-axis
	 */
	double get_angle() const;

	/**
	 * x1y1+x2y2
	 * get dot product
	 */
	double dot(const vector_2d_c & other) const;

	/**
	 * get x1y2-x2y1
	 * cross product in 2 dimensions
	 * same as norm of cross product in 3 dimensions
	 */
	double cross(const vector_2d_c & other) const;

public:
	double x;
	double y;
};


vector_2d_c operator *(double, const vector_2d_c &);


/**
 * creates a new vector_2d_c object, but with the length and the angle with the x-axis as parameters
 */
vector_2d_c vector_2d_radian(double radius, double angle);


/**
 * sends vector to ostream
 */
std::ostream & operator << (std::ostream &, const vector_2d_c &);


} // namespace wif_core

#endif // __VECTOR_2D_C_HPP_INCLUDED__

