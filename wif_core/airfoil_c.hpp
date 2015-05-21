#ifndef __AIRFOIL_C_HPP_INCLUDED__
#define __AIRFOIL_C_HPP_INCLUDED__


#include <string>
#include <vector>

#include "line_2d_c.hpp"


namespace wif_core
{


class airfoil_c
{
public:
	/**
	 * Empty constructor. Gives an airfoil without points nor lines.
	 */
	airfoil_c();

	/**
	 * Generates a polygon with gives center, raduis and corners. Use to simulate a circles for large amount of corners
	 * @midpoint: center point of circle/polygon
	 * @radius:: radius of center/polygon
	 * @corners: amount of corners
	 */
	airfoil_c(const vector_2d_c & midpoint, double radius, unsigned int corners);

	/**
	 * Read the file's data
	 */
	airfoil_c(const std::string & filename);

	/**
	 * Create airfoil with given name and vector
	 */
	airfoil_c(const std::vector<vector_2d_c> & points, const std::string & name);

	/**
	 * New airfoil named old name + " circle projected with $n subdivisions centered on $projection_center with radius $radius".
	 * @param n subdivisions
	 * @param projection_center centerpoint of the circle that is projected from
	 * @param radius raduis of the circle taht is prjected from
	 * @param epsilon see intersect
	 */
	airfoil_c get_circle_projection(uint32_t n, const vector_2d_c & projection_center, double radius, double epsilon = 0.0001) const;

	airfoil_c get_circle_projection(uint32_t n, double epsilon = 0.0001) const;

	bool check_lengths() const;

	typename std::vector<vector_2d_c>::size_type get_index_of_last_upper_panel() const;

	typename std::vector<vector_2d_c>::size_type get_index_of_first_lower_panel() const;

	std::vector<double> get_upper_panels_x() const;

	std::vector<double> get_lower_panels_x() const;

	/**
	 * Check if the last points is the same as the first, within a certain error
	 * @epsilon maximum distance as with points are considered the same
	 */
	bool is_closed(double epsilon = 0.0001) const;

	/**
	* returns a closed airfoil by adding a new point inbetween the start and endpoint
	*/
	airfoil_c closed_merge(double epsilon = 0.0001) const;

	/**
	* returns a closed airfoil  by adding a new point on the intersection of the first and last frame
	*/
	airfoil_c closed_intersect(double epsilon = 0.0001) const;

	/**
	 * Check if the airfoil contains at least 1 point.
	 */
	bool is_valid() const;

	/**
	 * Returns the points
	 */
	const std::vector<vector_2d_c> & get_points() const;

	/**
	 * Makes pairss ([i], [i+1]) lines, with i = 0 to points.size()-1.
	 */
	std::vector<line_2d_c> get_lines() const;

	/**
	 * Same as get_lines() but in the reversed order
	 */
	std::vector<line_2d_c> get_lines_reversed() const;

	/**
	 * Gives the name
	 */
	std::string get_name() const;

	void set_name(const std::string & new_name);

private:
	/**
	 * returns first intersection with line
	 */
	vector_2d_c get_intersection_first(const line_2d_c & line) const;

	/**
	 * returns last intersection with line
	 */
	vector_2d_c get_intersection_last(const line_2d_c & line) const;

private:
	std::string name;
	std::vector<vector_2d_c> points;

	friend std::ostream & operator << (std::ostream & output, const airfoil_c & airfoil);

};


std::ostream & operator << (std::ostream & output, const airfoil_c & airfoil);


} // namespace wif_core

#endif // __AIRFOIL_C_HPP_INCLUDED__
