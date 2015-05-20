#ifndef __FLOW_C_HPP_INCLUDED__
#define __FLOW_C_HPP_INCLUDED__

#include "vector_2d_c.hpp"


namespace wif_core
{

class flow_c
{
public:
	flow_c();

	virtual ~flow_c();

	/**
	 * gives the streamfunction at a given point
	 */
	virtual double get_psi(const vector_2d_c & position) const = 0;

	/**
	 * gives the streampotential at a given point
	 */
	virtual double get_phi(const vector_2d_c & position) const = 0;

	/**
	 * gives the velocity at a given point
	 */
	virtual vector_2d_c get_velocity(const vector_2d_c & position) const = 0;

protected:

private:
};


} // namespace wif_core


#endif // __FLOW_SHEET_HPP_INCLUDED__
