#ifndef __VISUALIZATION_HPP_INCLUDED__
#define __VISUALIZATION_HPP_INCLUDED__

#include <memory>
#include <wif_core/wif_core.hpp>
#include <iostream>

namespace wif_viz
{

using namespace wif_core;

enum E_SCALAR_DRAW_STYLE
{
	ESDS_GRADIENT = 0x1,
	ESDS_DISCRETE = 0x2,
	ESDS_CONTOURS = 0x4
};

enum E_VECTOR_DRAW_STYLE
{
	EVDS_ARROWS      = 0x1,
	EVDS_STREAMLINES = 0x2
};

class global_settings_c
{
public:
	global_settings_c() :
		output_to_file(false),
		draw_scale(true),
		autoguess_stagnation_points(false),
		stagnation_tolerance(0.001)
	{
		//
	}

	bool output_to_file;
	bool draw_scale;
	bool autoguess_stagnation_points;
	double stagnation_tolerance;
};

class field_c
{
public:
	field_c() :
		bins(0, 0)
	{
		//
	}

	vector_2d_c bins;
};

class scalar_field_c : public field_c
{
public:
	scalar_field_c() :
		gradient_colors(10),
		contour_locations(),
		style(ESDS_GRADIENT)
	{
		//
	}

	uint32_t gradient_colors;
	std::vector<double> contour_locations;
	E_SCALAR_DRAW_STYLE style;
};

class vector_field_c : public field_c
{
public:
	vector_field_c() :
		arrow_scale(0.1),
		streamline_seeds(0, 0, 0, 0),
		streamline_resolution(100),
		style(EVDS_ARROWS)
	{
		//
	}

	double arrow_scale;
	line_2d_c streamline_seeds;
	uint32_t streamline_resolution;
	E_VECTOR_DRAW_STYLE style;
};

class visualization_c
{
public:

	visualization_c(std::shared_ptr<flow_c> flow, const vector_2d_c & min_range, const vector_2d_c & max_range);

	virtual ~visualization_c();

	/**
	 * Redelijk voor de hand liggend, de minima van het plot-bereik zitten in min,
	 * de maxima in max.
	 */
	void set_range(const vector_2d_c & new_min_range, const vector_2d_c & new_max_range);
	void default_contour_locations();

	/**
	 * Standaard zijn de bins == (0, 0)
	 * Als de bins == (0, 0), moet niks getekend worden van dat veld.
	 * Als de bins != (0, 0), worden de bins gevonden door round(abs(x)), round(abs(y))
	 * Ge rond de absolute waarden van de vector af.
	 */
	void set_psi_bins(const vector_2d_c & bins);
	void set_phi_bins(const vector_2d_c & bins);
	void set_velocity_bins(const vector_2d_c & bins);

	/**
	 * Tekent de velden/airfoil/extra stroomlijnen/stagnatiepunten.
	 *
	 * Pas in deze method worden de grids opgevuld met punten.
	 *
	 * Als de filename == "", print naar het scherm, anders naar het
	 * bestand dat gegeven is door filename.c_str()
	 */

	void set_contours(const std::vector<double> & contours);
	void set_contours(uint32_t contours);
	void set_output_to_file(bool file_output);
	void set_stagnation_tolerance(double epsilon);

	void set_airfoil(wif_core::airfoil_c * new_airfoil);

	virtual void draw(const std::string & name = "") = 0;

	virtual void plotVectors(std::vector<std::vector<double>>, std::vector<double>, std::vector<std::string>, std::string, std::string, std::string, std::string);

	//

	scalar_field_c & get_psi_field()
	{
		return this->psi_field;
	}

	scalar_field_c & get_phi_field()
	{
		return this->phi_field;
	}

	vector_field_c & get_v_field()
	{
		return this->v_field;
	}

	global_settings_c & get_global_settings()
	{
		return this->global_settings;
	}

protected:
	std::shared_ptr<flow_c> flow;
	vector_2d_c min_range;
	vector_2d_c max_range;

	vector_2d_c psi_bins;
	vector_2d_c phi_bins;
	vector_2d_c velocity_bins;

	std::vector<double> contour_locations;
	bool output_to_file;

	double arrow_scale;

	mutable std::vector<vector_2d_c> stagnation_point;
	double stagnation_tolerance;

	wif_core::airfoil_c * airfoil;

	wif_core::line_2d_c streamline_seeds;
	uint32_t streamline_resolution;

protected:
	scalar_field_c psi_field;
	scalar_field_c phi_field;
	vector_field_c v_field;

	global_settings_c global_settings;
};


} // namespace wif_viz

#endif // __VISUALIZATION_HPP_INCLUDED__

