#ifndef __VISUALIZATION_HPP_INCLUDED__
#define __VISUALIZATION_HPP_INCLUDED__

#include <memory>
#include <wif_core/wif_core.hpp>


namespace wif_viz
{


using wif_core::vector_2d_c;
using wif_core::flow_c;


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
	//virtual void set_velocityarrows(const vector_2d_c & bins) = 0;

	void set_contours(const std::vector<double> & contours);
	void set_contours(uint32_t contours);
	void set_clip_range(double min, double max);
	void set_output_to_file(bool file_output);
	void set_stagnation_tolerance(double epsilon);

	double clip_value(double value) const;

	void set_airfoil(wif_core::airfoil_c * new_airfoil);

	virtual void draw(const std::string & filename = "") = 0;

	virtual void draw_ivo(const std::string & filename = "")
	{
		draw(filename);
	}

	void set_streamline_seeds(const wif_core::line_2d_c & streamline_seeds);

	void set_streamline_resolution(uint32_t streamline_resolution);

	//

protected:
	std::shared_ptr<flow_c> flow;
	vector_2d_c min_range;
	vector_2d_c max_range;

	vector_2d_c psi_bins;
	vector_2d_c phi_bins;
	vector_2d_c velocity_bins;

	std::vector<double> contour_locations;
	double clip_min;
	double clip_max;
	bool output_to_file;

	mutable std::vector<vector_2d_c> stagnation_point;
	double stagnation_tolerance;

	wif_core::airfoil_c * airfoil;

	wif_core::line_2d_c streamline_seeds;
	uint32_t streamline_resolution;

private:

};


} // namespace wif_viz

#endif // __VISUALIZATION_HPP_INCLUDED__

