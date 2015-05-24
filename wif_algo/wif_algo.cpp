#include "wif_algo.hpp"
#include <gsl/gsl_linalg.h>
#include <algorithm>
#include <gsl/gsl_integration.h>

namespace wif_algo
{


uint32_t get_version()
{
	wif_core::get_version();
	return 1;
}

struct integration_function_parameters
{
	double beta;
	double betaj;
	double xc;
	double yc;
	double xa;
	double ya;
};

double xj(double xc, double xa, double s, double beta)
{
	return xc - (xa - s * std::sin(beta));
}

double yj(double yc, double ya, double s, double beta)
{
	return yc - (ya + s * std::cos(beta));
}

double source_sheet_function(double s, void * parameters)
{
	struct integration_function_parameters * params = (struct integration_function_parameters *)parameters;
	double beta = (params->beta);
	double betaj = (params->betaj);
	double xc = (params->xc);
	double yc = (params->yc);
	double xa = (params->xa);
	double ya = (params->ya);

	double a = cos(beta) * xj(xc, xa, s, betaj) + sin(beta) * yj(yc, ya, s, betaj);
	double b = pow(xj(xc, xa, s, betaj), 2.) + pow(yj(yc, ya, s, betaj), 2.);

	return a / b;
}

double vortex_sheet_function_1(double s, void * parameters)
{
	struct integration_function_parameters * params = (struct integration_function_parameters *)parameters;
	double beta = (params->beta);
	double betaj = (params->betaj);
	double xc = (params->xc);
	double yc = (params->yc);
	double xa = (params->xa);
	double ya = (params->ya);

	double a = sin(beta) * xj(xc, xa, s, betaj) - cos(beta) * yj(yc, ya, s, betaj);
	double b = pow(xj(xc, xa, s, betaj), 2) + pow(yj(yc, ya, s, betaj), 2);

	return a / b;
}

double vortex_sheet_function_lastrow(double s, void * parameters)
{
	struct integration_function_parameters * params = (struct integration_function_parameters *)parameters;
	double beta = (params->beta);
	double betaj = (params->betaj);
	double xc = (params->xc);
	double yc = (params->yc);
	double xa = (params->xa);
	double ya = (params->ya);

	double a = -sin(beta) * xj(xc, xa, s, betaj) + cos(beta) * yj(yc, ya, s, betaj);
	double b = pow(xj(xc, xa, s, betaj), 2) + pow(yj(yc, ya, s, betaj), 2);

	return a / b;
}

double vortex_sheet_function_lastelement(double s, void * parameters)
{
	struct integration_function_parameters * params = (struct integration_function_parameters *)parameters;
	double beta = (params->beta);
	double betaj = (params->betaj);
	double xc = (params->xc);
	double yc = (params->yc);
	double xa = (params->xa);
	double ya = (params->ya);

	double a = -xj(xc, xa, s, betaj) * cos(beta) - yj(yc, ya, s, betaj) * sin(beta);
	double b = pow(xj(xc, xa, s, betaj), 2) + pow(yj(yc, ya, s, betaj), 2);

	return a / b;
}

double v_t_source_function(double s, void * parameters)
{
	struct integration_function_parameters * params = (struct integration_function_parameters *)parameters;
	double beta = (params->beta);
	double betaj = (params->betaj);
	double xc = (params->xc);
	double yc = (params->yc);
	double xa = (params->xa);
	double ya = (params->ya);
	double a = (-xj(xc, xa, s, betaj) * sin(beta) + yj(yc, ya, s, betaj) * cos(beta)) ;
	double b = (pow(xj(xc, xa, s, betaj), 2.) + pow(yj(yc, ya, s, betaj), 2.));

	return a / b;
}

double v_t_vortex_function(double s, void * parameters)
{
	struct integration_function_parameters * params = (struct integration_function_parameters *)parameters;
	double beta = (params->beta);
	double betaj = (params->betaj);
	double xc = (params->xc);
	double yc = (params->yc);
	double xa = (params->xa);
	double ya = (params->ya);
	double a = xj(xc, xa, s, betaj) * cos(beta) + yj(yc, ya, s, betaj) * sin(beta);
	double b = (pow(xj(xc, xa, s, betaj), 2.) + pow(yj(yc, ya, s, betaj), 2.));

	return a / b;
}

double calculate_cl(const wif_core::airfoil_c & airfoil, double U_inf, double gamma)
{
	if(!airfoil.is_valid())
	{
		return 0.0;
	}

	std::vector<wif_core::line_2d_c> lines = airfoil.get_lines();
	std::vector<wif_core::vector_2d_c> points = airfoil.get_points();

	double length_sum = std::accumulate(lines.begin(), lines.end(), 0.0, [](double init, const wif_core::line_2d_c & b)
	{
		return init + b.get_length();
	});

	double x_max = (*std::max_element(points.begin(), points.end(), [](const wif_core::vector_2d_c & lhs, const wif_core::vector_2d_c & rhs)
	{
		return lhs.x < rhs.x;
	})).x;
	double x_min = (*std::min_element(points.begin(), points.end(), [](const wif_core::vector_2d_c & lhs, const wif_core::vector_2d_c & rhs)
	{
		return lhs.x < rhs.x;
	})).x;

	return length_sum * gamma * 2.0 / (U_inf * (x_max - x_min));
}

calculation_results_c calculate_flow(const wif_core::airfoil_c & myAirfoil, std::shared_ptr<wif_core::uniform_flow_c> myFlow, bool Kutta, double gamma)
{
	calculation_results_c c;
	double pi = M_PI;

	std::vector<wif_core::line_2d_c> mylines = myAirfoil.get_lines();
	int num_lines = mylines.size();




	std::vector<double> lengths(num_lines);
	std::vector<wif_core::vector_2d_c> centers(num_lines);
	std::vector<double> angles(num_lines);
	std::vector<wif_core::vector_2d_c> points_airfoil = myAirfoil.get_points();

	for(int i = 0; i < num_lines; i++)
	{
		wif_core::line_2d_c temp_line = mylines[i];
		lengths[i] = temp_line.get_length();
		centers[i] = temp_line.get_center_point();

		wif_core::vector_2d_c a((temp_line.end.y - temp_line.begin.y), -(temp_line.end.x - temp_line.begin.x));

		angles[i] = a.get_angle();
	}





	double U_inf = myFlow->get_strength();
	double angle_attack = myFlow->get_angle();
	double result, error;
	double s_0 = 0.;
	std::vector<double> c_p(num_lines);
	std::vector<double> v_ti(num_lines);
	std::shared_ptr<wif_core::flow_accumulate_c> accumulate_flow = std::make_shared<wif_core::flow_accumulate_c>();
	size_t nevals;
	double v_t_i;
	gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(100);

	std::vector<double> Sigma(num_lines);

	//Check if one uses Kutta condition or not
	if(!Kutta)
	{
		int num_rows = num_lines;
		int num_columns = num_lines;

		double matrix_A_data [num_rows * num_columns];
		double vector_b_data [num_columns];

		gsl_matrix_view matrix_A_view = gsl_matrix_view_array(matrix_A_data, num_rows, num_columns);
		struct integration_function_parameters parameters;
		gsl_function FUNC;
		FUNC.function = &source_sheet_function;
		FUNC.params = &parameters;

		//Setting matrix A and vector B to solve the system
		for(int i = 0; i < num_rows; i++)
		{
			vector_b_data[i] = -U_inf * (cos(angles[i]) * cos(angle_attack) + sin(angles[i]) * sin(angle_attack));

			for(int j = 0; j < num_columns; j++)
			{
				if(i == j)
				{
					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, 0.5);
				}
				else
				{
					parameters = {angles[i], angles[j], centers[i].x, centers[i].y, mylines[j].begin.x, mylines[j].begin.y};

					gsl_integration_cquad(&FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, result / (2.*pi));
				}
			}
		}

		//Solve the system
		gsl_vector_view vector_b_view
		    = gsl_vector_view_array(vector_b_data, num_columns);

		gsl_vector * x = gsl_vector_alloc(num_columns);

		int q;

		gsl_permutation * p = gsl_permutation_alloc(num_columns);

		gsl_linalg_LU_decomp(&matrix_A_view.matrix, p, &q);

		gsl_linalg_LU_solve(&matrix_A_view.matrix, p, &vector_b_view.vector, x);

		for(int i = 0; i < num_columns; i++)
		{
			Sigma[i] = gsl_vector_get(x, i);
		}

		gsl_permutation_free(p);
		gsl_vector_free(x);
	} // if (Kutta)
	else
	{
		int num_rows = num_lines + 1;
		int num_columns = num_lines + 1;
		double matrix_A_data [num_rows * num_columns];
		double vector_b_data [num_columns];
		int k = 0; //first panel
		int l = num_lines - 1; //last panel
		struct integration_function_parameters parameters;

		gsl_matrix_view matrix_A_view
		    = gsl_matrix_view_array(matrix_A_data, num_rows, num_columns);

		// CLEANED UP MATRIX CONSTRUCTION
		gsl_function FUNC;
		FUNC.params = &parameters;

		for(int i = 0; i < num_rows; i++)
		{
			if(i < num_lines)
			{
				vector_b_data[i] = - U_inf * cos(angle_attack - angles[i]);
			}
			else
			{
				vector_b_data[i] = - U_inf * sin(angle_attack - angles[k]) - U_inf * sin(angle_attack - angles[l]);
			}

			for(int j = 0; j < num_columns; j++)
			{
				if(i < num_lines && j < num_lines)
				{
					//SOURCE SHEET BLOK (Np x Np)
					FUNC.function = &source_sheet_function;

					if(i == j)
					{
						gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, 0.5);
					}
					else
					{
						parameters = {angles[i], angles[j], centers[i].x, centers[i].y, mylines[j].begin.x, mylines[j].begin.y};

						gsl_integration_cquad(&FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

						gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, result / (2.*pi));
					}

				}
				else if(i < num_lines && j == num_lines)
				{

					//VORTEX_SHEET_1 BLOK (EXTRA KOLOM) (: x Np+1 )
					FUNC.function = &vortex_sheet_function_1;
					double last_column_value = 0;

					for(int r = 0; r < num_lines; r++)
					{
						if(r != i)
						{
							parameters = {angles[i], angles[r], centers[i].x, centers[i].y, mylines[r].begin.x, mylines[r].begin.y};
							gsl_integration_cquad(&FUNC, s_0, lengths[r], 0., 1e-7, w, &result, &error, &nevals);
							last_column_value += result;
						}

					}

					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, -last_column_value / (2.0 * pi));

				}
				else if(i == num_lines && j < num_lines)
				{

					//LAATSTE RIJ BLOK (Np+1 x Np)
					double last_row_value = 0;
					FUNC.function = &vortex_sheet_function_lastrow;

					if(j != k)
					{
						parameters = {angles[k], angles[j], centers[k].x, centers[k].y, mylines[j].begin.x, mylines[j].begin.y};
						gsl_integration_cquad(&FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
						last_row_value += result;
					}

					if(j != l)
					{
						parameters = {angles[l], angles[j], centers[l].x, centers[l].y, mylines[j].begin.x, mylines[j].begin.y};
						gsl_integration_cquad(&FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
						last_row_value += result;
					}

					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, last_row_value / (2.*pi));

				}
				else if(i == num_lines && j == num_lines)
				{

					//LAATSTE ELEMENTJE IN RECHTER ONDER HOEK
					//NEED NICK (Np+1,Np+1)

					FUNC.function = &vortex_sheet_function_lastelement;
					double last_element_value = 0;

					for(int r = 0; r < num_lines; r++)
					{
						if(r != k)
						{
							parameters = {angles[k], angles[r], centers[k].x, centers[k].y, mylines[r].begin.x, mylines[r].begin.y};
							gsl_integration_cquad(&FUNC, s_0, lengths[r], 0., 1e-7, w, &result, &error, &nevals);
							last_element_value += result / (2.0 * pi);
						}
						else
						{
							last_element_value += -0.5;
						}

						if(r != l)
						{
							parameters = {angles[l], angles[r], centers[l].x, centers[l].y, mylines[r].begin.x, mylines[r].begin.y};
							gsl_integration_cquad(&FUNC, s_0, lengths[r], 0., 1e-7, w, &result, &error, &nevals);
							last_element_value += result / (2.0 * pi);
						}
						else
						{
							last_element_value += -0.5;
						}
					}

					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, last_element_value);
				}
			}
		}

		//Solve the system

		gsl_vector_view vector_b_view
		    = gsl_vector_view_array(vector_b_data, num_columns);

		gsl_vector * x = gsl_vector_alloc(num_columns);

		int q;

		gsl_permutation * p = gsl_permutation_alloc(num_columns);

		gsl_linalg_LU_decomp(&matrix_A_view.matrix, p, &q);

		gsl_linalg_LU_solve(&matrix_A_view.matrix, p, &vector_b_view.vector, x);

		for(int i = 0; i < num_lines; i++)
		{
			Sigma[i] = gsl_vector_get(x, i);
		}

		gamma = gsl_vector_get(x, num_lines);

		std::cout << "Gamma is " << gamma << std::endl;

		gsl_permutation_free(p);
		gsl_vector_free(x);
	} // else kutta

	//Calculating c_p
	for(int i = 0; i < num_lines; i++)
	{
		v_t_i = U_inf * (- sin(angles[i]) * cos(angle_attack) + cos(angles[i]) * sin(angle_attack));

		for(int j = 0; j < num_lines; j++)
		{
			if(i != j)
			{
				struct integration_function_parameters parameters_2;
				gsl_function V_FUNC;
				V_FUNC.params = &parameters_2;
				parameters_2 = {angles[i], angles[j], centers[i].x, centers[i].y, mylines[j].begin.x, mylines[j].begin.y};

				V_FUNC.function = &v_t_source_function;
				gsl_integration_cquad(&V_FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
				v_t_i += (Sigma[j] / (2.0 * M_PI)) * result;

				V_FUNC.function = &v_t_vortex_function;
				gsl_integration_cquad(&V_FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
				v_t_i -= (gamma / (2.0 * M_PI)) * result;
			}
		}

		v_ti[i] = v_t_i;
		c_p[i] = 1.0 - pow(v_t_i / U_inf, 2);
	}

	accumulate_flow->add_flow(myFlow);
	accumulate_flow->add_source_sheets(Sigma, myAirfoil);

	if(gamma != 0)
	{
		accumulate_flow->add_vortex_sheets(gamma, myAirfoil);
	}

	c.airfoil = myAirfoil;
	c.c_p = c_p;
	c.c_l = calculate_cl(myAirfoil, U_inf, gamma);
	c.flow = accumulate_flow;
	c.v_t = v_ti;
	c.closed_body_check = std::inner_product(lengths.begin(), lengths.end(), Sigma.begin(), 0.0);

	std::cout << "Closed body check: " << c.closed_body_check << std::endl;

	gsl_integration_cquad_workspace_free(w);
	return c;
}

} // namespace wif_algo
