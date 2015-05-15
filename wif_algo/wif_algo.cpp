#include "wif_algo.hpp"
#include <gsl/gsl_linalg.h>
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

double source_sheet_function(double s, void * parameters)
{
	struct integration_function_parameters * params = (struct integration_function_parameters *)parameters;
	double beta = (params->beta);
	double betaj = (params->betaj);
	double xc = (params->xc);
	double yc = (params->yc);
	double xa = (params->xa);
	double ya = (params->ya);
	double source_sheet_function = ((xc - (xa - s * sin(betaj))) * cos(beta) + (yc - (ya + s * cos(betaj))) * sin(beta)) / (pow((xc - (xa - s * sin(betaj))), 2.) + pow((yc - (ya + s * cos(betaj))), 2.));
	return source_sheet_function;
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

	double a = sin(beta) * (xc - (xa - s * sin(betaj))) - cos(beta) * (yc - (ya + s * cos(betaj)));
	double b = pow(xc - (xa - s * sin(betaj)), 2) + pow(yc - (ya + s * cos(betaj)), 2);

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

	double a = -cos(beta) * (xc - (xa - s * sin(betaj))) - sin(beta) * (yc + (ya + s * cos(betaj)));
	double b = pow(xc - (xa - s * sin(betaj)), 2) + pow(yc - (ya + s * cos(betaj)), 2);

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

	double a = -(xc - (xa - s * sin(betaj))) * cos(beta)  - (yc + (ya + s * cos(betaj))) * sin(beta);
	double b = pow(xc - (xa - s * sin(betaj)), 2) + pow(yc - (ya + s * cos(betaj)), 2);

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
	double a = (-(xc - (xa - s * sin(betaj))) * sin(beta) + (yc - (ya + s * cos(betaj))) * cos(beta)) ;
	double b = (pow((xc - (xa - s * sin(betaj))), 2.) + pow((yc - (ya + s * cos(betaj))), 2.));

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
	double a = ((xc - (xa - s * sin(betaj))) * cos(beta) + (yc - (ya + s * cos(betaj))) * sin(beta)) ;
	double b = (pow((xc - (xa - s * sin(betaj))), 2.) + pow((yc - (ya + s * cos(betaj))), 2.));

	return a / b;
}

calculation_results_c calculate_flow(const wif_core::airfoil_c & myAirfoil, std::shared_ptr<wif_core::uniform_flow_c> myFlow, bool Kutta, double gamma)
{
	calculation_results_c c;
	double pi = 3.1415;

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

		if(temp_line.begin.x > temp_line.end.x && temp_line.begin.y > temp_line.end.y)
		{
			angles[i] = temp_line.get_angle() - pi / 2;
		}
		else if(temp_line.begin.x > temp_line.end.x && temp_line.begin.y < temp_line.end.y)
		{
			angles[i] = temp_line.get_angle() - pi / 2;
		}
		else if(temp_line.begin.x < temp_line.end.x && temp_line.begin.y > temp_line.end.y)
		{
			angles[i] = temp_line.get_angle() - pi / 2;
		}
		else if(temp_line.begin.x < temp_line.end.x && temp_line.begin.y < temp_line.end.y)
		{
			angles[i] = temp_line.get_angle() + (3 * pi) / 2;
		}
		else if(temp_line.begin.x == temp_line.end.x)
		{
			angles[i] = 0;
		}
		else if(temp_line.begin.y == temp_line.end.y)
		{
			angles[i] = pi / 2;
		}
		else
		{
			std::cerr << i << "  Not in any of these categories" << std::endl;
		}

	}





	double U_inf = myFlow->get_strength();
	double angle_attack = myFlow->get_angle();
	double result, error;
	double s_0 = 0.;
	std::vector<double> c_p(num_lines);
	double c_l;
	std::shared_ptr<wif_core::flow_accumulate_c> accumulate_flow = std::make_shared<wif_core::flow_accumulate_c>();
	size_t nevals;
	double v_t_i;
	gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(100);


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
					parameters = {angles[i], angles[j], centers[i].x, centers[i].y, points_airfoil[j].x, points_airfoil[j].y};

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

		std::vector<double> Sigma(num_columns);

		for(int i = 0; i < num_columns; i++)
		{
			Sigma[i] = gsl_vector_get(x, i);
		}

		gsl_permutation_free(p);
		gsl_vector_free(x);

		//Calculating c_p
		gsl_function V_FUNC;
		V_FUNC.params = &parameters;

		for(int i = 0; i < num_lines; i++)
		{
			v_t_i = 0.;

			for(int j = 0; j < num_lines; j++)
			{
				if(i != j)
				{
					parameters = {angles[i], angles[j], centers[i].x, centers[i].y, points_airfoil[j].x, points_airfoil[j].y};

					V_FUNC.function = &v_t_source_function;
					gsl_integration_cquad(&V_FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
					v_t_i += ((Sigma[j] / (2 * pi)) * result);

					V_FUNC.function = &v_t_vortex_function;
					gsl_integration_cquad(&V_FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
					v_t_i -= gamma / (2 * pi) * result;
				}

			}

			c_p[i] = 1.0 - pow((U_inf * (-sin(angles[i]) * cos(angle_attack) + cos(angles[i]) * sin(angle_attack)) + v_t_i) / U_inf, 2);
		}

		//Calculate c_l
		double ll = std::accumulate(lengths.begin(), lengths.end(), 0);
		double x_max = points_airfoil[0].x;
		double x_min = points_airfoil[0].x;

		for(int j = 1; j < num_lines; j++)
		{
			if(points_airfoil[j].x > x_max)
			{
				x_max = points_airfoil[j].x;
			}
			else if(points_airfoil[j].x < x_min)
			{
				x_min = points_airfoil[j].x;
			}
		}

		c_l = ll * gamma * 2 / U_inf / (x_max - x_min);


		//Building return

		accumulate_flow->add_flow(myFlow);
		accumulate_flow->add_source_sheets(Sigma, myAirfoil);
		accumulate_flow->add_vortex_sheets(gamma, myAirfoil);
		c.airfoil = myAirfoil;
		c.flow = accumulate_flow;
		c.c_p = c_p;
		c.c_l = c_l;

		////
		double E = 0;

		for(int i = 0; i < num_lines; i++)
		{
			E += lengths[i] * Sigma[i];
		}

		c.closed_body_check = E;

	} // if (Kutta)
	else
	{
		int num_rows = num_lines + 1;
		int num_columns = num_lines + 1;
		double matrix_A_data [num_rows * num_columns];
		double vector_b_data [num_columns];
		int k = 0; //first panel
		int l = num_lines-1; //last panel
		struct integration_function_parameters parameters;

		gsl_matrix_view matrix_A_view
		    = gsl_matrix_view_array(matrix_A_data, num_rows, num_columns);

		// CLEANED UP MATRIX CONSTRUCTION
		gsl_function FUNC;
		FUNC.params = &parameters;

		for(int i = 0; i < num_rows; i++)
		{
			for(int j = 0; j < num_columns; j++)
			{
				if(i < num_lines && j < num_lines)
				{
					//SOURCE SHEET BLOK (Np x Np)
					FUNC.function = &source_sheet_function;

					vector_b_data[i] = -U_inf * cos(angle_attack - angles[i]);

					if(i == j)
					{
						gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, 0.5);
					}
					else
					{
						parameters = {angles[i], angles[j], centers[i].x, centers[i].y, points_airfoil[j].x, points_airfoil[j].y};

						gsl_integration_cquad(&FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

						gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, result / (2.*pi));
					}

				}
				else if(i < num_lines && j == num_lines)
				{

					//VORTEX_SHEET_1 BLOK (EXTRA KOLOM) (: x Np+1 )
					vector_b_data[i] = -U_inf * cos(angle_attack - angles[i]);
					FUNC.function = &vortex_sheet_function_1;
					double last_column_value = 0;

					for(int r = 0; r < num_lines; r++)
					{
						if(r != i)
						{
							parameters = {angles[i], angles[r], centers[i].x, centers[i].y, points_airfoil[r].x, points_airfoil[r].y};
							gsl_integration_cquad(&FUNC, s_0, lengths[r], 0., 1e-7, w, &result, &error, &nevals);
							last_column_value += result;
						}

					}

					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, -result / (2.*pi));

				}
				else if(i == num_lines && j < num_lines)
				{

					//LAATSTE RIJ BLOK (Np+1 x Np)
					//vector_b_data[i] = -U_inf * cos(angles[k] - angle_attack) - U_inf * cos(angles[l] - angle_attack);
					vector_b_data[i] = -U_inf * sin(angle_attack - angles[k]) - U_inf * sin(angle_attack - angles[l]);
					double last_row_value = 0;
					FUNC.function = &vortex_sheet_function_lastrow;

					if(j != k)
					{
						parameters = {angles[k], angles[j], centers[k].x, centers[k].y, points_airfoil[j].x, points_airfoil[j].y};
						gsl_integration_cquad(&FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
						last_row_value += result;
					}

					if(j != l)
					{
						parameters = {angles[l], angles[j], centers[l].x, centers[l].y, points_airfoil[j].x, points_airfoil[j].y};
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
							parameters = {angles[k], angles[r], centers[k].x, centers[k].y, points_airfoil[r].x, points_airfoil[r].y};
							gsl_integration_cquad(&FUNC, s_0, lengths[r], 0., 1e-7, w, &result, &error, &nevals);
							last_element_value += result / (2 * pi);
						}
						else
						{
							last_element_value += -0.5;
						}

						if(r != l)
						{
							parameters = {angles[l], angles[r], centers[l].x, centers[l].y, points_airfoil[r].x, points_airfoil[r].y};
							gsl_integration_cquad(&FUNC, s_0, lengths[r], 0., 1e-7, w, &result, &error, &nevals);
							last_element_value += result / (2 * pi);
						}
						else
						{
							last_element_value += -0.5;
						}
					}


					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, last_element_value);

				}
				else
				{

					std::cout << i << "  " << j << "  Deze indices zitten er niet bij ?" << std::endl;

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

		std::vector<double> Sigma(num_lines);

		for(int i = 0; i < num_lines; i++)
		{
			Sigma[i] = gsl_vector_get(x, i);
		}

		double Gamma = gsl_vector_get(x, num_lines);
		//printf("a = \n");
		//gsl_matrix_fprintf(stdout, &matrix_A_view.matrix, "%g");

		gsl_permutation_free(p);
		gsl_vector_free(x);



		//Calculating c_p
		gsl_function V_FUNC;
		V_FUNC.params = &parameters;

		for(int i = 0; i < num_lines; i++)
		{
			v_t_i = 0.;

			for(int j = 0; j < num_lines; j++)
			{
				if(i != j)
				{
					parameters = {angles[i], angles[j], centers[i].x, centers[i].y, points_airfoil[j].x, points_airfoil[j].y};

					V_FUNC.function = &v_t_source_function;
					gsl_integration_cquad(&V_FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
					v_t_i += ((Sigma[j] / (2 * pi)) * result);

					V_FUNC.function = &v_t_vortex_function;
					gsl_integration_cquad(&V_FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
					v_t_i -= (Gamma / (2 * pi)) * result;

				}

			}

			c_p[i] = 1.0 - pow((U_inf * (-sin(angles[i]) * cos(angle_attack) + cos(angles[i]) * sin(angle_attack)) + v_t_i) / U_inf, 2);
		}

		//Calculate c_l
		double ll = std::accumulate(lengths.begin(), lengths.end(), 0);
		double x_max = points_airfoil[0].x;
		double x_min = points_airfoil[0].x;

		for(int j = 1; j < num_lines; j++)
		{
			if(points_airfoil[j].x > x_max)
			{
				x_max = points_airfoil[j].x;
			}
			else if(points_airfoil[j].x < x_min)
			{
				x_min = points_airfoil[j].x;
			}
		}

		c_l = ll * Gamma * 2 / U_inf / (x_max - x_min);

		//Building return
		accumulate_flow->add_flow(myFlow);
		accumulate_flow->add_source_sheets(Sigma, myAirfoil);
		accumulate_flow->add_vortex_sheets(Gamma, myAirfoil);

		std::cout << "Gamma is " << Gamma << std::endl;

		c.airfoil = myAirfoil;
		c.flow = accumulate_flow;
		c.c_p = c_p;
		c.c_l = c_l;

		////
		double E = 0;

		for(int i = 0; i < num_lines; i++)
		{
			E += lengths[i] * Sigma[i];
		}

		c.closed_body_check = E;
	} // else kutta



	gsl_integration_cquad_workspace_free(w);
	return c;
}

} // namespace wif_algo
