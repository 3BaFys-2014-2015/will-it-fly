#include "wif_algo.hpp"

#include <wif_core/wif_core.hpp>

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

double vortex_sheet_function_2(double s, void * parameters)
{
	struct integration_function_parameters * params = (struct integration_function_parameters *)parameters;
	double beta = (params->beta);
	double betaj = (params->betaj);
	double xc = (params->xc);
	double yc = (params->yc);
	double xa = (params->xa);
	double ya = (params->ya);

	double a = cos(beta) * (xc - (xa - s * sin(betaj))) - sin(beta) * (yc - (ya + s * cos(betaj)));
	double b = pow(xc - (xa - s * sin(betaj)), 2) + pow(yc - (ya + s * cos(betaj)), 2);

	return a / b;
}

double v_t_function(double s, void * parameters)
{
	struct integration_function_parameters * params = (struct integration_function_parameters *)parameters;
	double beta = (params->beta);
	double betaj = (params->betaj);
	double xc = (params->xc);
	double yc = (params->yc);
	double xa = (params->xa);
	double ya = (params->ya);
	double a = ((xc - (xa - s * sin(betaj))) * -sin(beta) + (yc - (ya + s * cos(betaj))) * cos(beta)) ;
	double b = (pow((xc - (xa - s * sin(betaj))), 2.) + pow((yc - (ya + s * cos(betaj))), 2.));

	return a / b;
}

calculation_results_c calculate_flow(const wif_core::airfoil_c & myAirfoil, std::shared_ptr<wif_core::uniform_flow_c> myFlow, bool Kutta)
{
	calculation_results_c c;
	double pi = 3.1415;

	std::vector<wif_core::line_2d_c> mylines = myAirfoil.get_lines();
	unsigned int num_lines = mylines.size();


	std::vector<double> lengths(num_lines);
	std::vector<wif_core::vector_2d_c> centers(num_lines);
	std::vector<double> angles(num_lines);
	std::vector<wif_core::vector_2d_c> points_airfoil = myAirfoil.get_points();

	for(unsigned int i = 0; i < num_lines; i++)
	{
		wif_core::line_2d_c temp_line = mylines[i];
		lengths[i] = temp_line.get_length();
		centers[i] = temp_line.get_center_point();

		if(centers[i].y > 0)
		{
			angles[i] = atan2(centers[i].y, centers[i].x);
		}
		else
		{
			angles[i] = atan2(centers[i].y, centers[i].x) + 2 * pi;
		}

	}





	double U_inf = myFlow->get_strength();
	double angle_attack = myFlow->get_angle();
	double result, error;
	double s_0 = 0.;
	vector<double> c_p(num_lines);
	double c_l;
	std::shared_ptr<wif_core::flow_accumulate_c> accumulate_flow = std::make_shared<wif_core::flow_accumulate_c>();
	size_t nevals;
	double v_t_i;
	gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(100);


	//Check if one uses Kutta condition or not
	if(!Kutta)
	{
		unsigned int num_rows = num_lines;
		unsigned int num_columns = num_lines;

		double matrix_A_data [num_rows * num_columns];
		double vector_b_data [num_columns];


		gsl_matrix_view matrix_A_view = gsl_matrix_view_array(matrix_A_data, num_rows, num_columns);

		gsl_function FUNC;
		FUNC.function = &source_sheet_function;

		//Setting matrix A and vector B to solve the system
		for(unsigned int i = 0; i < num_rows; i++)
		{
			vector_b_data[i] = -U_inf * cos(angles[i] - angle_attack);

			for(unsigned int j = 0; j < num_columns; j++)
			{
				if(i == j)
				{
					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, 0.5);

				}
				else
				{
					struct integration_function_parameters parameters = {angles[i], angles[j], centers[i].x, centers[i].y, points_airfoil[j].x, points_airfoil[j].y};

					FUNC.params = &parameters;

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

		for(unsigned int i = 0; i < num_columns; i++)
		{
			Sigma[i] = gsl_vector_get(x, i);
		}



		gsl_permutation_free(p);
		gsl_vector_free(x);

		//Calculating c_p
		gsl_function V_FUNC;
		V_FUNC.function = &v_t_function;

		for(unsigned int i = 0; i < num_lines; i++)
		{
			v_t_i = 0.;

			for(unsigned int j = 0; j < num_lines; j++)
			{
				if(i == j)
				{
				}
				else
				{
					struct integration_function_parameters parameters = {angles[i], angles[j], centers[i].x, centers[i].y, points_airfoil[j].x, points_airfoil[j].y};

					V_FUNC.params = &parameters;

					gsl_integration_cquad(&V_FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

					v_t_i = v_t_i + ((Sigma[j] / (2 * pi)) * result);
					//std::cout << Sigma[j]*result << "  " << i << "  " << j << std::endl;
				}

			}

			std::cout << v_t_i << std::endl;
			c_p[i] = 1 - pow((-U_inf * cos(angles[i] - angle_attack) + v_t_i) / U_inf, 2);
		}



		//calculating c_l
		c_l = 0;

		//Building return

		accumulate_flow->add_flow(myFlow);

		accumulate_flow->add_source_sheets(Sigma, myAirfoil);

		c.airfoil = myAirfoil;
		c.flow = accumulate_flow;
		c.c_p = c_p;
		c.c_l = c_l;

		////
		double E = 0;

		auto i = lengths.size();

		for(i = 0; i < lengths.size(); i++)
		{
			E = E + lengths[i] * Sigma[i];
		}

		std::cout << "Voor een gesloten lichaam moet de som van alle source sterktes gelijk zijn aan nul, vergelijking (31):";
		std::cout << E << std::endl;

	} // if (Kutta)
	else
	{
		unsigned int num_rows = num_lines + 1;
		unsigned int num_columns = num_lines + 1;
		double matrix_A_data [num_rows * num_columns];
		double vector_b_data [num_columns];
		int k = 0; //first panel
		int l = num_rows - 1; //last panel

		gsl_matrix_view matrix_A_view
		    = gsl_matrix_view_array(matrix_A_data, num_rows, num_columns);

		gsl_function FUNC;
		FUNC.function = &source_sheet_function;

		//First set matrix a and vector b
		for(unsigned int i = 0; i < num_rows - 1; i++)
		{
			vector_b_data[i] = -U_inf * cos(angles[i] - angle_attack);

			for(unsigned int j = 0; j < num_columns - 1; j++)
			{
				if(i == j)
				{
					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, 0.5);

				}
				else
				{
					struct integration_function_parameters parameters = {angles[i], angles[j], centers[i].x, centers[i].y, points_airfoil[j].x, points_airfoil[j].y};

					FUNC.params = &parameters;


					gsl_integration_cquad(&FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

					gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, result / (2.*pi));

				}
			}
		}

		gsl_function VORTEX1;
		VORTEX1.function = &vortex_sheet_function_1;

		//Set last collumn of matrix A
		for(unsigned int i = 0; i < num_rows; i++)
		{
			unsigned int j = num_columns - 1;

			struct integration_function_parameters parameters = {angles[i], angles[j], centers[i].x, centers[i].y, points_airfoil[j].x, points_airfoil[j].y};

			VORTEX1.params = &parameters;

			gsl_integration_cquad(&VORTEX1, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

			gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, -0.5 / pi * result);
		}

		//Set last row of matrix A and set vector B
		for(unsigned int j = 0; j < num_columns; j++)
		{
			unsigned int i = num_rows - 1;

			vector_b_data[i] = -U_inf * cos(angle_attack - angles[k]) - U_inf * cos(angle_attack - angles[l]);

			struct integration_function_parameters parameters = {angles[k], angles[j], centers[k].x, centers[k].y, points_airfoil[j].x, points_airfoil[j].y};

			VORTEX1.params = &parameters;
			gsl_integration_cquad(&VORTEX1, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

			double result_temp = -0.5 / pi * result;

			parameters = {angles[l], angles[j], centers[l].x, centers[l].y, points_airfoil[j].x, points_airfoil[j].y};

			VORTEX1.params = &parameters;

			size_t nevals;
			gsl_integration_cquad(&VORTEX1, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

			result_temp = result_temp - 0.5 / (pi) * result;

			gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, result_temp);

		}

		//Set last ellement of A

		int i = num_rows - 1;
		int j = num_columns - 1;


		gsl_function VORTEX2;
		VORTEX2.function = &vortex_sheet_function_2;

		struct integration_function_parameters parameters = {angles[k], angles[j], centers[k].x, centers[k].y, points_airfoil[j].x, points_airfoil[j].y};

		VORTEX2.params = &parameters;

		gsl_integration_cquad(&VORTEX2, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

		double result_temp = -0.5 / (pi) * result;

		parameters = {angles[l], angles[j], centers[l].x, centers[l].y, points_airfoil[j].x, points_airfoil[j].y};

		VORTEX2.params = &parameters;

		size_t nevals;
		gsl_integration_cquad(&VORTEX2, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);

		result_temp = result_temp - 0.5 / (pi) * result;

		gsl_matrix_set(&matrix_A_view.matrix, (size_t) i, (size_t) j, result_temp);

		//Solve the system
		gsl_vector_view vector_b_view
		    = gsl_vector_view_array(vector_b_data, num_columns);

		gsl_vector * x = gsl_vector_alloc(num_columns);

		int q;

		gsl_permutation * p = gsl_permutation_alloc(num_columns);

		gsl_linalg_LU_decomp(&matrix_A_view.matrix, p, &q);

		gsl_linalg_LU_solve(&matrix_A_view.matrix, p, &vector_b_view.vector, x);

		std::vector<double> Sigma(num_columns - 1);

		for(unsigned int i = 0; i < num_columns - 1; i++)
		{
			Sigma[i] = gsl_vector_get(x, i);
		}

		double Gamma;
		Gamma = gsl_vector_get(x, num_columns - 1);

		gsl_permutation_free(p);
		gsl_vector_free(x);

		//Calculate c_p
		gsl_function V_FUNC;
		V_FUNC.function = &v_t_function;

		double v_t_i;

		for(unsigned int i = 0; i < num_lines; j++)
		{
			v_t_i = 0;

			for(unsigned int j = 0; j < num_lines; j++)
			{
				struct integration_function_parameters parameters = {angles[i], angles[j], centers[i].x, centers[i].y, points_airfoil[j].x, points_airfoil[j].y};

				V_FUNC.params = &parameters;

				size_t nevals;
				gsl_integration_cquad(&V_FUNC, s_0, lengths[j], 0., 1e-7, w, &result, &error, &nevals);
				v_t_i = v_t_i + Sigma[j] / (2 * pi) * result;

			}

			c_p[i] = 1 - pow((-U_inf * sin(angles[i]) + v_t_i) / U_inf, 2);
		}

		//Calculate c_l
		double ll = std::accumulate(lengths.begin(), lengths.end(), 0);
		double x_max = points_airfoil[0].x;
		double x_min = points_airfoil[0].x;

		for(unsigned int j = 1; j < points_airfoil.size(); j++)
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
		c.airfoil = myAirfoil;
		c.flow = accumulate_flow;
		c.c_p = c_p;
		c.c_l = c_l;
	} // else kutta



	gsl_integration_cquad_workspace_free(w);
	return c;
}

} // namespace wif_algo
