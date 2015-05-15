#include <iostream>
#include <cmath>
#include <fstream>

#include <wif_core/wif_core.hpp>
#include <wif_algo/wif_algo.hpp>
#include <wif_viz/wif_viz.hpp>

int main()
{
	double pi = 3.1415;
	//opdracht 8
	/*double radius = 2;
	unsigned int num_lines = 100;

	wif_core::vector_2d_c midpoint(0, 0);
	wif_core::airfoil_c myAirfoil(midpoint, radius, num_lines);
	std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);
	std::vector<wif_core::line_2d_c> mylines = myAirfoil.get_lines();
	std::vector<wif_core::vector_2d_c> centers(num_lines);
	std::vector<double> angles(num_lines);
	bool Kutta = 0;
	int imax = 20;
	int stepsize = 5;
	int j = 0;
	std::vector<vector<double>> data_to_plot_gamma(imax / stepsize, std::vector<double>(num_lines));
	std::cout << imax / stepsize << std::endl;
	std::vector<std::string> Legend(imax / stepsize);
	std::vector<double> X_as(num_lines);

	for(unsigned int i = 0; i < num_lines; i++)
	{
		wif_core::line_2d_c temp_line = mylines[i];
		centers[i] = temp_line.get_center_point();
		X_as[i] = i;

		if(centers[i].y > 0)
		{
			angles[i] = atan2(centers[i].y, centers[i].x);
		}
		else
		{
			angles[i] = atan2(centers[i].y, centers[i].x) + 2 * pi;
		}

	}




	for(int i = 0 ; i < imax ; i = i + stepsize)
	{
		double gamma = i * 1.;

		std::string filename = std::to_string(i) + ".txt";
		wif_algo::calculation_results_c calculate_flow2 = wif_algo::calculate_flow(myAirfoil, myFlow, Kutta, gamma);


		data_to_plot_gamma[j] = calculate_flow2.c_p;

		Legend[j] = std::to_string(gamma);

		j++;

	}

	std::shared_ptr<wif_viz::visualization_c> myRoot = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
	myRoot->plotVectors(data_to_plot_gamma, X_as, Legend, "test.pdf", "paneel", "cp", "cp(gamma)");*/


	//Opdracht 14
	/*double radius = 2;
	wif_core::vector_2d_c midpoint(0, 0);
	wif_core::airfoil_c myAirfoil("wif_core/airfoils/n0012-il.dat");
	wif_core::airfoil_c myClosedAirfoil = myAirfoil.closed_merge();
	double gamma = 0;

	bool Kutta = 0;
	std::shared_ptr<wif_core::uniform_flow_c> myFlow2 = std::make_shared<wif_core::uniform_flow_c>(0, 1);

	int jj = 0;

		std::vector<wif_core::line_2d_c> mylines = myClosedAirfoil.get_lines();
		int num_lines = mylines.size();
		std::vector<wif_core::vector_2d_c> centers(num_lines);
		std::vector<double> angles(num_lines);
		std::vector<double> upperfoil_x(num_lines);
		std::vector<double> underfoil_x(num_lines);
		std::vector<double> centers_x(num_lines);
		std::vector<std::vector<double>> data_to_plot_under_upper(2, std::vector<double>(num_lines));
		std::vector<std::vector<double>> data_to_plot_under(5, std::vector<double>(num_lines));
		std::vector<std::vector<double>> data_to_plot_upper(5, std::vector<double>(num_lines));
		std::vector<std::string> Legend(2);
		std::vector<std::string> Legend2(5);
		std::vector<double> X_as(num_lines);
		for(int i = 0 ; i < 5; i++)
		{
			double angle_of_attack = i*(pi/8);
			std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(angle_of_attack, 1);
			wif_algo::calculation_results_c calculate_flow2 = wif_algo::calculate_flow(myClosedAirfoil, myFlow, Kutta);


			for(unsigned int j = 0; j < num_lines; j++)
			{
				wif_core::line_2d_c temp_line = mylines[j];
				centers[j] = temp_line.get_center_point();
				X_as[j]=centers[j].x;

				if(centers[j].y > 0)
				{
					angles[j] = atan2(centers[i].y, centers[i].x);
					upperfoil_x[j]=calculate_flow2.c_p[i];
				}
				else
				{
					angles[j] = atan2(centers[j].y, centers[j].x) + 2 * pi;
					underfoil_x[j]=calculate_flow2.c_p[j];
				}

			}

			data_to_plot_under_upper[0]=underfoil_x;
			data_to_plot_under_upper[1]=upperfoil_x;
			Legend[0]="Under foil with angle of attack " + std::to_string(angle_of_attack);
			Legend[1]="Upper foil with angle of attack " + std::to_string(angle_of_attack);
			if(jj==0)
			{
				std::shared_ptr<wif_viz::visualization_c> myRoot = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
				myRoot->plotVectors(data_to_plot_under_upper, X_as, Legend, "test4.pdf", "x-positie van paneel", "cp", "Cp of panels under the foil");
			}
			data_to_plot_under[jj]=underfoil_x;
			data_to_plot_upper[jj]=upperfoil_x;
			Legend2[jj]=std::to_string(angle_of_attack);
			jj++;

		}
	std::shared_ptr<wif_viz::visualization_c> myRoot = wif_viz::create_visualization_root(myFlow2, midpoint, midpoint);
	myRoot->plotVectors(data_to_plot_under, X_as, Legend2, "test2.pdf", "x-positie van paneel", "cp", "Cp of panels under the foil");
	myRoot->plotVectors(data_to_plot_upper, X_as, Legend2, "test3.pdf", "x-positie van paneel", "cp", "Cp of panels on top of the foil");*/
	//Opdracht 15
	/*double radius = 2;
	wif_core::vector_2d_c midpoint(0, 0);
	wif_core::airfoil_c myAirfoil("wif_core/airfoils/n0012-il.dat");
	wif_core::airfoil_c myClosedAirfoil = myAirfoil.closed_merge();
	double gamma = 0;

	bool Kutta = 1;
	std::shared_ptr<wif_core::uniform_flow_c> myFlow2 = std::make_shared<wif_core::uniform_flow_c>(0, 1);

	int jj = 0;

		std::vector<wif_core::line_2d_c> mylines = myClosedAirfoil.get_lines();
		int num_lines = mylines.size();
		std::vector<wif_core::vector_2d_c> centers(num_lines);
		std::vector<double> angles(num_lines);
		std::vector<double> upperfoil_x(num_lines);
		std::vector<double> underfoil_x(num_lines);
		std::vector<double> centers_x(num_lines);
		std::vector<std::vector<double>> data_to_plot_under_upper(2, std::vector<double>(num_lines));
		std::vector<std::vector<double>> data_to_plot_under(5, std::vector<double>(num_lines));
		std::vector<std::vector<double>> data_to_plot_upper(5, std::vector<double>(num_lines));
		std::vector<std::string> Legend(2);
		std::vector<std::string> Legend2(5);
		std::vector<double> X_as(num_lines);
		for(int i = 0 ; i < 5; i++)
		{
			double angle_of_attack = i*(pi/8);
			std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(angle_of_attack, 1);
			wif_algo::calculation_results_c calculate_flow2 = wif_algo::calculate_flow(myClosedAirfoil, myFlow, Kutta);


			for(unsigned int j = 0; j < num_lines; j++)
			{
				wif_core::line_2d_c temp_line = mylines[j];
				centers[j] = temp_line.get_center_point();
				X_as[j]=centers[j].x;

				if(centers[j].y > 0)
				{
					angles[j] = atan2(centers[i].y, centers[i].x);
					upperfoil_x[j]=calculate_flow2.c_p[i];
				}
				else
				{
					angles[j] = atan2(centers[j].y, centers[j].x) + 2 * pi;
					underfoil_x[j]=calculate_flow2.c_p[j];
				}

			}

			data_to_plot_under_upper[0]=underfoil_x;
			data_to_plot_under_upper[1]=upperfoil_x;
			Legend[0]="Under foil with angle of attack " + std::to_string(angle_of_attack);
			Legend[1]="Upper foil with angle of attack " + std::to_string(angle_of_attack);
			if(jj==0)
			{
				std::shared_ptr<wif_viz::visualization_c> myRoot = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
				myRoot->plotVectors(data_to_plot_under_upper, X_as, Legend, "test5.pdf", "x-positie van paneel", "cp", "Cp of panels under the foil with Kutta");
			}
			data_to_plot_under[jj]=underfoil_x;
			data_to_plot_upper[jj]=upperfoil_x;
			Legend2[jj]=std::to_string(angle_of_attack);
			jj++;

		}
	std::shared_ptr<wif_viz::visualization_c> myRoot = wif_viz::create_visualization_root(myFlow2, midpoint, midpoint);
	myRoot->plotVectors(data_to_plot_under, X_as, Legend2, "test6.pdf", "x-positie van paneel", "cp", "Cp of panels under the foil with Kutta");
	myRoot->plotVectors(data_to_plot_upper, X_as, Legend2, "test7.pdf", "x-positie van paneel", "cp", "Cp of panels on top of the foil witfh Kutta");*/

	//Opdracht 16
	double radius = 2;
	wif_core::vector_2d_c midpoint(0, 0);
	wif_core::airfoil_c myAirfoil("wif_core/airfoils/n0012-il.dat");
	wif_core::airfoil_c myAirfoil2("wif_core/airfoils/lednicer.dat");
	wif_core::airfoil_c myAirfoil3("wif_core/airfoils/selig.dat");
	double gamma = 0;

	bool Kutta = 1;
	std::shared_ptr<wif_core::uniform_flow_c> myFlow2 = std::make_shared<wif_core::uniform_flow_c>(0, 1);
		std::vector<wif_core::line_2d_c> mylines = myAirfoil.get_lines();
		int num_lines = mylines.size();
		std::vector<wif_core::vector_2d_c> centers(num_lines);
		std::vector<double> angles(num_lines);
		std::vector<double> centers_x(num_lines);
		std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);
		wif_algo::calculation_results_c calculate_flow2 = wif_algo::calculate_flow(myAirfoil, myFlow, Kutta);
		wif_algo::calculation_results_c calculate_flow3 = wif_algo::calculate_flow(myAirfoil2, myFlow, Kutta);
		wif_algo::calculation_results_c calculate_flow4 = wif_algo::calculate_flow(myAirfoil3, myFlow, Kutta);

		std::cout<< calculate_flow2.c_l <<"   "<< calculate_flow3.c_l<< "   " << calculate_flow4.c_l << std::endl;





	return 0;
}
