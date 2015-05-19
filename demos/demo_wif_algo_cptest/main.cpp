#include <iostream>
#include <cmath>
#include <fstream>

#include <wif_core/wif_core.hpp>
#include <wif_algo/wif_algo.hpp>
#include <wif_viz/wif_viz.hpp>

/*
#include <TCanvas.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLegend.h>
#include <TGraph.h>
*/

int main()
{
	double pi = 3.1415;
	double radius = 1;
	unsigned int num_lines = 100;
	wif_core::vector_2d_c midpoint(0, 0);
	wif_core::airfoil_c myAirfoil(midpoint, radius, num_lines);


	std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);
	bool Kutta = 0;
	wif_algo::calculation_results_c calculate_circle = wif_algo::calculate_flow(myAirfoil, myFlow, Kutta);


	std::vector<wif_core::line_2d_c> mylines = myAirfoil.get_lines();
	std::vector<wif_core::vector_2d_c> centers(num_lines);
	std::vector<double> angles(num_lines);
	std::vector<double> c_p_boven;
	std::vector<double> x_as_boven;
	std::vector<double> c_p_onder;
	std::vector<double> x_as_onder;

	std::cout << "c_p algoritme" << "\t" << "c_p polair theoretisch" << "\t" << "c_p cartesisch theoretisch" << std::endl;

	for(unsigned int i = 0; i < num_lines; i++)
	{
		wif_core::line_2d_c temp_line = mylines[i];
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
			std::cout << i << "  Not in any of these categories" << std::endl;
		}

		std::cout << calculate_circle.c_p[i] << "\t" << 1 - 4 * pow(sin(angles[i]), 2) << "\t" << 1 - 4 * pow(centers[i].y / radius, 2) << std::endl;
	}


	return 0;
}
