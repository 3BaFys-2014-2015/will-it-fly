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
	wif_core::airfoil_c myCircle(midpoint, radius, num_lines);

	wif_core::airfoil_c myNaca("wif_core/airfoils/n0012-il.dat");
	myNaca = myNaca.get_circle_projection(num_lines, midpoint, radius);
	myNaca = myNaca.closed_merge();

	std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);
	bool Kutta = 0;
	wif_algo::calculation_results_c calculate_circle = wif_algo::calculate_flow(myCircle, myFlow, Kutta);
	wif_algo::calculation_results_c calculate_naca = wif_algo::calculate_flow(myNaca, myFlow, Kutta);

	std::cout << "Closed body check for circle: " << calculate_circle.closed_body_check << std::endl;
	std::cout << "Closed body check for naca0012: " << calculate_naca.closed_body_check << std::endl;


	return 0;
}
