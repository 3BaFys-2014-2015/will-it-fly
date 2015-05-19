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
	double radius = 2;
	unsigned int num_lines = 100;
	wif_core::vector_2d_c midpoint(0, 0);
	wif_core::airfoil_c myAirfoil(midpoint, radius, num_lines);
	std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);
	bool Kutta = 1;

	wif_algo::calculation_results_c calculate_flow2 = wif_algo::calculate_flow(myAirfoil, myFlow, Kutta);
	//double gamma = 0;

	std::vector<wif_core::line_2d_c> mylines = myAirfoil.get_lines();
	std::vector<wif_core::vector_2d_c> centers(num_lines);
	std::cout << "y" << "   " << "cp calc" << "   " << "cp theor pol" << "   " << "cp theor cart" << std::endl;
	std::vector<double> angles(num_lines);
	//std::ofstream myfile;
	//myfile.open("data2.txt");

	for(unsigned int i = 0; i < num_lines; i++)
	{
		wif_core::line_2d_c temp_line = mylines[i];
		centers[i] = temp_line.get_center_point();

		if(centers[i].y > 0)
		{
			angles[i] = atan2(centers[i].y, centers[i].x);
		}
		else
		{
			angles[i] = atan2(centers[i].y, centers[i].x) + 2 * pi;
		}

		std::cout << centers[i].y << "   " << calculate_flow2.c_p[i] << "   " << 1 - 4 * pow(sin(angles[i]), 2) << "   " << 1 - 4 * pow((centers[i].y / radius), 2)   << std::endl;
	}

	//myfile.close();


	////////////////////////Plot////////////////////////


	/*
	TCanvas * canvas = new TCanvas("canvas", "Name", 0, 0, 800, 800);
	canvas->SetRightMargin(0.15);
	canvas->SetTopMargin(0.15);
	canvas->SetFillColor(18);

	std::vector<double> y(num_lines);

	for(int i = 0; i < num_lines; i++)
	{
		y[i] = centers[i].y;
	}

	TGraph * cp = new TGraph(100, y, calculate_flow2.c_p);

	gStyle->SetTitleW(0.99);
	gStyle->SetTitleH(0.08);

	cp->Draw("AC*");
	cp->GetXaxis()->SetTitle("y");
	cp->GetYaxis()->SetTitle("Drukcoefficient");
	cp->SetTitle("Drukcoefficient, r=a");


	std::vector<double> druk(num_lines);

	for(int i = 0; i < 100; i++)
	{
		druk[i] = 1 - 4 * pow(centers[i].y / radius, 2);
	}

	TGraph * cp_t = new TGraph(100, centers.y, druk);
	cp_t->SetLineColor(3);
	cp_t->Draw("SAME");


	TLegend * leg = new TLegend(0.12, 0.14, 0.3, 0.33);
	leg->AddEntry("cp", "Drukcoefficient", "l");
	leg->AddEntry("druk2", "1-4(y/R)^2", "l");
	leg->Draw("SAME");

	cc->Update();
	cc->Print("drukcoefficient.pdf");
	*/
	return 0;
}
