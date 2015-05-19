#include <iostream>
#include <cmath>
#include <fstream>

#include <wif_core/wif_core.hpp>
#include <wif_algo/wif_algo.hpp>
#include <wif_viz/wif_viz.hpp>

#include <ctime>
int main()
{

	double pi = 3.1415;
	/*double radius = 2;


	unsigned int n=20;
	std::clock_t start;
	std::vector<double> times(n);
	std::vector<double> fout(n);
	std::vector<double> som_sigma_l(n);
	std::vector<double> num_lines_all(n);

	for(int q = 1; q < n+1; q++)
	{
		start = std::clock();

		unsigned int num_lines = 10 * q;
		wif_core::vector_2d_c midpoint(0, 0);
		wif_core::airfoil_c myAirfoil(midpoint, radius, num_lines);
		std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);
		bool Kutta = 0;
		wif_algo::calculation_results_c calculate_flow2 = wif_algo::calculate_flow(myAirfoil, myFlow, Kutta);
		double gamma = 0;

		std::vector<wif_core::line_2d_c> mylines = myAirfoil.get_lines();
		std::vector<wif_core::vector_2d_c> centers(num_lines);
		std::vector<double> angles(num_lines);
	    double fout_accum=0;
	    std::vector<double> x_plot(num_lines);
	    std::vector<double> cp_theory(num_lines);

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

			cp_theory[i]=1- 4 * pow((centers[i].y / radius), 2);
	        fout_accum+=pow(1 - 4 * pow((centers[i].y / radius), 2)-calculate_flow2.c_p[i],2);
	        x_plot[i]=i+1;

		}

	    fout[q-1]= log10((fout_accum/num_lines));
		times[q-1]=(std::clock() - start) / (double)CLOCKS_PER_SEC ;
		som_sigma_l[q-1]=(calculate_flow2.closed_body_check)*1000;
		num_lines_all[q-1]=num_lines;

	    std::cout<< fout[q-1] <<std::endl;

		std::vector<std::vector<double>> y_plot;
		y_plot.push_back(calculate_flow2.c_p);
		y_plot.push_back(cp_theory);

	    std::vector<std::string> legende_cp;
	    legende_cp.push_back("Numeriek bepaalde waarden");
	    legende_cp.push_back("Theoretische waarden");
	    std::string title_cp ="Cp: "+ std::to_string(num_lines)+ "panelen gebruikt";
	    std::string filename_cp= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/plot_cp "+std::to_string(num_lines)+".png";

	    std::shared_ptr<wif_viz::visualization_c> myRoot_1 = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
	    myRoot_1->plotVectors(y_plot,x_plot, legende_cp, filename_cp, "Index paneel", "Cp",title_cp);

	}


	//Dummy myFlow en midpoint:
	wif_core::vector_2d_c midpoint(0, 0);
	wif_core::airfoil_c myAirfoil(midpoint, 2, 10);
	std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);

	//Plots:
	std::vector<std::string> legende_error;
	legende_error.push_back("Fout");
	std::string title_error ="Fout op cp als functie van het aantal panelen";
	std::string filename_error= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/fout op cp.png";
	std::shared_ptr<wif_viz::visualization_c> myRoot_2 = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
	std::vector<std::vector<double>> y_plot_fout;
	y_plot_fout.push_back(fout);
	myRoot_2 ->plotVectors(y_plot_fout,num_lines_all, legende_error, filename_error, "Aantal panelen", "log10(Fout op Cp)",title_error);

	std::vector<std::string> legende_time;
	legende_time.push_back("Tijdsduur");
	std::string title_time ="Tijdsduur als functie van het aantal panelen";
	std::string filename_time= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/tijd cp.png";
	std::shared_ptr<wif_viz::visualization_c> myRoot_3 = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
	std::vector<std::vector<double>> y_plot_time;
	y_plot_time.push_back(times);
	myRoot_3 ->plotVectors(y_plot_time,num_lines_all,legende_time, filename_time, "Aantal panelen", "Tijd (s)",title_time);

	std::vector<std::string> legende_sigma;
	legende_sigma.push_back("Som sigma * l");
	std::string title_sigma ="Som sigma * l als functie van het aantal panelen";
	std::string filename_sigma= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/som sigma l cp.png";
	std::shared_ptr<wif_viz::visualization_c> myRoot_4 = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
	std::vector<std::vector<double>> y_plot_sigma;
	y_plot_sigma.push_back(som_sigma_l);
	myRoot_4 ->plotVectors(y_plot_sigma,num_lines_all,legende_sigma, filename_sigma, "Aantal panelen", "Som sigma *l (10^3)",title_sigma);
	*/


	///////////////////////////////////////////////////
	double radius_airfoils = 1;
	unsigned int num_lines_airfoils = 200;
	double angle_in_degrees = 0;

	bool Kutta_airfoils = 1;
	std::string path_airfoils = "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/";

	wif_core::vector_2d_c midpoint_airfoils(0, 0);
	std::shared_ptr<wif_core::uniform_flow_c> UFlow = std::make_shared<wif_core::uniform_flow_c>((angle_in_degrees / 180) * pi, 1);

	std::cout << "radius = " << radius_airfoils << std::endl;
	std::cout << "angle of attack in degrees = " << angle_in_degrees << std::endl;

	//Verschillende airfoils:
	//A18 http://airfoiltools.com/airfoil/details?airfoil=a18-il
	std::cout << "A18 (original)" << std::endl;
	wif_core::airfoil_c Airfoil_1(path_airfoils + "A18.dat");
	Airfoil_1 = Airfoil_1.get_circle_projection(num_lines_airfoils, midpoint_airfoils, radius_airfoils, 0.0001);
	Airfoil_1 = Airfoil_1.closed_merge();
	wif_algo::calculation_results_c results_Airfoil_1 = wif_algo::calculate_flow(Airfoil_1, UFlow, Kutta_airfoils);
	std::cout << "c_l=" << results_Airfoil_1.c_l << std::endl;

	//NACA 0012 http://airfoiltools.com/airfoil/details?airfoil=n0012-il
	std::cout << "NACA 0012" << std::endl;
	wif_core::airfoil_c Airfoil_2(path_airfoils + "NACA 0012.dat");
	Airfoil_2 = Airfoil_2.get_circle_projection(num_lines_airfoils, midpoint_airfoils, radius_airfoils, 0.0001);
	Airfoil_2 = Airfoil_2.closed_merge();
	wif_algo::calculation_results_c results_Airfoil_2 = wif_algo::calculate_flow(Airfoil_2, UFlow, Kutta_airfoils);
	std::cout << "c_l=" << results_Airfoil_2.c_l << std::endl;


	//DAVIS BASIC B-24 WING http://airfoiltools.com/plotter/index?airfoil=davis_corrected-il
	std::cout << "DAVIS BASIC B-24 WING" << std::endl;
	wif_core::airfoil_c Airfoil_3(path_airfoils + "DAVIS BASIC B-24 WING.dat");
	Airfoil_3 = Airfoil_3.get_circle_projection(num_lines_airfoils, midpoint_airfoils, radius_airfoils, 0.0001);
	Airfoil_3 = Airfoil_3.closed_merge();
	wif_algo::calculation_results_c results_Airfoil_3 = wif_algo::calculate_flow(Airfoil_3, UFlow, Kutta_airfoils);
	std::cout << "c_l=" << results_Airfoil_3.c_l << std::endl;

	return 0;
}
