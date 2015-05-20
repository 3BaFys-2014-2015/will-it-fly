#include <iostream>

#include <wif_core/wif_core.hpp>
#include <wif_viz/wif_viz.hpp>
#include <wif_algo/wif_algo.hpp>

int main()
{
	std::shared_ptr<wif_core::flow_c> unifl = std::make_shared<wif_core::uniform_flow_c>(0.0, 0.10);

	std::cout << "1.1" << std::endl;
	std::shared_ptr<wif_core::flow_accumulate_c> ff = std::make_shared<wif_core::flow_accumulate_c>();//unifl);

	/*
	#if 1
		wif_core::airfoil_c airfoil({0, 0}, 1, 30);
	#else
		wif_core::airfoil_c airfoil("../../wif_core/airfoils/selig.dat");
	#endif*/

	wif_core::airfoil_c airfoil({0, 0}, 1, 30);
	std::cout << "2" << std::endl;
	//std::shared_ptr<wif_core::flow_accumulate_c> ff = std::make_shared<wif_core::flow_accumulate_c>();
	std::cout << "2.1" << std::endl;
	ff->add_source_sheets(std::vector<double>(airfoil.get_lines().size(), 1), airfoil);
	std::cout << "3" << std::endl;

	wif_core::vector_2d_c min, max, bins;
	min.x = -2;
	min.y = -2;
	max.x = 2;
	max.y = 2;
	bins.x = 20;
	bins.y = 20;
	//int binsx = 20, binsy = 20;

	std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(unifl, min, max);
	//vizy->set_velocityarrows(bins);
	vizy->set_psi_bins(bins);
	vizy->set_phi_bins(bins);
	vizy->set_velocity_bins(bins);
	std::cout << "4" << std::endl;
	vizy->draw("circle");
	std::cout << "5" << std::endl;

	return 0;

}
