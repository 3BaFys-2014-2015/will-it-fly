#include <iostream>

#include <wif_core/wif_core.hpp>
#include <wif_algo/wif_algo.hpp>
#include <wif_viz/wif_viz.hpp>

int main()
{
	std::cout << wif_core::get_version() << std::endl;
	std::cout << wif_algo::get_version() << std::endl;
	std::cout << wif_viz::get_version()  << std::endl;

	std::shared_ptr<wif_core::flow_c> unifl = std::make_shared<wif_core::uniform_flow_c>();

	std::shared_ptr<wif_core::flow_accumulate_c> flow = std::make_shared<wif_core::flow_accumulate_c>();//unifl);
	std::shared_ptr<wif_core::source_sheet_c> ss = std::make_shared<wif_core::source_sheet_c>();
	
	flow->add_flow(ss);//wif_core::line_2d_c(-1,-1,1,1), 1));
	//flow->add_flow(unifl);
	
	std::cout << ss->get_psi({1, 1});

	wif_core::vector_2d_c min, max, bins;
	min.x = -2;
	min.y = -2;
	max.x = 2;
	max.y = 2;
	bins.x = 100;
	bins.y = 100;
	//int binsx = 20, binsy = 20;

	std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(flow, min, max);

	//vizy->set_velocityarrows(bins);
	vizy->set_psi_bins(bins);
	vizy->set_phi_bins(bins);

	vizy->draw("test.png");

}
