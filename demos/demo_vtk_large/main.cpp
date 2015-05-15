#include <iostream>

#include <wif_core/wif_core.hpp>
#include <wif_algo/wif_algo.hpp>
#include <wif_viz/wif_viz.hpp>

#include <thread>

void print_all_fields(const std::string & name, std::shared_ptr<wif_core::flow_c>  flow,
                      const wif_core::vector_2d_c                 &                min,
                      const wif_core::vector_2d_c                 &                max,
                      const wif_core::vector_2d_c                 &                binning)
{
	std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(flow, min, max);

	vizy->set_psi_bins(binning);
	vizy->set_phi_bins(binning);
	vizy->set_velocity_bins(binning);

	vizy->draw_ivo(name);
}

void print_psi(const std::string        &        name,
               std::shared_ptr<wif_core::flow_c> flow,
               const wif_core::vector_2d_c   &   min,
               const wif_core::vector_2d_c   &   max,
               const wif_core::vector_2d_c   &   binning)
{
	std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(flow, min, max);
	vizy->set_clip_range(-0.1, 0.1);
	vizy->set_psi_bins(binning);

	vizy->draw_ivo("");
}

//

//



//

void print_phi(const std::string        &        name,
               std::shared_ptr<wif_core::flow_c> flow,
               const wif_core::vector_2d_c   &   min,
               const wif_core::vector_2d_c   &   max,
               const wif_core::vector_2d_c   &   binning)
{
	std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(flow, min, max);
	vizy->set_clip_range(-0.1, 0.1);
	vizy->set_phi_bins(binning);

	vizy->draw_ivo("");
}


void print_velocity(const std::string        &        name,
                    std::shared_ptr<wif_core::flow_c> flow,
                    const wif_core::vector_2d_c   &   min,
                    const wif_core::vector_2d_c   &   max,
                    const wif_core::vector_2d_c   &   binning)
{
	std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(flow, min, max);

	vizy->set_velocity_bins(binning);


	vizy->draw_ivo("");
}



void visualize_all(bool                               screen,
                   const std::string         &        name,
                   std::shared_ptr<wif_core::flow_c>  flow,
                   const wif_core::vector_2d_c    &   min,
                   const wif_core::vector_2d_c    &   max,
                   const wif_core::vector_2d_c    &   binning)
{
	if(screen)
	{
		print_psi(name, flow, min, max, binning);
		print_phi(name, flow, min, max, binning);
		print_velocity(name, flow, min, max, binning);
	}
	else
	{
		print_all_fields(name, flow, min, max, binning);
	}
}


std::string get_name(bool screen, const std::string & filename)
{
	return screen ? "" : filename;
}

void test_airfoil(wif_core::airfoil_c & foil)
{
	if(!foil.is_valid())
	{
		return;
	}

	std::shared_ptr<wif_core::uniform_flow_c> flow = std::make_shared<wif_core::uniform_flow_c>(0.0 * M_PI / 180.0, 1.0);

	auto result = wif_algo::calculate_flow(foil, flow, false, 0.0);

	{
		std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(result.flow, { -0.5, -1.0}, {1.5, 1});
		vizy->set_psi_bins({101, 101});
		vizy->set_phi_bins({101, 101});
		vizy->set_contours(20);
		vizy->set_airfoil(&foil);

		vizy->draw_ivo("");
	}

	{
		std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(result.flow, { -0.5, -1.0}, {1.5, 1});
		vizy->set_velocity_bins({101, 101});
		vizy->set_streamline_resolution(100);
		vizy->set_airfoil(&foil);

		vizy->draw_ivo("");
	}
}



void test_uniflow(bool screen)
{
	std::shared_ptr<wif_core::uniform_flow_c> flow = std::make_shared<wif_core::uniform_flow_c>(M_PI / 4.0, 1.0);

	std::cout << flow->get_angle() / M_PI << std::endl;

	{
		std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(flow, { -1, -1}, {1, 1});
		vizy->set_psi_bins({101, 101});
		vizy->set_contours(20);

		vizy->draw_ivo("");
	}

	{
		std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(flow, { -1, -1}, {1, 1});
		vizy->set_phi_bins({101, 101});

		vizy->draw_ivo("");
	}

	//visualize_all(screen, "test-uniflow", flow, { -2, -2}, {2, 2}, {31, 31});

	//std::cout << flow->get_psi({ -2.0, -2.0}) << std::endl;
}


void test_circle(bool screen)
{
	wif_core::airfoil_c airfoil({0, 0}, 1, 30);

	std::shared_ptr<wif_core::flow_accumulate_c> flow = std::make_shared<wif_core::flow_accumulate_c>();
	flow->add_source_sheets(std::vector<double>(airfoil.get_lines().size(), 1), airfoil);

	visualize_all(screen, "test-circle", flow, { -2, -2}, {2, 2}, {101, 101});
}


void test_circle_flow(bool screen)
{
	wif_core::airfoil_c airfoil({0, 0}, 1, 4);

	std::shared_ptr<wif_core::flow_c> uniflow = std::make_shared<wif_core::uniform_flow_c>(0.0, -1.0);
	std::shared_ptr<wif_core::flow_accumulate_c> flow = std::make_shared<wif_core::flow_accumulate_c>();
	flow->add_source_sheets(std::vector<double>(airfoil.get_lines().size(), 1.0), airfoil);

	visualize_all(screen, "test-circle-flow", flow, { -2, -2}, {2, 2}, {101, 101});
}

void test_sheet(bool screen)
{
	std::shared_ptr<wif_core::flow_c> uniflow = std::make_shared<wif_core::uniform_flow_c>(0.0, 0.0);
	std::shared_ptr<wif_core::flow_accumulate_c> flow = std::make_shared<wif_core::flow_accumulate_c>();

	wif_core::line_2d_c l({0, -1.0}, {0, 1.0});

	flow->add_flow(std::make_shared<wif_core::vortex_sheet_c>(l, 1.0));
	//flow->add_source_sheets(std::vector<double>(airfoil.get_lines().size(), 1.0), airfoil);

	visualize_all(screen, "test-circle-flow", flow, { -2, -2}, {2, 2}, {101, 101});
}

using namespace wif_core;
using namespace wif_viz;

void test_source(bool screen)
{
	vector_2d_c min(-2.0, -2.0);
	vector_2d_c max(2.0, 2.0);

	auto source_sink = std::make_shared<source_sink_c>(vector_2d_c(0.0, 0.0), 1.0);

	{
		auto viz = create_visualization_vtk(source_sink, min, max);

		viz->set_clip_range(-0.1, 0.1);
		viz->set_psi_bins({100, 100});
		viz->draw();
	}

	{
		auto viz = create_visualization_vtk(source_sink, min, max);

		viz->set_clip_range(-1.0, 1.0);
		viz->set_phi_bins({100, 100});
		viz->draw();
	}


	visualize_all(screen, "test-source", source_sink, { -2, -2}, {2, 2}, {101, 101});

	auto source_sink_2 = std::make_shared<source_sink_c>(vector_2d_c(1.0, 1.0), 1.0);

	visualize_all(screen, "test-source", source_sink_2, { -2, -2}, {2, 2}, {101, 101});

	auto both = std::make_shared<flow_accumulate_c>(source_sink);

	both->add_flow(source_sink_2);

	visualize_all(screen, "test-source", both, { -2, -2}, {2, 2}, {101, 101});
}

void test_airfoil()
{

#if 1
	wif_core::airfoil_c airfoil("../../wif_core/airfoils/selig.dat");
	wif_core::airfoil_c n_airfoil = airfoil.closed_intersect(0);

	auto lines = n_airfoil.get_points();

	for(const auto & l : lines)
	{
		std::cout << l << std::endl;
	}

	//std::cout << n_airfoil << std::endl;*/
#else
	wif_core::airfoil_c airfoil({0, 0}, 1, 100);
	//wif_core::airfoil_c n_airfoil = airfoil.get_circle_projection(100, {0.0, 0.0}, 2.0);

	wif_core::airfoil_c n_airfoil = airfoil;
#endif // 1

	if(!airfoil.is_valid())
	{
		return;
	}

	std::shared_ptr<wif_core::uniform_flow_c> uniflow = std::make_shared<wif_core::uniform_flow_c>(0.2, 1.0);
	std::shared_ptr<wif_core::flow_accumulate_c> flow = std::make_shared<wif_core::flow_accumulate_c>(uniflow);

	auto f = wif_algo::calculate_flow(n_airfoil, uniflow, false, 0.0);

	//flow->add_flow(f);
	//flow->add_source_sheets(std::vector<double>(airfoil.get_lines().size(), 1), airfoil);

	std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(f.flow, { -2, -3}, {2, 3});
	//vizy->set_clip_range(-1, 1);
	vizy->set_velocity_bins({51, 51});
	vizy->set_airfoil(&n_airfoil);
	vizy->set_contours(200);
	vizy->set_streamline_resolution(200);

	vizy->draw_ivo("");

	//visualize_all(screen, "test-circle", flow, { -2, -2}, {2, 2}, {101, 101});
}

void tests()
{
	bool screen = true;

	//test_uniflow(screen);
	//test_circle(screen);
	//test_circle_flow(screen);
	//test_sheet(screen);
	//test_source(screen);

	wif_core::airfoil_c a = wif_core::airfoil_c("../../wif_core/airfoils/selig.dat").closed_intersect(0)/*.get_circle_projection(10, {0.5, 0.0}, 0.5)*/;

	test_airfoil(a);
}


int main(int argc, char ** argv)
{
#if 0

	if(argc != 2)
	{
		std::cout << "No filename" << std::endl;
		return 0;
	}

	std::string filename = argv[1];
#else
	//std::string filename = "../../../coord_seligFmt/mh114.dat";
	//std::string filename = "../../wif_core/airfoils/lednicer.dat";
	std::string filename = "../../wif_core/airfoils/n0012-il.dat";
#endif
	wif_core::airfoil_c a = wif_core::airfoil_c(filename);

	if(!a.is_valid())
	{
		std::cout << "Error loading datafile" << std::endl;
		return 0;
	}

	//.closed_intersect(0)/*.get_circle_projection(10, {0.5, 0.0}, 0.5)*/;

	a = a.closed_intersect(0);

	test_airfoil(a);

	return 0;
}
