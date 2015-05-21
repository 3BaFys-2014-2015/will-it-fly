#include <iostream>

#include <wif_core/wif_core.hpp>
#include <wif_algo/wif_algo.hpp>
#include <wif_viz/wif_viz.hpp>

#include <algorithm>

void test_airfoil(wif_core::airfoil_c & foil)
{
	if(!foil.is_valid())
	{
		return;
	}

	if(!foil.check_lengths())
	{
		std::cout << "LENGTH ERROR" << std::endl;
		return;
	}

	std::shared_ptr<wif_core::uniform_flow_c> flow = std::make_shared<wif_core::uniform_flow_c>(0.0 * M_PI / 180.0, 1.0);

	auto result = wif_algo::calculate_flow(foil, flow, true, 0.0);

	std::cout << "Calculated flow from: " << foil.get_name() << std::endl;
	std::cout << "In a uniform flow with strength " << flow->get_strength() << " and AoA " << flow->get_angle() * 180.0 / M_PI << std::endl;
	std::cout << std::endl;
	std::cout << "C_l = " << result.c_l << std::endl;
	std::cout << "CBC = " << result.closed_body_check << std::endl;
	std::cout << "First v_t: " << result.v_t.front() << std::endl;
	std::cout << "Last  v_t: " << result.v_t.back() << std::endl;
	std::cout << std::endl;

	{
		uint32_t i = foil.get_index_of_first_lower_panel();
		std::cout << i << std::endl;

		std::vector<std::vector<double>> data;

		std::vector<double> c_p(result.c_p.begin(), result.c_p.begin() + i);
		std::reverse(c_p.begin(), c_p.end());

		data.push_back(c_p);

		std::vector<std::string> Legend(1);
		Legend[0] = "CP Bovenkant";
		//Legend[1] = "Onderkant";
		//Legend[2] = "Verschil";

		foil.get_upper_panels_x();
		foil.get_lower_panels_x();

		std::vector<double> x_as;

		for(int j = 0; j < i; j++)
		{
			std::cout << foil.get_lines()[j].get_center_point().x << std::endl;

			x_as.push_back(foil.get_lines()[j].get_center_point().x);
		}

		std::reverse(x_as.begin(), x_as.end());

		std::shared_ptr<wif_viz::visualization_c> root = wif_viz::create_visualization_root(0, {0, 0}, {0, 0});
		root->plotVectors(data, x_as, Legend, "lel.png", "x", "cp", "Aantal panelen = 100, Alpha = 45, Met Kutta");
	}

	{
		std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(result.flow, { -0.5, -1.0}, {1.5, 1});
		vizy->get_psi_field().bins = {201, 201};
		vizy->get_psi_field().style = (wif_viz::E_SCALAR_DRAW_STYLE)(wif_viz::ESDS_GRADIENT | wif_viz::ESDS_CONTOURS);
		vizy->get_phi_field().bins = {201, 201};
		vizy->get_phi_field().style = wif_viz::ESDS_CONTOURS;

		vizy->set_airfoil(&foil);

		vizy->draw("");
	}

	{
		std::shared_ptr<wif_viz::visualization_c> vizy = wif_viz::create_visualization_vtk(result.flow, { -0.5, -0.5}, {1.5, 0.5});
		vizy->get_v_field().bins = {101, 101};
		vizy->get_v_field().style = wif_viz::EVDS_STREAMLINES;
		vizy->get_v_field().arrow_scale = 0.001;
		vizy->get_v_field().streamline_resolution = 50;
		vizy->get_v_field().streamline_seeds = { -0.5, 0.5, -0.5, -0.5};

		vizy->set_airfoil(&foil);

		vizy->draw("41584864");
		std::cout << "DFSDF" << std::endl;
	}
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
	std::string filename = "../../../coord_seligFmt/mh70.dat";
	//std::string filename = "../../wif_core/airfoils/lednicer.dat";
	//std::string filename = "../../wif_core/airfoils/n0012-il.dat";
	//std::string filename = "../../wif_core/airfoils/selig.dat";
#endif

	if(true)
	{
		wif_core::airfoil_c a = wif_core::airfoil_c(filename);

		if(!a.is_valid())
		{
			std::cout << "Error loading datafile" << std::endl;
			return 0;
		}


		//.closed_intersect(0)/*)*/;

		for(const auto & p : a.get_points())
		{
			std::cout << p << std::endl;
		}

		std::cout << std::endl << std::endl;

		a = a.closed_intersect(0.0).get_circle_projection(40);

		for(const auto & p : a.get_points())
		{
			std::cout << p << std::endl;
		}

		test_airfoil(a);
	}
	else
	{
		wif_core::airfoil_c x({0.5, 0.0}, 0.5, 10);

		test_airfoil(x);
	}


	return 0;
}
