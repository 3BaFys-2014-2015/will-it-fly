#include <iostream>

#include <wif_core/wif_core.hpp>
#include <wif_algo/wif_algo.hpp>
#include <wif_viz/wif_viz.hpp>
#include <fstream>

int main(int argc, char * argv[])
{

	/*system("cd /home/uauser/Shared/coord_seligFmt.tar/coord_seligFmt");
	system("ls > files.txt");
	std::ifstream filelist("/home/uauser/Shared/coord_seligFmt.tar/coord_seligFmt/files.txt");*/

	if(argc >= 1)
	{
		std::string filename1;
		filename1 = argv[1];
		std::cout << filename1 << std::endl;
		std::cout << filename1 << std::endl;
		wif_core::airfoil_c airfoil1(filename1);
		std::cout << airfoil1;
		std::cout << airfoil1.closed_merge();
		std::cout << airfoil1.closed_intersect(0);
		std::cout << airfoil1.get_circle_projection(8, wif_core::vector_2d_c(0.5, 0), 0.5, 0.0001);

		for(wif_core::line_2d_c line : airfoil1.get_lines())
		{
			std::cout << line.get_length() << ",";
		}

		std::cout << std::endl;

		for(wif_core::line_2d_c line : airfoil1.get_circle_projection(8, wif_core::vector_2d_c(0.5, 0), 0.5, 0.0001).get_lines())
		{
			std::cout << line.get_length() << ",";
		}
	}

	//break;
	//std::cout << airfoil1.get_circle_projection(9, wif_core::vector_2d_c(0.5, 0.0) , 0.5);


	/* Deze code geeft een gigantische memory leak. Het crasht mijn vb... Fixt dit aub
	Gebruik ook geen absolute paden...
	kan niet anders, want een relatief pad kan elleen verterekken in de build/bin directory.

	std::string filename1 = "/home/uauser/eindopdracht/will-it-fly/wif_core/airfoils/lednicer.dat";
	wif_core::airfoil_c airfoil1(filename1);
	std::cout << airfoil1;
	std::cout << airfoil1.get_circle_projection(9, wif_core::vector_2d_c(0.5, 0.0) , 0.5);
	std::string filename2 = "/home/uauser/eindopdracht/will-it-fly/wif_core/airfoils/selig.dat";
	wif_core::airfoil_c airfoil2(filename2);
	std::cout << airfoil2;
	std::cout << airfoil2.get_circle_projection(9, wif_core::vector_2d_c(0.5, 0.0), 0.5);

	std::vector<wif_core::vector_2d_c> test;
	test.push_back(wif_core::vector_2d_c(1.0, 0.02));
	test.push_back(wif_core::vector_2d_c(0.0, 0.02));
	wif_core::airfoil_c airfoil(test, "flat");
	std::cout << airfoil;
	std::cout << airfoil.get_circle_projection(10, wif_core::vector_2d_c(0.5, 0.0), 0.5);

	wif_core::airfoil_c airfoil3(wif_core::vector_2d_c(0.5, 0.0), 0.5, 8);
	std::cout << airfoil3;*/


	return 0; // Komaan mensen? return code vergeten????
}
