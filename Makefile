# Makefile
# Compilatie objecten (*.o) gaan naar ./obj/, binaries gaan naar ./bin/, library headers (*.hpp) gaan naar ./include/ en libraries (*.so) gaan naar ./lib/

SRC_PATH=$(wildcard ./will_it_fly/*.cpp)
SRC_WIF=$(notdir $(SRC_PATH))
OBJ_WIF=$(SRC_WIF:.cpp=.o)
OBJ_PATH=$(addprefix ./build/obj/wif/, $(OBJ_WIF))

PATH:=/home/uauser/Software/root-v5.34.25/builddir/bin:$(PATH)
LD_LIBRARY_PATH:=/home/uauser/Software/root-v5.34.25/builddir/lib:$(LD_LIBRARY_PATH)

export CXX=g++ -std=c++11
export CXX_FLAGS=-Wall
export BUILD_SUFFIX=

INCLUDE=-I./build/include
LIBDIR=-L./build/lib
LIBS=-lwif_core -lwif_algo -lwif_viz
EXECUTABLE=./build/bin/wif

INSTALL_PATH=/usr/

all : demos wif
	echo "Building all"

debug : export CXX_FLAGS=-Wall -g
debug : export BUILD_SUFFIX=_d
debug : all

optimized : export CXX_FLAGS=-Wall -O3
optimized : export BUILD_SUFFIX=_o
optimized : all

doxygen : builddir
	doxygen
	ln -sLf ./html/index.html ./build/doc/index.html

libs : wifcore wifalgo wifviz

# Makefile demos

builddir :
	mkdir -p ./build
	mkdir -p ./build/doc
	mkdir -p ./build/obj
	mkdir -p ./build/include
	mkdir -p ./build/bin
	mkdir -p ./build/lib
	mkdir -p ./build/exec

demo_% : export DEMO_NAME = $@
demo_% : builddir libs
	mkdir -p ./build/obj/demos
	mkdir -p ./build/obj/demos/$@
	$(MAKE) -C ./demos/$@ demo

demos : viz_demos algo_demos core_demos team_demos

team_demos : team_viz_demos team_core_demos team_algo_demos

# De demos van de viz library
viz_demos : demo_vtk

# De demos van de algo library
algo_demos :

# De demos van de core library
core_demos : demo_core_vector_line demo_core_flow_elements demo_core_airfoil demo_solve_234

# De demos van de teams, die eigenlijk niks met onze libraries te maken hebben.
# Hieronder demo's toevoegen als dependancies van de respectievelijke targets!
# Gebruik dezelfde naam als uw directory!

team_viz_demos : vtk_compilatie_demos demo_vtk demo_vtk_large demo_root

team_core_demos :

team_algo_demos : root_algo_demos gsl_demos demo_algo_diagonalizationspeed demo_algo_integrationspeed demo_wif_algo demo_wif_algo_amy  demo_wif_algo_closedbodycheck demo_wif_algo_cptest

vtk_compilatie_demos : demo_vtk_comp_dolf demo_vtk_comp_kristof demo_vtk_comp_merel demo_vtk_comp_senne demo_vtk_merel demo_vtk_wtf

gsl_demos : demo_gsl_integral demo_gsl_diag demo_gsl_andy

root_algo_demos : demo_root_algo_diag demo_root_algo_integral

wif : libs
	mkdir -p ./build/obj/wif
	$(MAKE) -C ./will_it_fly willitfly

wifcore : builddir
	mkdir -p ./build/obj/wif_core
	mkdir -p ./build/include/wif_core
	cat ./wif_core/install_headers.txt | xargs -I % cp -u ./wif_core/% ./build/include/wif_core
	cat ./wif_core/install_headers.txt | xargs -I % chmod 755 ./build/include/wif_core/%
	$(MAKE) -C ./wif_core core

wifalgo : wif_core
	mkdir -p ./build/obj/wif_algo
	mkdir -p ./build/include/wif_algo
	cat ./wif_algo/install_headers.txt | xargs -I % cp -u ./wif_algo/% ./build/include/wif_algo
	cat ./wif_algo/install_headers.txt | xargs -I % chmod 755 ./build/include/wif_algo/%
	$(MAKE) -C ./wif_algo algo

wifviz : wif_core
	mkdir -p ./build/obj/wif_viz
	mkdir -p ./build/include/wif_viz
	cat ./wif_viz/install_headers.txt | xargs -I % cp -u ./wif_viz/% ./build/include/wif_viz
	cat ./wif_viz/install_headers.txt | xargs -I % chmod 755 ./build/include/wif_viz/%
	$(MAKE) -C ./wif_viz viz

clean :
	echo "Cleaning build directory"
	rm -rf ./build

install:
	echo "Installing... Please wait."
	cp -ru ./build/lib/* $(INSTALL_PATH)/lib/
	chmod 755 $(INSTALL_PATH)/lib/*
	cp -ru ./build/include/* $(INSTALL_PATH)/include/
	cp -ru ./build/bin/wif $(INSTALL_PATH)/bin/
	chmod 755 $(INSTALL_PATH)/bin/wif
	#$(MAKE) clean
	echo "Done"

uninstall:
	echo "Uninstalling... Please wait."
	rm -f $(INSTALL_PATH)/lib/libwif_core.so
	rm -f $(INSTALL_PATH)/lib/libwif_algo.so
	rm -f $(INSTALL_PATH)/lib/libwif_viz.so
	rm -rf $(INSTALL_PATH)/include/wif_core
	rm -rf $(INSTALL_PATH)/include/wif_algo
	rm -rf $(INSTALL_PATH)/include/wif_viz
	rm -f $(INSTALL_PATH)/bin/wif
	echo "Done"

