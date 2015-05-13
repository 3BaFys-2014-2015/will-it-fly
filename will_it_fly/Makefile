# Makefile will-it-fly

SRC_WIF=$(wildcard *.cpp)
OBJ_WIF=$(SRC_WIF:.cpp=$(BUILD_SUFFIX).o)
OBJ_PATH=$(addprefix ../build/obj/wif/,  $(OBJ_WIF))

CXX=g++ -std=c++11

INCLUDE=-I../build/include -I/usr/include/gsl
LIBDIR=-L../build/lib -L/usr/lib/ -L/usr/lib/gsl
LIBS=-lwif_core -lwif_algo -lwif_viz -lgsl -lgslcblas -lm

willitfly : $(OBJ_PATH)
	echo "Building will-it-fly"
	$(CXX) -o ../build/bin/wif $(OBJ_PATH) $(LIBDIR) $(LIBS)

../build/obj/wif/%$(BUILD_SUFFIX).o : %.cpp
	$(CXX) -c $(CXX_FLAGS) -o $@ $< $(INCLUDE)
