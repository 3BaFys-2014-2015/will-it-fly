# Makefile

SRC_DEMO=$(wildcard *.cpp)
OBJ_DEMO=$(SRC_DEMO:.cpp=$(BUILD_SUFFIX).o)
OBJ_PATH=$(addprefix ../../build/obj/demos/$(DEMO_NAME)/, $(OBJ_DEMO))

INCLUDE=-I../../build/include -I/usr/include/gsl
LIBDIR=-L../../build/lib -L/usr/lib/ -L/usr/lib/gsl
LIBS=-lwif_core -lwif_algo -lwif_viz -lgsl -lgslcblas -lm

demo : ../../build/bin/wif_$(DEMO_NAME)

../../build/bin/wif_$(DEMO_NAME) : $(OBJ_PATH)
	echo "Building ${DEMO_NAME}"
	$(CXX) $(CXX_FLAGS) -o ../../build/bin/wif_$(DEMO_NAME) $(OBJ_PATH) $(INCLUDE) $(LIBDIR) $(LIBS)

../../build/obj/demos/$(DEMO_NAME)/%$(BUILD_SUFFIX).o : %.cpp
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE) -o $@ $<
