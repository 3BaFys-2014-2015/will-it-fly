# Makefile wif_algo

SRC_ALGO=$(wildcard *.cpp)
OBJ_ALGO=$(SRC_ALGO:.cpp=$(BUILD_SUFFIX).o)
OBJ_PATH=$(addprefix ../build/obj/wif_algo/,  $(OBJ_ALGO))

INCLUDE=-I../build/include/

#Vul GSL_LIBS aan met de GSL libraries die je gebruikt. bv. -lgsl.
GSL_LIBDIR=-L/usr/lib
GSL_LIBS=-lgsl -lgslcblas -lm
GSL_INCLUDE=-I/usr/include/gsl
ROOT_LIBS=$(shell root-config --libs)
ROOT_INCLUDE=$(shell root-config --cflags)

algo: ../build/lib/libwif_algo.so 

../build/lib/libwif_algo.so : $(OBJ_PATH)
	echo "Building wif-algo"
	$(CXX) -shared -o ../build/lib/libwif_algo.so $(OBJ_PATH) $(GSL_LIBDIR) $(GSL_LIBS) $(ROOT_INCLUDE) $(ROOT_LIBDIR) $(ROOT_LIBS)
	ln -sLf ../lib/libwif_algo.so ../build/exec/libwif_algo.so

../build/obj/wif_algo/%$(BUILD_SUFFIX).o : %.cpp
	$(CXX) -c -fPIC $(CXX_FLAGS) -o $@ $< $(INCLUDE) $(GSL_INCLUDE) $(GSL_LIBDIR) $(GSL_LIBS) $(ROOT_INCLUDE) $(ROOT_LIBDIR) $(ROOT_LIBS)

