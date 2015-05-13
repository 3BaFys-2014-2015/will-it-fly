# Makefile wif_viz

SRC_VIZ=$(wildcard *.cpp)
OBJ_VIZ=$(SRC_VIZ:.cpp=$(BUILD_SUFFIX).o)
OBJ_PATH=$(addprefix ../build/obj/wif_viz/, $(OBJ_VIZ))

INCLUDE=-I../build/include/
LIBDIR=-L../build/lib

#Vul VTK_LIBS aan met de vtk libraries die je gebruikt. bv. -lvtkCommon.
VTK_LIBDIR=-L/usr/lib/
VTK_LIBS=-lwif_core -lwif_algo -lvtkRendering -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtkIO -lvtkHybrid -lvtkWidgets
VTK_INCLUDE=-I/usr/include/vtk-5.8/
ROOT_LIBS=$(shell root-config --libs)
ROOT_INCLUDE=$(shell root-config --cflags)

viz: ../build/lib/libwif_viz.so

../build/lib/libwif_viz.so : $(OBJ_PATH)
	echo "Building wif-viz"
	$(CXX)  -shared -o ../build/lib/libwif_viz.so $(OBJ_PATH) $(LIBDIR) $(VTK_LIBDIR) $(VTK_LIBS) $(ROOT_LIBDIR) $(ROOT_LIBS)
	ln -sLf ../lib/libwif_viz.so ../build/exec/libwif_viz.so

../build/obj/wif_viz/%$(BUILD_SUFFIX).o : %.cpp
	$(CXX) -c -fPIC $(CXX_FLAGS) -o $@ $< $(INCLUDE) $(VTK_INCLUDE) $(ROOT_INCLUDE) $(VTK_LIBDIR) $(ROOT_LIBS) $(VTK_LIBS)

