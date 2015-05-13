# Makefile wif_core

SRC_CORE=$(wildcard *.cpp)
OBJ_CORE=$(SRC_CORE:.cpp=$(BUILD_SUFFIX).o)
OBJ_PATH=$(addprefix ../build/obj/wif_core/,  $(OBJ_CORE))

core : ../build/lib/libwif_core.so

../build/lib/libwif_core.so : $(OBJ_PATH)
	echo "Building wif-core"
	$(CXX) -shared -o ../build/lib/libwif_core.so $(OBJ_PATH)
	ln -sLf ../lib/libwif_core.so ../build/exec/libwif_core.so

../build/obj/wif_core/%$(BUILD_SUFFIX).o : %.cpp
	$(CXX) -c -fPIC $(CXX_FLAGS) -o $@  $<


	
