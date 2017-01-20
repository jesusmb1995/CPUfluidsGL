#
# Make for the CPU version of 
#  NVidia's fluids demo 
#
#           By Jesus Martin Berlanga
#
#####################################
#####################################
#
# Minimal requirements: 
#    - FFTW library: http://www.fftw.org/
#          (Install at Ubuntu: `sudo apt-get install libfftw3-dev`)
#    - OpenGL library
#        + GLEW library
#	   (Install at Ubuntu: `sudo apt-get install libglew-dev`)
#    - g++ compiler
#
# Suggested (optional):
#     - gprof
#     - Intel C++ compiler with a 64 bits Intel CPU
#          (Student versions https://software.intel.com/en-us/qualify-for-free-software/student
#	    You can use Parallel Studio XE and only install the Intel compiler.)
#
# Tested with:
#     - OS: Ubuntu 14.04.4 LTS 	64 bits
#         + kernel: Linux 3.13.0-57-generic (x86_64)
#     - CPU: 8x Intel(R) Core(TM) i7 CPU 	960 @ 3.20GHz
#     - RAM: 18488 MB
#     - GPU: GeForce GTX 780
#         + CUDA Capability 3.5
#         + NVIDIA Driver Version: 361.93.02
#     - Compilers:
#         + g++ version 4.8.4 (Ubuntu 4.8.4-2ubuntu1~14.04.3)
#         + icpc 17.0.1
#
#####################################

INTEL_COMPILERS_FOLDER = /opt/intel/compilers_and_libraries/linux/bin/intel64
INTEL_ICPC = "$(INTEL_COMPILERS_FOLDER)/icpc"

BINARY_DIR = bin
FPS_RESULTS_DIR = fps_tests
PROFILE_RESULTS_DIR = profiling

CXXFLAGS = -Wall

all: build
	make clean-objects; make build-debug
	make clean-objects; make build-profile 
	make clean-objects

build: 
	make clean-objects; make build-gcc-O0
	make clean-objects; make build-gcc-O1
	make clean-objects; make build-gcc-O2
	make clean-objects; make build-gcc-O3
	make clean-objects; make build-icc-O0
	make clean-objects; make build-icc-O1
	make clean-objects; make build-icc-O2
	make clean-objects; make build-icc-O3
	make clean-objects; make build-icc-O3plus
	make clean-objects

build-profile:
	make clean-objects; make build-profile-gcc-O0
	make clean-objects; make build-profile-gcc-O1
	make clean-objects; make build-profile-icc-O3plus
	make clean-objects

build-gcc-O0: CXX = g++
build-gcc-O0: BIN_POST_NAME = -gcc-O0
build-gcc-O0: OPT_FLAGS = -O0
build-gcc-O0: fluidsGL

build-debug: DBG_NAME_APPEND = -dbg
build-debug: CXXFLAGS += -g
build-debug: build-gcc-O0

build-profile-gcc-O0: DBG_NAME_APPEND = -dbg-gprof
build-profile-gcc-O0: CXXFLAGS += -g -pg
build-profile-gcc-O0: build-gcc-O0

build-profile-gcc-O1: DBG_NAME_APPEND = -dbg-gprof
build-profile-gcc-O1: CXXFLAGS += -g -pg
build-profile-gcc-O1: build-gcc-O1

build-profile-icc-O3plus: DBG_NAME_APPEND = -dbg-gprof
build-profile-icc-O3plus: CXXFLAGS += -g -pg
build-profile-icc-O3plus: build-icc-O3plus-final

build-gcc-O1: CXX = g++
build-gcc-O1: BIN_POST_NAME = -gcc-O1
build-gcc-O1: OPT_FLAGS = -O1
build-gcc-O1: fluidsGL

build-gcc-O2: CXX = g++
build-gcc-O2: BIN_POST_NAME = -gcc-O2
build-gcc-O2: OPT_FLAGS = -O2
build-gcc-O2: fluidsGL

build-gcc-O3: CXX = g++
build-gcc-O3: BIN_POST_NAME = -gcc-O3
build-gcc-O3: OPT_FLAGS = -O3
build-gcc-O3: fluidsGL

build-icc-O0: CXX = $(INTEL_ICPC)
build-icc-O0: BIN_POST_NAME = -icc-O0
build-icc-O0: OPT_FLAGS = -O0
build-icc-O0: fluidsGL

build-icc-O1: CXX = $(INTEL_ICPC)
build-icc-O1: BIN_POST_NAME = -icc-O1
build-icc-O1: OPT_FLAGS = -O1
build-icc-O1: fluidsGL

build-icc-O2: CXX = $(INTEL_ICPC)
build-icc-O2: BIN_POST_NAME = -icc-O2
build-icc-O2: OPT_FLAGS = -O2
build-icc-O2: fluidsGL

build-icc-O3: CXX = $(INTEL_ICPC)
build-icc-O3: BIN_POST_NAME = -icc-O3
build-icc-O3: OPT_FLAGS = -O3
build-icc-O3: fluidsGL

build-icc-O3plus: CXX = $(INTEL_ICPC)
build-icc-O3plus: BIN_POST_NAME = -icc-O3plus
# First build pass (with modest flags. e.g O0) to generate profiling files 
build-icc-O3plus: OPT_FLAGS = -prof-gen -O0 -qopt-malloc-options=3 -ansi-alias -qopt-prefetch -prof-dir=$(PROFILE_RESULTS_DIR) 
build-icc-O3plus: fluidsGL
# Run test to generate profiling files
build-icc-O3plus: fluidsGL-icc-O3plus.stress_test.txt  
build-icc-O3plus:
	rm $(BINARY_DIR)/fluidsGL-icc-O3plus
	make clean-objects; make build-icc-O3plus-final
# And then compile to optimize with the profiling files
#  use all possible optimization flags for the final build
build-icc-O3plus-final: CXX = $(INTEL_ICPC)
build-icc-O3plus-final: BIN_POST_NAME = -icc-O3plus
####################################################################
# SSE optimization Compiler options
# for Nehalem use -xSSE4.2
# for processors prior to dunnington, replace -xSSE4.1 with -xSSSE3
####################################################################
SSE = -xSSSE3
build-icc-O3plus-final: OPT_FLAGS = $(SSE) -O3 -no-prec-div -prof-use -qopt-malloc-options=3 -ansi-alias -qopt-prefetch -prof-dir=$(PROFILE_RESULTS_DIR)  
build-icc-O3plus-final: fluidsGL

defines.o: defines.c
	$(CXX) $(CXXFLAGS) -c defines.c -o $(BINARY_DIR)/$@
bilinear_interpolation.o: bilinear_interpolation.cpp
	$(CXX) $(CXXFLAGS) $(OPT_FLAGS) -c bilinear_interpolation.cpp -o $(BINARY_DIR)/$@
fluidsGL.o: fluidsGL.cpp
	$(CXX) $(CXXFLAGS) -c fluidsGL.cpp -o $(BINARY_DIR)/$@
fluidsGL_cpu.o: fluidsGL_cpu.cpp
	$(CXX) $(CXXFLAGS) $(OPT_FLAGS) -c fluidsGL_cpu.cpp -o $(BINARY_DIR)/$@
fluidsGL: fluidsGL.o fluidsGL_cpu.o bilinear_interpolation.o defines.o
	$(CXX) $(CXXFLAGS) $(BINARY_DIR)/defines.o $(BINARY_DIR)/bilinear_interpolation.o $(BINARY_DIR)/fluidsGL_cpu.o $(BINARY_DIR)/fluidsGL.o -o $(BINARY_DIR)/$@$(BIN_POST_NAME)$(DBG_NAME_APPEND) -lGL -lGLU -lGLEW -lglut -lfftw3f


clean-objects:
	rm -f $(BINARY_DIR)/*.o
clean: clean-objects
	rm -f $(BINARY_DIR)/*gcc* $(BINARY_DIR)/*icc*

stress_tests: fluidsGL-gcc-O0.stress_test.txt fluidsGL-gcc-O1.stress_test.txt fluidsGL-gcc-O2.stress_test.txt fluidsGL-gcc-O3.stress_test.txt fluidsGL-icc-O0.stress_test.txt fluidsGL-icc-O1.stress_test.txt fluidsGL-icc-O2.stress_test.txt fluidsGL-icc-O3.stress_test.txt fluidsGL-icc-O3plus.stress_test.txt
fluidsGL-gcc-O0.stress_test.txt:
	$(BINARY_DIR)/fluidsGL-gcc-O0 stress_test > $(FPS_RESULTS_DIR)/$@
fluidsGL-gcc-O1.stress_test.txt:
	$(BINARY_DIR)/fluidsGL-gcc-O1 stress_test > $(FPS_RESULTS_DIR)/$@
fluidsGL-gcc-O2.stress_test.txt:
	$(BINARY_DIR)/fluidsGL-gcc-O2 stress_test > $(FPS_RESULTS_DIR)/$@
fluidsGL-gcc-O3.stress_test.txt:
	$(BINARY_DIR)/fluidsGL-gcc-O3 stress_test > $(FPS_RESULTS_DIR)/$@
fluidsGL-icc-O0.stress_test.txt:
	$(BINARY_DIR)/fluidsGL-icc-O0 stress_test > $(FPS_RESULTS_DIR)/$@
fluidsGL-icc-O1.stress_test.txt:
	$(BINARY_DIR)/fluidsGL-icc-O1 stress_test > $(FPS_RESULTS_DIR)/$@
fluidsGL-icc-O2.stress_test.txt:
	$(BINARY_DIR)/fluidsGL-icc-O2 stress_test > $(FPS_RESULTS_DIR)/$@
fluidsGL-icc-O3.stress_test.txt:
	$(BINARY_DIR)/fluidsGL-icc-O3 stress_test > $(FPS_RESULTS_DIR)/$@
fluidsGL-icc-O3plus.stress_test.txt:
	$(BINARY_DIR)/fluidsGL-icc-O3plus stress_test > $(FPS_RESULTS_DIR)/$@

profile:
	make profile-gcc-O0;
	make profile-gcc-O1;
	make profile-icc-O3plus;

profile-gcc-O0:
	$(BINARY_DIR)/fluidsGL-gcc-O0-dbg-gprof stress_test
	mv gmon.out $(PROFILE_RESULTS_DIR)/gmon-gcc-O0-dbg.out
	gprof $(BINARY_DIR)/fluidsGL-gcc-O0-dbg-gprof $(PROFILE_RESULTS_DIR)/gmon-gcc-O0-dbg.out > $(PROFILE_RESULTS_DIR)/fluidsGL-gcc-O0-dbg-gprof.gprof.txt

profile-gcc-O1:
	$(BINARY_DIR)/fluidsGL-gcc-O1-dbg-gprof stress_test
	mv gmon.out $(PROFILE_RESULTS_DIR)/gmon-gcc-O1-dbg.out
	gprof $(BINARY_DIR)/fluidsGL-gcc-O1-dbg-gprof $(PROFILE_RESULTS_DIR)/gmon-gcc-O1-dbg.out > $(PROFILE_RESULTS_DIR)/fluidsGL-gcc-O1-dbg-gprof.gprof.txt

profile-icc-O3Plus:
	$(BINARY_DIR)/fluidsGL-icc-O3plus-dbg-gprof stress_test
	mv gmon.out $(PROFILE_RESULTS_DIR)/gmon-gcc-icc-O3Plus-dbg.out
	gprof $(BINARY_DIR)/fluidsGL-icc-O3plus-dbg-gprof $(PROFILE_RESULTS_DIR)/gmon-gcc-icc-O3Plus-dbg.out > $(PROFILE_RESULTS_DIR)/fluidsGL-icc-O3Plus-dbg-gprof.gprof.txt

reall: clean
reall: all

rebuild: clean
rebuild: build

redebug: clean
redebug: build-debug

reprofile: clean 
reprofile: build-profile

