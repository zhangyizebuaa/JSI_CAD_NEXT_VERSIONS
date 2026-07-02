CXX             ?= g++
CASE            ?= 1
THREADS         ?=
DIST_THREADS    ?= 51
INTER_THREADS   ?= 55
USE_COMPUTE_K   ?= 1
DIST_RUN_THREADS = $(if $(strip $(THREADS)),$(THREADS),$(DIST_THREADS))
INTER_RUN_THREADS = $(if $(strip $(THREADS)),$(THREADS),$(INTER_THREADS))
CXXFLAGS        ?= -O3 -std=c++17 -I include/ -I src/
LDFLAGS         ?=
USE_ASCEND_ACL  ?= 0
ASCEND_HOME     ?= /usr/local/Ascend/ascend-toolkit/latest
PERF_CXXFLAGS   ?= -O3 -g -fno-omit-frame-pointer -std=c++17 -I include/ -I src/

ifeq ($(USE_ASCEND_ACL),1)
	CXXFLAGS += -DUSE_ASCEND_ACL -I $(ASCEND_HOME)/include
	LDFLAGS += -L $(ASCEND_HOME)/lib64 -lacl
endif

all: build/test_eval build/test_intersect build/test_dist

prepare:
	@mkdir -p build

build/data_i11.o: include/nurbs_data.hpp examples/data/inter_case1-1_data.cu | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c examples/data/inter_case1-1_data.cu -o build/data_i11.o

build/data_i12.o: include/nurbs_data.hpp examples/data/inter_case1-2_data.cu | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c examples/data/inter_case1-2_data.cu -o build/data_i12.o

build/data_i21.o: include/nurbs_data.hpp examples/data/inter_case2-1_data.cu | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c examples/data/inter_case2-1_data.cu -o build/data_i21.o

build/data_i22.o: include/nurbs_data.hpp examples/data/inter_case2-2_data.cu | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c examples/data/inter_case2-2_data.cu -o build/data_i22.o

build/data_d1.o: include/nurbs_data.hpp examples/data/dist_case1.cu | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c examples/data/dist_case1.cu -o build/data_d1.o

build/data_d2.o: include/nurbs_data.hpp examples/data/dist_case2.cu | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c examples/data/dist_case2.cu -o build/data_d2.o

build/data_d3.o: include/nurbs_data.hpp examples/data/dist_case3.cu | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c examples/data/dist_case3.cu -o build/data_d3.o

build/data_d4.o: include/nurbs_data.hpp examples/data/dist_case4.cu | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c examples/data/dist_case4.cu -o build/data_d4.o

build/common.o: src/common/common.cu src/common/common.cuh | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c src/common/common.cu -o build/common.o

build/evaluation.o: src/evaluation/evaluation.cu src/evaluation/evaluation.cuh src/common/common.cuh | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c src/evaluation/evaluation.cu -o build/evaluation.o

build/distance_analysis.o: src/distance_analysis/distance_analysis.cu src/distance_analysis/distance_analysis.cuh src/evaluation/evaluation.cuh src/common/common.cuh | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c src/distance_analysis/distance_analysis.cu -o build/distance_analysis.o

build/intersection.o: src/intersection/intersection.cu src/intersection/intersection.cuh src/evaluation/evaluation.cuh src/common/common.cuh | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c src/intersection/intersection.cu -o build/intersection.o

build/ascend_hotspots.o: src/ascend/ascend_hotspots.cpp src/ascend/ascend_hotspots.hpp | prepare
	$(CXX) $(CXXFLAGS) -x c++ -c src/ascend/ascend_hotspots.cpp -o build/ascend_hotspots.o

lib: build/common.o build/evaluation.o build/distance_analysis.o build/intersection.o build/ascend_hotspots.o

data: build/data_i11.o build/data_i12.o build/data_i21.o build/data_i22.o build/data_d1.o build/data_d2.o build/data_d3.o build/data_d4.o

build/test_dist: lib data examples/distance_analysis_example.cu include/nurbs_data.hpp include/jsi_cad.hpp
	$(CXX) $(CXXFLAGS) -x c++ examples/distance_analysis_example.cu -x none build/common.o build/evaluation.o build/distance_analysis.o build/intersection.o build/ascend_hotspots.o build/data_d1.o build/data_d2.o build/data_d3.o build/data_d4.o -o build/test_dist $(LDFLAGS)

run_dist: build/test_dist
	OMP_NUM_THREADS=$(DIST_RUN_THREADS) ./build/test_dist $(CASE)

build/test_intersect: lib data examples/intersection_example.cu include/nurbs_data.hpp include/jsi_cad.hpp
	$(CXX) $(CXXFLAGS) -x c++ examples/intersection_example.cu -x none build/common.o build/evaluation.o build/distance_analysis.o build/intersection.o build/ascend_hotspots.o build/data_i11.o build/data_i12.o build/data_i21.o build/data_i22.o -o build/test_intersect $(LDFLAGS)

run_intersect: build/test_intersect
	CAD_INTERSECT_USE_K=$(USE_COMPUTE_K) OMP_NUM_THREADS=$(INTER_RUN_THREADS) ./build/test_intersect $(CASE)

run_all_cases: build/test_dist build/test_intersect
	OMP_NUM_THREADS=$(DIST_RUN_THREADS) ./build/test_dist 1
	OMP_NUM_THREADS=$(DIST_RUN_THREADS) ./build/test_dist 2
	OMP_NUM_THREADS=$(DIST_RUN_THREADS) ./build/test_dist 3
	OMP_NUM_THREADS=$(DIST_RUN_THREADS) ./build/test_dist 4
	CAD_INTERSECT_USE_K=$(USE_COMPUTE_K) OMP_NUM_THREADS=$(INTER_RUN_THREADS) ./build/test_intersect 1
	CAD_INTERSECT_USE_K=$(USE_COMPUTE_K) OMP_NUM_THREADS=$(INTER_RUN_THREADS) ./build/test_intersect 2

build/test_eval: lib data examples/evaluation_example.cu include/nurbs_data.hpp include/jsi_cad.hpp
	$(CXX) $(CXXFLAGS) -x c++ examples/evaluation_example.cu -x none build/common.o build/evaluation.o build/distance_analysis.o build/intersection.o build/ascend_hotspots.o build/data_d1.o -o build/test_eval $(LDFLAGS)

run_eval: build/test_eval
	./build/test_eval

clean:
	@if [ -d build ]; then \
		stale="build.clean.$$$$"; \
		mv build "$$stale"; \
		rm -rf "$$stale" || echo "warning: left $$stale for later cleanup"; \
	fi

profile-gprof:
	$(MAKE) clean
	$(MAKE) all CXXFLAGS="$(PERF_CXXFLAGS) -pg" LDFLAGS="$(LDFLAGS) -pg"

profile-perf:
	$(MAKE) clean
	$(MAKE) all CXXFLAGS="$(PERF_CXXFLAGS)" LDFLAGS="$(LDFLAGS)"

perf-stat-eval: profile-perf
	perf stat -d ./build/test_eval

perf-stat-dist: profile-perf
	perf stat -d ./build/test_dist

perf-stat-intersect: profile-perf
	perf stat -d ./build/test_intersect
