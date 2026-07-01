cmake_minimum_required(VERSION 3.16)
project(jsi_ascend_hotspots LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(ASCEND_HOME $ENV{ASCEND_HOME})
if(NOT ASCEND_HOME)
  set(ASCEND_HOME "/usr/local/Ascend/ascend-toolkit/latest")
endif()

add_library(jsi_ascend_hotspots SHARED src/jsi_ascend_hotspots.cpp)
target_include_directories(jsi_ascend_hotspots PRIVATE include "${ASCEND_HOME}/include")
target_compile_options(jsi_ascend_hotspots PRIVATE -O3 -march=armv8.2-a)
target_link_directories(jsi_ascend_hotspots PRIVATE "${ASCEND_HOME}/lib64")
target_link_libraries(jsi_ascend_hotspots PRIVATE acl)
set_target_properties(jsi_ascend_hotspots PROPERTIES OUTPUT_NAME "jsi_ascend_hotspots")
