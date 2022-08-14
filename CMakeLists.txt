cmake_minimum_required(VERSION 3.23.3)
project(_nrm)

add_subdirectory(pybind11)
pybind11_add_module(_nrm src/nrm.cpp)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(_nrm PRIVATE Eigen3::Eigen)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(_nrm
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
