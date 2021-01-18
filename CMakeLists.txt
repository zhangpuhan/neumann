cmake_minimum_required(VERSION 3.17)
project(neumann)

set(CMAKE_CXX_STANDARD 20)

find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(neumann main.cpp model.cpp model.h util.h)

target_link_libraries(neumann ${LAPACK_LIBRARIES} ${BLAS_LIBARIES} ${ARMADILLO_LIBRARIES})