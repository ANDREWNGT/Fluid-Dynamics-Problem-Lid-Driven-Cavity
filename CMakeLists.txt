cmake_minimum_required(VERSION 3.1)
project (FullLidSolver)


set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


add_executable(FullLidSolver LidDrivenCavitySolver.cpp LidDrivenCavity.cpp
	PoissonSolver.cpp ProgramOptions.cpp MpiConfig.cpp)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIR})
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

target_link_libraries(FullLidSolver ${BLAS_LIBRARIES} ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES}) 
