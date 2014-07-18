project(main)
cmake_minimum_required(VERSION 2.8)

include_directories( ${CMAKE_SOURCE_DIR}/exp_utilities
					 ${CMAKE_SOURCE_DIR}/utilities
					 ${CMAKE_SOURCE_DIR}/startgrids
					 ${CMAKE_SOURCE_DIR}/propagators
					 ${CMAKE_SOURCE_DIR}/observables)

set(SOURCESMAIN main.cpp ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_itp.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_rte.cpp					 
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_tools.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_evaluation.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/plot_with_mgl.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/bh3binaryfile.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/complexgrid.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/realgrid.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/gauss_random.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/realgrid.cpp
					 ${CMAKE_SOURCE_DIR}/startgrids/bh3defaultgrid.cpp
					 ${CMAKE_SOURCE_DIR}/observables/bh3observables.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_Contour.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_binaryfile.cpp
					 )

set(SOURCESPLOTTER evaluator.cpp ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_itp.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_rte.cpp					 
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_tools.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_evaluation.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/plot_with_mgl.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/bh3binaryfile.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/complexgrid.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/realgrid.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/gauss_random.cpp
					 ${CMAKE_SOURCE_DIR}/utilities/realgrid.cpp
					 ${CMAKE_SOURCE_DIR}/startgrids/bh3defaultgrid.cpp
					 ${CMAKE_SOURCE_DIR}/observables/bh3observables.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_Contour.cpp
					 ${CMAKE_SOURCE_DIR}/exp_utilities/EXP2D_binaryfile.cpp
					 )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fopenmp -Ofast -funroll-loops -march=native ") #-ftree-vectorizer-verbose=1")

add_executable(main ${SOURCESMAIN})
add_executable(evaluator ${SOURCESPLOTTER})

target_link_libraries(main ${OTHER_LIBRARIES} -lmgl -lhdf5 -lconfig++ -lgsl -lgslcblas -lgomp -lpthread -lfftw3 -lfftw3_threads -lm -lboost_program_options -Wnarrowing)
target_link_libraries(evaluator ${OTHER_LIBRARIES} -lmgl -lhdf5 -lconfig++ -lgsl -lgslcblas -lgomp -lpthread -lfftw3 -lfftw3_threads -lm -lboost_program_options -Wnarrowing)

install(TARGETS main RUNTIME DESTINATION bin)
install(TARGETS evaluator RUNTIME DESTINATION bin)
