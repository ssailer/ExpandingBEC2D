project(main)
cmake_minimum_required(VERSION 2.8.4)

include_directories( ${CMAKE_SOURCE_DIR}/src
					 /usr/include/eigen3					 
					 #/home/nhrjr/libraries/libconfig-1.5/lib
					 				 
					 )

set(EXP2DSOURCES     ${CMAKE_SOURCE_DIR}/src/rk4.cpp					 					 
					 ${CMAKE_SOURCE_DIR}/src/tools.cpp
					 ${CMAKE_SOURCE_DIR}/src/evaluation.cpp
					 ${CMAKE_SOURCE_DIR}/src/plot_with_mgl.cpp	
					 ${CMAKE_SOURCE_DIR}/src/contour.cpp
					 ${CMAKE_SOURCE_DIR}/src/binaryfile.cpp
					 ${CMAKE_SOURCE_DIR}/src/splitstep.cpp					 
					 ${CMAKE_SOURCE_DIR}/src/plotter.cpp
					 ${CMAKE_SOURCE_DIR}/src/hydro.cpp
					 ${CMAKE_SOURCE_DIR}/src/lmfitter.cpp
					 ${CMAKE_SOURCE_DIR}/src/gauss_random.cpp

					 ## DEPRECATED SOURCEFILES
					 )

add_library(EXP2DLIBRARY OBJECT ${EXP2DSOURCES})

set(SOURCESMAIN main.cpp ) #${EXP2DSOURCES})

#set(SOURCESTESTER tester.cpp ) #${EXP2DSOURCES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -fopenmp -Ofast -funroll-loops -march=native -pthread -m64 -g")# -g   -ftree-vectorizer-verbose=1")  

add_executable(main ${SOURCESMAIN} $<TARGET_OBJECTS:EXP2DLIBRARY>)
#add_executable(tester ${SOURCESTESTER} $<TARGET_OBJECTS:EXP2DLIBRARY>)

target_link_libraries(main ${OTHER_LIBRARIES} -lmgl2 -lhdf5 -lconfig++ -lgsl -lgomp -lfftw3 -lboost_program_options -lpthread -lm -Wnarrowing)
	#-L/usr/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread  -lgslcblas -ldl -pthread -rdynamic  -lfftw3_thread
	   

#target_link_libraries(tester ${OTHER_LIBRARIES} -lmgl2 -lhdf5 -lconfig++ -lgsl -lgslcblas  -lgomp -lpthread -lfftw3 -lfftw3_threads -lm -lboost_program_options 
#	-L/usr/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -Wnarrowing)

install(TARGETS main RUNTIME DESTINATION bin)
#install(TARGETS tester RUNTIME DESTINATION bin)
