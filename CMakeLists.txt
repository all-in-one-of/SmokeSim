cmake_minimum_required(VERSION 3.5)
# Find eigen in system directory and treat it as an interface library
# find_package(Eigen3 REQUIRED)
# ## Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH CISPBA_RELATIVE_PATH ${CISPBA_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${CISPBA_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

add_cispba_executable(SmokeSim main.cpp
                    globalIncludes.h
                    smokeSolver.h
                    smokeSolver.cpp
                    smokeQuantity.h
                    smokeQuantity.cpp
                    gridStash.h
                    gridStash.cpp
                    smokeParticles.h
                    smokeParticles.cpp
                    )

# INCLUDE_DIRECTORIES ( "$ENV{Eigen}" )
        
# target_include_directories(SmokeSim SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
# target_include_directories(SmokeSim SYSTEM PUBLIC /usr/local/include/eigen3/)
target_include_directories(SmokeSim SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Eigen)
set(CMAKE_CXX_FLAGS "-O3")
target_link_libraries(SmokeSim partio)
