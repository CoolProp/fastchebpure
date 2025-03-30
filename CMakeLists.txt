cmake_minimum_required(VERSION 3.12)
project(fastcheb)
set(CMAKE_CXX_STANDARD 20)

set(TEQP_DISABLED_FACTORIES "CPA,SAFTVRMIE,GENERICSAFT,SQUAREWELL,EXP6,2CLJF,MIE,MULTIFLUIDACTIVITY,PCSAFT,GERG200X,MULTIFLUIDASSOCIATION,LKP,ECSHUBERELY1994,AMMONIAWATERTILLNERROTH,  ADVANCEDMIXCUBIC,SOFTSAFT" CACHE STRING "Disable some unused models")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp")

set(CHEBTOOLS_STATIC_LIBRARY true CACHE INTERNAL "Build ChebTools as a static library")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/ChebTools")

add_executable(fitcheb fastcheb.cpp main.cpp)
target_link_libraries(fitcheb teqpinterface ChebTools teqpcpp)
target_include_directories(fitcheb
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp/externals/autodiff"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/externals"
)
target_compile_definitions(fitcheb PUBLIC -DTEQP_MULTIPRECISION_ENABLED)
target_compile_definitions(fitcheb PUBLIC -DTEQP_MULTICOMPLEX_ENABLED)

add_executable(bench bench.cpp)
target_link_libraries(bench teqpinterface ChebTools teqpcpp Catch2WithMain)
target_include_directories(bench
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp/externals/autodiff"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/externals"
)
target_compile_definitions(bench PUBLIC -DTEQP_MULTIPRECISION_ENABLED)
target_compile_definitions(bench PUBLIC -DTEQP_MULTICOMPLEX_ENABLED)
