cmake_minimum_required(VERSION 3.16)
project(GSLODE)

set(CMAKE_C_FLAGS_RELEASE "-O3 -Wall")

list(FILTER SOURCES EXCLUDE REGEX "gsl/*")
file(GLOB_RECURSE SOURCES src/*.c)

# # GSL - adding /usr/local/include so all targets that use GSL_INCLUDES get omp.h
# set(TARGET_NAME_GSL gsl)
# file(GLOB_RECURSE GSL_SOURCES ${PROJECT_SOURCE_DIR}/gsl/*.c ${PROJECT_SOURCE_DIR}/gsl/*/*.c)
# set(GSL_INCLUDES ${PROJECT_SOURCE_DIR}/gsl ${PROJECT_SOURCE_DIR}/gsl/specfunc ${PROJECT_SOURCE_DIR}/gsl/blas ${PROJECT_SOURCE_DIR}/gsl/rng ${PROJECT_SOURCE_DIR}/gsl/cdf ${PROJECT_SOURCE_DIR}/gsl/vector ${PROJECT_SOURCE_DIR}/gsl/err ${PROJECT_SOURCE_DIR}/gsl/sys ${PROJECT_SOURCE_DIR}/gsl/randist ${PROJECT_SOURCE_DIR}/gsl/matrix ${PROJECT_SOURCE_DIR}/gsl/cblas ${PROJECT_SOURCE_DIR}/gsl/complex ${PROJECT_SOURCE_DIR}/gsl/block ${PROJECT_SOURCE_DIR}/gsl/interpolation ${PROJECT_SOURCE_DIR}/gsl/linalg ${PROJECT_SOURCE_DIR}/gsl/permutation ${PROJECT_SOURCE_DIR}/gsl/ode-initval2 /usr/local/include)
# add_library(${TARGET_NAME_GSL} STATIC ${GSL_SOURCES})
# target_include_directories(${TARGET_NAME_GSL} PUBLIC ${GSL_INCLUDES})

find_package(GSL REQUIRED)

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(${PROJECT_NAME} PRIVATE m GSL::gsl)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

# Test on some data
# add testing so it can be called using make test
include(CTest)
enable_testing()
message(STATUS "${PROJECT_SOURCE_DIR}")

# NAR tests
add_test(
  NAME test_nar_rk4
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_nar_input.csv" -s rk4
)

add_test(
  NAME test_nar_msbdf
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_nar_input.csv" -s bsimp -f 0.1
)

add_test(
  NAME test_nar_rkf45
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_nar_input.csv" -s rkf45
)

add_test(
  NAME test_nar_msadams
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_nar_input.csv" -s msadams
)

add_test(
  NAME test_nar_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_nar_input.csv" -f 1
)

add_test(
  NAME test_nar_rk4_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_nar_input.csv" -s rk4 -f 0.1
)

add_test(
  NAME test_nar_rkf45_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_nar_input.csv" -s rkf45 -f 0.1
)

add_test(
  NAME test_nar_msadams_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_nar_input.csv" -s msadams -f 0.1
)


# Van der Pol
add_test(
  NAME test_vdp_rk4_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_vdp_input.csv" -s rk4 -f 0.1 -o VanDerPol
)

add_test(
  NAME test_vdp_rkf45_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_vdp_input.csv" -s rkf45 -f 0.1 -o VanDerPol
)

add_test(
  NAME test_vdp_msadams_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_vdp_input.csv" -s msadams -f 0.1 -o VanDerPol
)

add_test(
  NAME test_vdp_msbdf_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_vdp_input.csv" -s bsimp -f 0.1 -o VanDerPol
)


# Lorenz
add_test(
  NAME test_lorenz_rk4_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_lorenz_input.csv" -s rk4 -f 0.1 -o Lorenz
)

add_test(
  NAME test_lorenz_rkf45_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_lorenz_input.csv" -s rkf45 -f 0.1 -o Lorenz
)

add_test(
  NAME test_lorenz_msadams_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_lorenz_input.csv" -s msadams -f 0.1 -o Lorenz
)

add_test(
  NAME test_lorenz_msbdf_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_lorenz_input.csv" -s msbdf -f 0.1 -o Lorenz
)

# Robertson: this is a stiff system, should not work except with implicit solvers
add_test(
  NAME test_rob_rk4_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_rob_input.csv" -s rk4 -f 0.1 -o Robertson
)

add_test(
  NAME test_rob_rkf45_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_rob_input.csv" -s rkf45 -f 0.1 -o Robertson
)

add_test(
  NAME test_rob_msbdf_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_rob_input.csv" -s msbdf -f 0.1 -o Robertson
)

add_test(
  NAME test_errors
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_rob_input.csv" -s msbdf -f 0.1 -o Robertson -a 1e-8 -r 1e-4
)

add_test(
  NAME test_benchmark
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_rob_input.csv" -s msbdf -f 0.1 -o Robertson -b
)