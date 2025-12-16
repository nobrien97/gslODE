cmake_minimum_required(VERSION 3.16)
project(GSLODE)

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

add_test(
  NAME test_rk4
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_input.csv" -s rk4
)

add_test(
  NAME test_rkf45
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_input.csv" -s rkf45
)

add_test(
  NAME test_msadams
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_input.csv" -s msadams
)

add_test(
  NAME test_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_input.csv" -f 1
)

add_test(
  NAME test_rk4_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_input.csv" -s rk4 -f 0.1
)

add_test(
  NAME test_rkf45_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_input.csv" -s rkf45 -f 0.1
)

add_test(
  NAME test_msadams_f
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_input.csv" -s msadams -f 0.1
)
