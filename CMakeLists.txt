cmake_minimum_required(VERSION 3.16)
project(GSLODE)

file(GLOB_RECURSE SOURCES src/*.c)

# GSL - adding /usr/local/include so all targets that use GSL_INCLUDES get omp.h
set(TARGET_NAME_GSL gsl)
file(GLOB_RECURSE GSL_SOURCES ${PROJECT_SOURCE_DIR}/gsl/*.c ${PROJECT_SOURCE_DIR}/gsl/*/*.c)
set(GSL_INCLUDES ${PROJECT_SOURCE_DIR}/gsl ${PROJECT_SOURCE_DIR}/gsl/specfunc ${PROJECT_SOURCE_DIR}/gsl/blas ${PROJECT_SOURCE_DIR}/gsl/rng ${PROJECT_SOURCE_DIR}/gsl/cdf ${PROJECT_SOURCE_DIR}/gsl/vector ${PROJECT_SOURCE_DIR}/gsl/err ${PROJECT_SOURCE_DIR}/gsl/sys ${PROJECT_SOURCE_DIR}/gsl/randist ${PROJECT_SOURCE_DIR}/gsl/matrix ${PROJECT_SOURCE_DIR}/gsl/cblas ${PROJECT_SOURCE_DIR}/gsl/complex ${PROJECT_SOURCE_DIR}/gsl/block ${PROJECT_SOURCE_DIR}/gsl/interpolation ${PROJECT_SOURCE_DIR}/gsl/linalg ${PROJECT_SOURCE_DIR}/gsl/permutation ${PROJECT_SOURCE_DIR}/gsl/ode-initval2 /usr/local/include)
add_library(${TARGET_NAME_GSL} STATIC ${GSL_SOURCES})
target_include_directories(${TARGET_NAME_GSL} PUBLIC ${GSL_INCLUDES})

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDES} PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(${PROJECT_NAME} PUBLIC gsl PRIVATE m)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)