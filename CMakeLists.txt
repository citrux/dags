cmake_minimum_required(VERSION 2.6)
project(dags C)
add_executable(spectrum_test tests/spectrum_test.c src/wavelet.c)
target_link_libraries(spectrum_test m)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/input.dat
               ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)