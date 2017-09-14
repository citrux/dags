cmake_minimum_required(VERSION 2.6)
project(dags C)
add_executable(steps tests/steps.c src/wavelet.c)
target_link_libraries(steps m)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/input.dat
               ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/steps.py
${CMAKE_CURRENT_BINARY_DIR} COPYONLY)