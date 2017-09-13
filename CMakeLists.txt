cmake_minimum_required(VERSION 2.6)
project(dags C)
add_executable(steps tests/steps.c src/wavelet.c)
target_link_libraries(steps m)