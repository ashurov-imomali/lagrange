cmake_minimum_required(VERSION 3.26)
project(cl)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(cl
        main.cpp
        main.h)
