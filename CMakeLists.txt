cmake_minimum_required(VERSION 3.13)
project(difference_cover)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(ConsoleApplication1)
include_directories(ConsoleApplication1/ConsoleApplication1)

add_executable(difference_cover
        ConsoleApplication1/ConsoleApplication1/targetver.h
        main.cpp
        main.h
        temp.cpp)
