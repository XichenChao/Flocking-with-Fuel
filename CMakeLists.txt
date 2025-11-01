cmake_minimum_required(VERSION 3.30)
project(cpp_online)

set(CMAKE_CXX_STANDARD 20)

include_directories(/opt/homebrew/include/eigen3)

add_executable(cpp_online main.cpp)
