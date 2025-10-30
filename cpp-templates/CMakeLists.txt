# Set the name of the project and target:
SET(TARGET "templates")

CMAKE_MINIMUM_REQUIRED(VERSION 3.2)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

project(${TARGET})
add_executable(vector-norm vector-norm.cc)
add_executable(square square.cc)
add_executable(functiontemplate functiontemplate.cc)
add_executable(explicit explicit.cc)
