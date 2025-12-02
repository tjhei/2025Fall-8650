# Set the name of the project and target:
SET(TARGET "main")

# Declare all source files the target consists of. Here, this is only:
SET(TARGET_SRC simd.cpp)

CMAKE_MINIMUM_REQUIRED(VERSION 3.2)

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

project(${TARGET})
add_executable(${TARGET} ${TARGET_SRC})

# Enable AVX2
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(${TARGET} PRIVATE -mavx2 -mfma)
endif()
