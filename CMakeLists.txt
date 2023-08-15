cmake_minimum_required(VERSION 3.16)
project(atar)

#--------------------------------------
# C++ standard
#--------------------------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#--------------------------------------
# Compiler Flags
#--------------------------------------
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CUSTOM_COMPILER_FLAGS}")


#--------------------------------------
# Includes
#--------------------------------------
include_directories(
    src
    include
)

#--------------------------------------
# Wavelib Library
#--------------------------------------
add_library(wavelib OBJECT 
    src/conv.c
    src/hsfft.c
    src/real.c
    src/wavefilt.c
    src/wavelib.c
    src/wtmath.c
)
target_compile_options(wavelib PRIVATE -O3 -fPIC)
