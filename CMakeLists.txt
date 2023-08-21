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
add_library(wavelib SHARED 
    src/conv.c
    src/conv.h
    src/hsfft.c
    src/hsfft.h
    src/real.c
    src/real.h
    src/wavefilt.c
    src/wavefilt.h
    src/wavelib.c
    src/wtmath.c
    src/wtmath.h
)
target_include_directories(wavelib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(wavelib PRIVATE -O3 -fPIC)
