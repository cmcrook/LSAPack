# only for cmake --version >= 3.5.1
cmake_minimum_required(VERSION 3.5.1)

# project name
project(Packing)

# creates the variable EXEC and sets it to hello
set(EXEC lsa)

# set the C++14 standard
set(CMAKE_CXX_STANDARD 17)

# I../includes
include_directories(src)

# puts all .cpp files inside src to the SOURCES variable
file(GLOB SOURCES src/*.cpp)

# compiles the files defined by SOURCES to generante the executable defined by EXEC
add_executable(${EXEC} ${SOURCES})

IF (WIN32)
message("This is windows")
    add_compile_definitions(strcasecmp=_stricmp)
ENDIF()