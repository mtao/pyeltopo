cmake_minimum_required(VERSION 3.11.1)
PROJECT(pyeltopo)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(PYBIND11_CPP_STANDARD -std=c++11)


include(FetchContent REQUIRED)

find_package(Eigen3 3.3.7)
if(NOT Eigen3_FOUND)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror
        GIT_TAG 3.37
        )
    FetchContent_Populate(eigen)
    add_subdirectory(${eigen_SOURCE_DIR} ${eigen_BINARY_DIR})
endif()

FetchContent_Declare(
    eltopo 
    GIT_REPOSITORY https://github.com/mtao/eltopo.git 
    GIT_TAG ad6e15c02945f48a243eb7d82cc166f810e87c86 
    )
FetchContent_Populate(eltopo)
add_subdirectory(${eltopo_SOURCE_DIR} ${eltopo_BINARY_DIR})

FetchContent_Declare(
    pybind11 
    GIT_REPOSITORY https://github.com/pybind/pybind11.git 
    GIT_TAG v2.4.3
    )
FetchContent_Populate(pybind11)
add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})


PYBIND11_ADD_MODULE(pyeltopo src/main.cpp
    src/eltopo.cpp
    include/eltopo.h
    )


TARGET_INCLUDE_DIRECTORIES(pyeltopo PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
TARGET_INCLUDE_DIRECTORIES(pyeltopo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(pyeltopo PRIVATE eltopo Eigen3::Eigen)
#target_link_libraries(pyeltopo PRIVATE pybind11::module eltopo)
#set_target_properties(pyeltopo PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
#                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")

MESSAGE(STATUS "Pyexec: ${PYTHON_EXECUTABLE}")

EXECUTE_PROCESS(
    COMMAND "${PYTHON_EXECUTABLE}" -m site --user-site
    OUTPUT_VARIABLE PYTHON_SITE
    OUTPUT_STRIP_TRAILING_WHITESPACE)
MESSAGE(STATUS "Pysites: ${PYTHON_SITE}")

MESSAGE(STATUS "TARGETS pyeltopo LIBRARY DESTINATION ${PYTHON_SITE}")
INSTALL(TARGETS pyeltopo LIBRARY DESTINATION ${PYTHON_SITE})
