cmake_minimum_required(VERSION 2.8.9)
PROJECT(pyeltopo)

SET(PYBIND11_CPP_STANDARD -std=c++11)

ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/extern/eltopo)


#ADD_LIBRARY(pyeltopo MODULE pyeltopo/main.cpp
#    pyeltopo/eltopo.cpp)
PYBIND11_ADD_MODULE(pyeltopo src/main.cpp
    src/eltopo.cpp
    include/eltopo.h
    )

SET(ELTOPO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/eltopo")

TARGET_INCLUDE_DIRECTORIES(pyeltopo PRIVATE .)
TARGET_INCLUDE_DIRECTORIES(pyeltopo PRIVATE include)
TARGET_INCLUDE_DIRECTORIES(pyeltopo PRIVATE /usr/include/eigen3)
TARGET_INCLUDE_DIRECTORIES(pyeltopo PRIVATE
    ${ELTOPO_DIR})

target_link_libraries(pyeltopo PRIVATE eltopo)
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
