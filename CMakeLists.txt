cmake_minimum_required(VERSION 4.2)
project(RBKinematicsHelpers)

set(CMAKE_CXX_STANDARD 23)

add_executable(RBKinematicsHelpers main.cpp
        matrix3x3.h
        integrators.cpp
        integrators.h
        utils.h
        PRV.h
        vector3.h)
