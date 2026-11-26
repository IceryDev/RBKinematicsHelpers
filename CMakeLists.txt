cmake_minimum_required(VERSION 3.20)
project(RBKinematicsHelpers)

set(CMAKE_CXX_STANDARD 23)

add_executable(RBKinematicsHelpers main.cpp
        matrix3x3.h
        integrators.cpp
        integrators.h
        utils.h
        PRV.h
        vector3.h
        quaternion.h
        CRP.h
        CRP.cpp
        MRP.h
        attitudeSolvers.h
        matrix.h)
