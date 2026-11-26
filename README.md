## Rigid Body Kinematics Helpers
C++ library of attitude descriptions, solvers, and numeric integrators to support my self-study of Spacecraft Dynamics and Control on Coursera.
The methods below are validated by the questions in the course.

---
Contains C++ classes for the following attitude descriptions:
- Directional Cosine Matrices (Matrix3x3)
- Euler Angles (Not a class but helpers for 321EA to DCM conversions)
- Principal Rotation Vectors
- Quaternions
- Classical and Modified Rodrigues Parameters

Along with integrators to solve differential kinematic equations using these descriptions.

---
Attitude Solvers:
- TRIAD Method
- Wahba's Problem Solvers (Davenport's q-Method & QUEST Method)
- OLAE Method

Also lots of matrix/vector operation functions, including inverses and determinants.

## Building

Requirements: CMake >= 3.20 and a C++23 compiler (i.e. GCC 12+ or Clang 15+).

```bash
cmake -S . -B build
cmake --build build
./build/RBKinematicsHelpers
```

## Example

`main.cpp` estimates a spacecraft's attitude from two vector observations (i.e. sun and magnetic-field directions), measured in both the body frame and the inertial frame. It solves for the rotation matrix using OLAE and prints the resulting 3x3 DCM.
