# CMakeList.txt : CMake project for cg_descent_1.1, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (CG_DESCENT-C_4.1   "cg_descent.h" "cg_descent.c" "driver1.c")
add_executable (CG_DESCENT-C_4.2   "cg_descent.h" "cg_descent.c" "driver2.c")
add_executable (CG_DESCENT-C_4.3   "cg_descent.h" "cg_descent.c" "driver3.c")
add_executable (CG_DESCENT-C_4.4   "cg_descent.h" "cg_descent.c" "driver4.c")
add_executable (CG_DESCENT-C_4.5   "cg_descent.h" "cg_descent.c" "driver5.c")

# TODO: Add tests and install targets if needed.
