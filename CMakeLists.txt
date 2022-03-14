# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("cg_descent")

# Include sub-projects.
add_subdirectory ("cg_descent_1.1")
add_subdirectory ("cg_descent_3.0")
add_subdirectory ("cg_descent_4.0")
add_subdirectory ("cg_descent_5.0")
add_subdirectory ("cg_descent_6.0")

# TODO: 如有需要，请添加测试并安装目标。

