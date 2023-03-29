cmake_minimum_required(VERSION 3.11-3.18)

file(GLOB openxfem_SRC CMAKE_CONFIGURE_DEPENDS "*.cpp")

project(openxfem)

add_executable(openxfem ${openxfem_SRC})

target_include_directories(openxfem PUBLIC ../inc ../defs)
