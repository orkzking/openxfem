cmake_minimum_required(VERSION 3.11-3.18)

project(openxfem)

add_executable(openxfem ${openxfem_SRC} ${openxfem_INC} ${openxfem_DEFS})

target_include_directories(openxfem PUBLIC ../inc ../defs)
