cmake_minimum_required(VERSION 3.11-3.18)

project(
    openxfem
    VERSION 1.0
    DESCRIPTION "eXtended Finite Element Tool"
)

file(GLOB openxfem_SRC CMAKE_CONFIGURE_DEPENDS "src/*.cpp")
file(GLOB openxfem_INC CMAKE_CONFIGURE_DEPENDS "inc/*.h")
file(GLOB openxfem_DEFS CMAKE_CONFIGURE_DEPENDS "defs/*.DEF")

add_subdirectory(src)