cmake_minimum_required(VERSION 3.1)
project(LearnFEM)

set(CMAKE_CXX_STANDARD 17)

set(LearnFEM_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Dependencies)
find_package(OpenGL REQUIRED)


# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)


file(GLOB SRCFILES
        "src/*.cpp"
        "src/*.h"
        "src/SVD/*.h"
        "src/SVD/*.cpp"
        )
FOREACH(item ${SRCFILES})
    IF(${item} MATCHES "main.cpp")
        LIST(REMOVE_ITEM SRCFILES ${item})
    ENDIF(${item} MATCHES "main.cpp")
ENDFOREACH(item)
#foreach(item ${SRCFILES})
#    message(${item})
#endforeach()
#find_package(LIBIGL REQUIRED QUIET)
#find_package(EIGEN REQUIRED)
add_compile_definitions(RESOURCE="${PROJECT_SOURCE_DIR}/data/")
add_subdirectory(src/)
