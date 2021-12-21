#file(GLOB SRCFILES "./*.cpp" "./*.h")
add_executable(LearnFEM main.cpp Cloth.h Cloth.cpp)
target_link_libraries(LearnFEM PUBLIC igl::core igl::opengl_glfw)