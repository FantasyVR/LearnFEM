#file(GLOB SRCFILES "./*.cpp" "./*.h")
add_executable(LearnFEM main.cpp Mesh.h Mesh.cpp Softbody.cpp Softbody.h)
target_link_libraries(LearnFEM PUBLIC igl::core igl::opengl_glfw)
