add_executable(LearnFEM main.cpp ${SRCFILES})
target_link_libraries(LearnFEM PUBLIC igl::core igl::opengl_glfw)
