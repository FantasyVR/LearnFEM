add_executable(LearnFEM main.cpp ${SRCFILES} RigidBody.h Rigidbody.cpp)
target_link_libraries(LearnFEM PUBLIC igl::core igl::opengl_glfw)
