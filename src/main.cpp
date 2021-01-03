#include "Softbody.h"
#include "World.h"
#include <string>
#include <igl/opengl/glfw/Viewer.h>
Softbody* softbody;
World* world;
bool preDrawFunc(igl::opengl::glfw::Viewer& viewer)
{
    viewer.core().align_camera_center(softbody->getMesh()->getX(), softbody->getMesh()->getFaces());
    return false;
}
int main() {
    world = new World();
    world->setDt(0.025);
    std::string mshFile = std::string(RESOURCE) + "tet.msh";
    Eigen::Vector3d trans(0.0,3.0,0.0);
    Eigen::Vector3d scale = Eigen::Vector3d::Ones();
    Eigen::Quaterniond orientation = Eigen::Quaterniond::Identity();
    Mesh* mesh = new Mesh(mshFile,trans,scale,orientation);
    softbody = new Softbody(world,mesh,1000, 100000,0.4);
    world->addSoftbody(softbody);
    world->step();
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_pre_draw = &preDrawFunc;
    viewer.core().orthographic = false;
    viewer.core().camera_zoom *= 0.4;
    viewer.data().set_mesh(softbody->getMesh()->getX(), softbody->getMesh()->getFaces());
    viewer.launch();
    return 0;
}
