#include "Softbody.h"
#include <string>
#include <igl/opengl/glfw/Viewer.h>
Softbody* softbody;
bool preDrawFunc(igl::opengl::glfw::Viewer& viewer)
{
    viewer.core().align_camera_center(softbody->getMesh()->getX(), softbody->getMesh()->getFaces());
    return false;
}
int main() {
    std::string mshFile = std::string(RESOURCE) + "tet.msh";
    Mesh* mesh = new Mesh(mshFile);
    softbody = new Softbody(mesh,1000, 100000,0.4);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_pre_draw = &preDrawFunc;
    viewer.core().orthographic = false;
    viewer.core().camera_zoom *= 0.4;
    viewer.data().set_mesh(softbody->getMesh()->getX(), softbody->getMesh()->getFaces());
    viewer.launch();
    return 0;
}
