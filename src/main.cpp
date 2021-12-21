#include "Cloth.h"
#include <igl/opengl/glfw/Viewer.h>
Cloth cloth(10);
double h = 0.01;
int Max_ite = 100;
bool preDrawFunc(igl::opengl::glfw::Viewer& viewer)
{
    //cloth.step(h,Max_ite);
    cloth.chebyshev_step(h, Max_ite);
//    viewer.core().align_camera_center(cloth.getX(), cloth.getFaces());
    viewer.data().clear();
    viewer.data().set_mesh(cloth.getX(), cloth.getFaces());
    return false;
}
int main() {
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_pre_draw = &preDrawFunc;
    viewer.core().orthographic = false;
    viewer.core().camera_zoom *= 0.4;
    viewer.data().set_mesh(cloth.getX(), cloth.getFaces());
    viewer.core().is_animating = true;
    viewer.launch();
    return 0;
}
