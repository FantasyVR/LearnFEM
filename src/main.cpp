#include "Softbody.h"
#include "World.h"
#include <string>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cat.h>
#include <chrono>
World* world;
static int globalstep = 1;
bool preDrawFunc(igl::opengl::glfw::Viewer& viewer)
{
    auto start = std::chrono::high_resolution_clock::now();
    world->step();
    auto end = std::chrono::high_resolution_clock::now();
    auto seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout<<"Time step "<< globalstep++ <<" : "<< seconds << "ms\n";
    Softbody* softbody = world->getSoftbody();
    Rigidbody* rigidbody = world->getRigidbody();
    auto sX = softbody->getMesh()->getX();
    auto rX = rigidbody->getMesh()->getX();
    auto sF = softbody->getMesh()->getFaces();
    auto rF = rigidbody->getMesh()->getFaces();
    // 将两个物体的顶点和拓扑放在一起去
    Eigen::MatrixXd  V;
    Eigen::MatrixXi  F;
    igl::cat(1, sX, rX, V);
    igl::cat(1, sF, Eigen::MatrixXi(rF.array() + sX.rows()), F);
    viewer.core().align_camera_center(V, F);
    viewer.data().set_mesh(V, F);
    return false;
}
int main() {
    world = new World();
    world->setDt(0.025);
    std::string mshFile = std::string(RESOURCE) + "tet.msh";
    Eigen::Vector3d trans(0.0,0.0,0.0);
    Eigen::Vector3d scale = Eigen::Vector3d::Ones();
    Eigen::Quaterniond orientation = Eigen::Quaterniond::Identity();
    Mesh* mesh = new Mesh(mshFile,trans,scale,orientation);
    Softbody* softbody = new Softbody(world,mesh,1000, 100000,0.4);

    // apply init change of mesh
    Eigen::MatrixXd &positions = softbody->getMesh()->getX();
    Eigen::Vector3d row3 = positions.row(3);
    positions.row(3) = Eigen::Vector3d(0,row3[1]-0.5,0);
    world->addSoftbody(softbody);

    std::string cubeFile = std::string(RESOURCE) + "cube.obj";
    trans.y() = 1.0;
    mesh = new Mesh(cubeFile, trans, scale, orientation);
    Rigidbody* rigidbody = new Rigidbody(world, mesh,1.0, orientation);
    world->addRigidbody(rigidbody);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_pre_draw = &preDrawFunc;
    viewer.core().orthographic = false;
    viewer.core().camera_zoom *= 0.4;
    viewer.launch();
    delete world;
    return 0;
}
