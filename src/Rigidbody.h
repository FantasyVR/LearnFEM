//
// Created by yupeng on 2021/2/1.
//

#ifndef LEARNFEM_RIGIDBODY_H
#define LEARNFEM_RIGIDBODY_H
#include "Mesh.h"
class World;
class Rigidbody{
private:
    World* world{};
    Mesh* mesh{};
    Eigen::Vector3d cos; // center of mass
    Eigen::Quaterniond orientation;
    double mass{};
    Eigen::Vector3d linear_vel, angular_vel;
    Eigen::Matrix3d local_inertia;

    // simulation variables
    int numSubSteps;
    int numIteration;
public:
    Rigidbody(){};
    Rigidbody(World* world, Mesh* mesh, double mass,const Eigen::Quaterniond &rotation):world(world),mesh(mesh),mass(mass), orientation(rotation)
    {
        //TODO: init general position: cos and orientation

        //TODO: init local_inertia

        //TODO: init general velocity: linear_vel and angular_vel


        numSubSteps = 5;
        numIteration = 1;
    }
    // Position based rigid body simulation
    void update();

    const Mesh* getMesh() const{return mesh;}
    Mesh* getMesh() {return mesh;}

};

#endif //LEARNFEM_RIGIDBODY_H
