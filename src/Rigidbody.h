//
// Created by yupeng on 2021/2/1.
//

#ifndef LEARNFEM_RIGIDBODY_H
#define LEARNFEM_RIGIDBODY_H
#include "Mesh.h"
class World;
class Rigidbody{
private:
    World* world;
    Mesh* mesh;
    Eigen::Vector3d cos, cos_pre; // center of mass
    Eigen::Quaterniond orientation, orientation_pre;
    double mass;
    Eigen::Vector3d vel_linear, vel_angular;
    Eigen::Matrix3d inertia_local;

    // simulation variables
    int numSubSteps;
    int numIteration;
public:
    Rigidbody(){};
    Rigidbody(World* world, Mesh* mesh, double mass,const Eigen::Quaterniond &rotation):world(world),mesh(mesh),mass(mass), orientation(rotation)
    {
        computeCenterOfMass(mesh, cos);
        cos_pre = cos;
        auto &x = mesh->getX();
        auto box = x.colwise().maxCoeff() - x.colwise().minCoeff();
        // TODO: local_inertia of triangle mesh rigid body
        computeInertiaTensorBox(mass,box[0],box[1],box[2]);

        // init general velocity: linear_vel and angular_vel
        vel_linear.setZero();
        vel_angular.setZero();

        numSubSteps = 5;
        numIteration = 1;
    }
    // Position based rigid body simulation
    void update();

    const Mesh* getMesh() const{return mesh;}
    Mesh* getMesh() {return mesh;}

private:
    void computeCenterOfMass(const Mesh* mesh,Eigen::Vector3d &cos);
    void computeInertiaTensorBox(const double mass, const double width, const double height, const double depth);
};

#endif //LEARNFEM_RIGIDBODY_H
