//
// Created by yupeng on 2021/2/1.
//
#include "Rigidbody.h"
#include "World.h"
void Rigidbody::update()
{
    //TODO: Collision detection

    double h = world->getDt() / numSubSteps;
    for(int step = 0; step< numSubSteps; step++)
    {
        //Semi-Euler update
        cos_pre  = cos;
        vel_linear += h * world->getGravity();
        cos += h * vel_linear;
        orientation_pre = orientation;
        auto inertia_world = orientation.toRotationMatrix() * inertia_local * orientation.toRotationMatrix().transpose(); // TODO: validation
        vel_angular += h * inertia_world.inverse() * (- vel_angular.cross(vel_angular)); // TODO: validation
        Eigen::Quaterniond omega(0.0, vel_angular[0],vel_angular[1],vel_angular[2]);
        orientation.coeffs() +=  0.5 * h *  (omega * orientation).coeffs();
        orientation.normalize();
        //TODO: solve constraints
        for(int ite=0; ite < numIteration; ite++)
        {
            if( cos.y() < 0.0)
            {
                cos.y() = 0.0;
            }
        }
        //update velocities
        vel_linear = (cos - cos_pre)/h;
        auto Delta_q = orientation * orientation_pre.inverse(); // Do we need normalization?
        Eigen::Vector3d dq(Delta_q.x(), Delta_q.y(),Delta_q.z());
        vel_angular = 2 * dq / h;
        vel_angular = Delta_q.w() > 0.0 ? vel_angular: -vel_angular;
        //TODO: Post processing of velocities
    }
}


void Rigidbody::computeCenterOfMass(const Mesh* mesh,Eigen::Vector3d &cos)
{
    const Eigen::MatrixXd &x = mesh->getX();
    int numVert = x.rows();
    cos.setZero();
    auto temp = x.colwise().sum();
    cos = temp.transpose()/numVert;
}

void Rigidbody::computeInertiaTensorBox(const double mass, const double width, const double height, const double depth)
{
    const double Ix = (mass / static_cast<double>(12.0)) * (height * height + depth * depth);
    const double Iy = (mass / static_cast<double>(12.0)) * (width * width + depth * depth);
    const double Iz = (mass / static_cast<double>(12.0)) * (width * width + height * height);
    inertia_local = Eigen::DiagonalMatrix<double,3,3>(Ix,Iy,Iz);
}
