#include "Softbody.h"
#include <iostream>
Softbody::Softbody(Mesh *mesh, double density, double YM, double PR)
{
    this->mesh = mesh;
    this->youngModule = YM;
    this->possionRatio = PR;
    this->density = density;
    this->mu = YM / 2.0 / (1.0 + PR);
    this->lambda = YM * PR / (1.0 + PR) / (1.0 - 2.0 * PR);
    int numTets = mesh->getNumTets();
    vels.setZero(numTets * 3);
    gradient.setZero(numTets * 3);
    hessian.setZero(numTets*3, numTets*3);
    computeB();
    computeMassMatrix();
}

void Softbody::computeB() {
    int numTets = mesh->getNumTets();
    const Eigen::MatrixXd &positions = mesh->getX();
    const Eigen::MatrixXi &tets = mesh->getTets();
    for(int i = 0; i<numTets; i++)
    {
        auto tetIdx = tets.row(i);
        const Eigen::Vector3d& p0 = positions.row(tetIdx[0]);
        const Eigen::Vector3d& p1 = positions.row(tetIdx[1]);
        const Eigen::Vector3d& p2 = positions.row(tetIdx[2]);
        const Eigen::Vector3d& p3 = positions.row(tetIdx[3]);
        Eigen::Matrix3d D;
        D.col(0) = p1 - p0;
        D.col(1) = p2 - p0;
        D.col(2) = p3 - p0;
        B.emplace_back(D.inverse());
    }
}

void Softbody::computeMassMatrix() {
    const Eigen::MatrixXd &positions = mesh->getX();
    const Eigen::VectorXd &rest_volumes = mesh->getRestVolume();
    const Eigen::MatrixXi &tets = mesh->getTets();
    mass.setZero(positions.rows());
    for(int i = 0; i<tets.rows(); i++)
    {
         auto tetIdx = tets.row(i);
         auto rest_volume = rest_volumes(i);
         mass(tetIdx[0]) += rest_volume/4.0;
         mass(tetIdx[1]) += rest_volume/4.0;
         mass(tetIdx[2]) += rest_volume/4.0;
         mass(tetIdx[3]) += rest_volume/4.0;
    }
    mass *= density;
    std::cout<<"Mass Matrix: \n"<< mass <<std::endl;
}

Softbody::~Softbody() {
    delete mesh;
}
