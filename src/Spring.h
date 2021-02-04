//
// Created by yupeng on 2021/2/3.
//

#ifndef LEARNFEM_SPRING_H
#define LEARNFEM_SPRING_H
#include "Mesh.h"
#include <vector>
class World;
struct Mass
{
    Mass(double mass, Eigen::Vector3d pos, bool isFix):mass(mass),position(pos),isFixed(isFix)
    {
        position_pre = pos;
        velocity = Eigen::Vector3d::Zero();
    }
    void resetForce(){force = Eigen::Vector3d::Zero();}
    void resetP(){P = Eigen::Matrix3d::Zero();}
    double mass;
    Eigen::Vector3d  position,position_pre;
    Eigen::Vector3d velocity, velocity_predict;
    Eigen::Vector3d force;
    Eigen::Matrix3d P;
    bool isFixed;
};
struct Spring
{
    Spring(Mass* m1, Mass* m2, double stiffness):m1(m1),m2(m2),k(stiffness)
    {
        rest_len = (m1->position-m2->position).norm();
        lambda = 0.0;
    }
    void resetLambda(){lambda = 0.0;}
    void resetDualGradient(){dualGradient = 0.0;}
    Mass *m1, *m2;
    double k;
    double rest_len;

    // Dual variables
    double lambda;
    double dualGradient;
};
class SpringSystem
{
public:
    SpringSystem(World* w, const Eigen::Vector3d& start, const Eigen::Vector3d& end, int numPoint,const std::vector<double> mass, double stiffness,  const std::vector<unsigned  int> fixList);
    ~SpringSystem();
    void primalUpdate();
    void dualUpdate();
    void getEdges(Eigen::MatrixXd& P1, Eigen::MatrixXd& p2);
private:
    World* world;
    std::vector<Spring*> springs;
    std::vector<Mass*> masses;
    void computeGradient();
    void computePreconditioner();
};
#endif //LEARNFEM_SPRING_H
