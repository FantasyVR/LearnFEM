//
// Created by yupeng on 2021/2/3.
//

#include "Spring.h"
#include "World.h"

SpringSystem::SpringSystem(World* w,const Eigen::Vector3d &start, const Eigen::Vector3d &end, int numPoint,
                           const std::vector<double> mass, double stiffness, const std::vector<unsigned int> fixList)
{
    world = w;
    Eigen::Vector3d step = (end - start)/ numPoint;
    for(int i =0; i< numPoint; i++)
    {
        Mass* m = new Mass(mass[i],start+step*i,false);
        masses.push_back(m);
    }
    for(auto fix: fixList)
    {
        masses[fix]->isFixed = true;
    }

    for(int i = 0; i<numPoint-1; i++)
    {
        Spring* s = new Spring(masses[i],masses[i+1],stiffness);
        springs.push_back(s);
    }
}
SpringSystem::~SpringSystem()
{
    for(auto &s: springs)
    {
        delete s;
    }
    for(auto &m: masses)
    {
        delete m;
    }
    springs.clear();
    masses.clear();
}

void SpringSystem::update()
{
    //TODO collision detection

    // Semi-Euler Update
    double h = world->getDt();
    for(auto &m: masses)
    {
        if(m->isFixed)
            continue;
        m->velocity += h * world->getGravity();
        m->velocity_predict = m->velocity;
        m->position_pre = m->position;
        m->position += h * m->velocity;
    }
    for(int i = 0; i < 500; i++)
    {
        // Initialize force f, and Jacobian diagonal p;
        for(auto &m: masses)
        {
            m->resetForce();
            m->resetP();
        }
        // Evaluate forces and derivatives
        for(auto &s: springs)
        {
            auto m1 = s->m1;
            auto m2 = s->m2;
            auto constraint = (m1->position - m2->position).norm() - s->rest_len;
            auto dir = (m1->position - m2->position) / (m1->position - m2->position).norm();
            m1->force -= s->k * constraint * dir;
            m2->force += s->k * constraint * dir;
            // Build Pre-conditioner
            double param = h*h*s->k;
            Eigen::Vector3d Jd2 = param * dir.cwiseProduct(dir);
            Eigen::Vector3d Jd2_m1 = m1->mass * Eigen::Vector3d::Ones() + Jd2;
            Eigen::Vector3d Jd2_m2 = m2->mass * Eigen::Vector3d::Ones() - Jd2;
            Jd2_m1 = Jd2_m1.cwiseInverse();
            Jd2_m2 = Jd2_m2.cwiseInverse();
            m1->P = Eigen::DiagonalMatrix<double,3,3> (Jd2_m1);
            m2->P = Eigen::DiagonalMatrix<double,3,3> (Jd2_m2);
        }
        // Compute gradient
        for(auto &m: masses)
        {
            if(m->isFixed)
                continue;
            Eigen::Vector3d d;
            d = m->mass * (m->velocity - m->velocity_predict) - h * m->force;
            // Update state
            double alpha = 0.1; // TODO: we need to backtrace and determine step size: alpha
            m->velocity += - alpha * m->P * d;
            m->position = m->position_pre + h * m->velocity;
        }
    }
}
void SpringSystem::getEdges(Eigen::MatrixXd &P1, Eigen::MatrixXd &P2)
{
    P1.resize(springs.size(),3);
    P2.resize(springs.size(),3);
    for(int s = 0; s < springs.size(); ++s)
    {
        P1.row(s) = springs[s]->m1->position;
        P2.row(s) = springs[s]->m2->position;
    }
}