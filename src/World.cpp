#include "World.h"
World::World()
{
    dt = 0.025;
    gravity = Eigen::Vector3d(0.0,0.0,0.0);
}
void World::step()
{
//    if(softbody)
//        softbody->update();
//    if(rigidbody)
//        rigidbody->update();
    if(springsystem)
        springsystem->update();
}
void World::setDt(double timestep)
{
    dt = timestep;
}

const double World::getDt() const {
    return dt;
}

const Eigen::Vector3d &World::getGravity() const{
    return gravity;
}
void World::setGravity(const Eigen::Vector3d &gv)
{
    gravity = gv;
}


