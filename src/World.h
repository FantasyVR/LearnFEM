#ifndef WORLD
#define  WORLD
#include "Softbody.h"
#include "Rigidbody.h"
#include "Spring.h"
class World
{
private:
    Softbody *softbody;
    Rigidbody *rigidbody;
    SpringSystem* springsystem;
    double dt;
    Eigen::Vector3d gravity;
public:
    World();
    World(Softbody* sb):softbody(sb){}
    ~World(){delete softbody;}
    void addSoftbody(Softbody *sb){softbody = sb;}
    Softbody* getSoftbody(){return softbody;}
    void addRigidbody(Rigidbody* rb){rigidbody = rb;}
    Rigidbody* getRigidbody(){return rigidbody;}
    void addSpringsystem(SpringSystem* ss){ springsystem= ss;}
    SpringSystem* getSpringsystem(){return springsystem;}
    void step();
    void setDt(double timestep);
    const double getDt () const;
    void setGravity(const Eigen::Vector3d &gv);
    const Eigen::Vector3d& getGravity() const;
};

#endif