#ifndef WORLD
#define  WORLD
#include "Softbody.h"
#include "Rigidbody.h"
class World
{
private:
    Softbody *softbody;
    Rigidbody *rigidbody;
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
    void step();
    void setDt(double timestep);
    const double getDt () const;
    const Eigen::Vector3d& getGravity() const;
};

#endif