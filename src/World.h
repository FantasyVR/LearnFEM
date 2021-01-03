#ifndef WORLD
#define  WORLD
#include "Softbody.h"
class World
{
private:
    Softbody *softbody;
    double dt;
    Eigen::Vector3d gravity;
public:
    World();
    World(Softbody* sb):softbody(sb){}
    ~World(){delete softbody;}
    void addSoftbody(Softbody *sb){softbody = sb;}
    void step();
    void setDt(double timestep);
    const double getDt () const;
    const Eigen::Vector3d& getGravity() const;
};

#endif