#include "World.h"

void World::step()
{
    softbody->update();
}
void World::setDt(double timestep)
{
    dt = timestep;
}
double World::getDt()
{
    return dt;
}
