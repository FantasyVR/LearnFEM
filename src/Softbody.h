#include "Mesh.h"
#include <vector>
class Softbody
{
private:
    Mesh* mesh;
    Eigen::VectorXd mass;
    std::vector<Eigen::Matrix3d> B;
    Eigen::VectorXd vels;
    Eigen::VectorXd gradient;
    Eigen::MatrixXd hessian;
    double energy;

    double density;
    double youngModule;
    double possionRatio;
    double mu;
    double lambda;
    void computeB();
    void computeMassMatrix();
public:
    Softbody() = default;
    Softbody(Mesh* mesh, double density, double young, double possion);
    ~Softbody();
    const Mesh* getMesh() const{return mesh;}
};