#ifndef SOFTBODY
#define SOFTBODY

#include "Mesh.h"

#include <vector>
class World;
class Softbody
{
private:
    World* world;
    Mesh* mesh;
    Eigen::VectorXd mass;
    std::vector<Eigen::Matrix3d> B;
    std::vector<Eigen::Matrix3d> F, Uv,Vv;
    std::vector<Eigen::Vector3d> Sigmav;
    Eigen::VectorXd vels;
    Eigen::VectorXd gradient;
    Eigen::MatrixXd hessian;
    double energy;

    double density;
    double youngModule;
    double possionRatio;
    double mu;
    double lambda;

public:
    Softbody() = default;
    Softbody(World* world,Mesh* mesh, double density, double young, double possion);
    ~Softbody();
    void update();
    void computeDeformationGradient();
    // compute corotaional energy
    void computeEnergy();
    void computeInternalForce();
    void computeStifnessMatrix();
    void computeStiffnessMatrixByEle(Eigen::Matrix<double,12,12> &hessian,const Eigen::Matrix3d& U,const Eigen::Matrix3d& V, const Eigen::Vector3d &sigma, const Eigen::Matrix3d &A, const double rVolume);
    void assambleGlobalMatrix();

    const Mesh* getMesh() const{return mesh;}

private:
    void computeB();
    void compute_dE_div_dF(const Eigen::Matrix3d &F, const Eigen::Matrix3d& U, const Eigen::Matrix3d V, const Eigen::Vector3d &sigma, Eigen::Matrix3d &P);
    void computeMassMatrix();
    void compute_dE_div_dsigma(const Eigen::Vector3d& sigma, Eigen::Vector3d& dE_div_dsigma);
    void compute_d2E_div_dsigma2(const Eigen::Vector3d& sigma, Eigen::Matrix3d& d2E_div_dsigma2);
    void makePD(Eigen::Matrix3d &A);
    void makePD2d(Eigen::Matrix2d &A);
    void compute_BLeftCoef(const Eigen::Vector3d& sigma, Eigen::Vector3d &BLeftCoef);
    template <int colSize>
    void dF_div_dx_mult(const Eigen::Matrix<double, 9 , colSize>& right,
                        const Eigen::Matrix<double, 3, 3>& A,
                        Eigen::Matrix<double, 12, colSize>& result);
};
#endif