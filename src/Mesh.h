#ifndef MESH
#define MESH
#include <igl/readMSH.h>
#include <string>
class Mesh{
private:
    Eigen::MatrixXd x_pre, x_rest, x;
    Eigen::MatrixXi tets, faces;
    Eigen::VectorXd  rest_volume;
public:
    Mesh()= default;
    ~Mesh()= default;
    explicit Mesh(std::string filepath);
    Mesh(std::string filepath, Eigen::Vector3d &translation, Eigen::Vector3d &scale, Eigen::Quaterniond &orentation);
    int getNumTets() {return tets.rows();}
    const Eigen::MatrixXd& getX() const {return x;};
    Eigen::MatrixXd& getX() {return x;}
    const Eigen::MatrixXd& getPreX() const {return x_pre;};
    Eigen::MatrixXd& getPreX() {return x_pre;}
    const Eigen::MatrixXi& getTets() const {return tets;}
    Eigen::MatrixXi& getTets() {return tets;}
    const Eigen::MatrixXi& getFaces() const {return faces;}
    Eigen::MatrixXi& getFaces() {return faces;}
    const Eigen::VectorXd& getRestVolume()const {return rest_volume;}
    Eigen::VectorXd& getRestVolume() {return rest_volume;}
private:
    bool laodMesh(std::string filepath);
    void computeRestVolume();
};
#endif