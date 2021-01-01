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
    int getNumTets() {return tets.rows();}
    const Eigen::MatrixXd& getX() const {return x;};
    Eigen::MatrixXd& getX() {return x;}
    const Eigen::MatrixXi getTets() const {return tets;}
    Eigen::MatrixXi getTets() {return tets;}
    const Eigen::MatrixXi getFaces() const {return faces;}
    Eigen::MatrixXi getFaces() {return faces;}
    const Eigen::MatrixXd getRestVolume()const {return rest_volume;}
    Eigen::MatrixXd getRestVolume() {return rest_volume;}
private:
    bool laodMesh(std::string filepath);
    void computeRestVolume();
};