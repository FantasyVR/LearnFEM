#include <igl/readMSH.h>
#include <string>
class Cloth{
private:
    Eigen::MatrixXd x_pre, x_rest,x_old, x, x_next;
    Eigen::MatrixXi faces, edges;
    Eigen::MatrixXd vel;
    Eigen::VectorXd inv_mass, rest_len;
    int num_pos{0}, num_edges{0}, num_tris{0}, N{0};
public:
    Cloth()= default;
    ~Cloth()= default;
    explicit Cloth(int N);
    const Eigen::MatrixXd& getX() const {return x;};
    Eigen::MatrixXd& getX() {return x;}
    const Eigen::MatrixXi getFaces() const {return faces;}
    Eigen::MatrixXi getFaces() {return faces;}

    void step(double h, int num_ite);

    void chebyshev_step(double h, int num_ite);

    void print_x();
private:
    void init_pos();
    void init_tri();
    void init_edge();

    void init_inv_mass();
    void init_vel();
    void init_rest_len();


    void semi_euler(double h);
    void solve_constraints();
    void update_vel(double h);
};