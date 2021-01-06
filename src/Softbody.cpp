#include "Softbody.h"
#include "World.h"
#include "SVD/ImplicitQRSVD.h"

#include <iostream>
Softbody::Softbody(World* world,Mesh *mesh,double density, double YM, double PR)
{
    this->world = world;
    this->mesh = mesh;
    this->youngModule = YM;
    this->possionRatio = PR;
    this->density = density;
    this->mu = YM / 2.0 / (1.0 + PR);
    this->lambda = YM * PR / (1.0 + PR) / (1.0 - 2.0 * PR);
    int numVerts = mesh->getX().rows();
    vels.setZero(numVerts * 3);
    gradient.setZero(numVerts * 3);
    hessian.setZero(numVerts*3, numVerts*3);
    computeB();
    computeMassMatrix();
}

void Softbody::computeB() {
    int numTets = mesh->getNumTets();
    const Eigen::MatrixXd &positions = mesh->getX();
    const Eigen::MatrixXi &tets = mesh->getTets();
    for(int i = 0; i<numTets; i++)
    {
        auto tetIdx = tets.row(i);
        const Eigen::Vector3d& p0 = positions.row(tetIdx[0]);
        const Eigen::Vector3d& p1 = positions.row(tetIdx[1]);
        const Eigen::Vector3d& p2 = positions.row(tetIdx[2]);
        const Eigen::Vector3d& p3 = positions.row(tetIdx[3]);
        Eigen::Matrix3d D;
        D.col(0) = p1 - p0;
        D.col(1) = p2 - p0;
        D.col(2) = p3 - p0;
        B.emplace_back(D.inverse());
    }
}

void Softbody::computeMassMatrix() {
    const Eigen::MatrixXd &positions = mesh->getX();
    const Eigen::VectorXd &rest_volumes = mesh->getRestVolume();
    const Eigen::MatrixXi &tets = mesh->getTets();
    mass.setZero(positions.rows());
    for(int i = 0; i<tets.rows(); i++)
    {
         auto tetIdx = tets.row(i);
         auto rest_volume = rest_volumes(i);
         mass(tetIdx[0]) += rest_volume/4.0;
         mass(tetIdx[1]) += rest_volume/4.0;
         mass(tetIdx[2]) += rest_volume/4.0;
         mass(tetIdx[3]) += rest_volume/4.0;
    }
    mass *= density;
    std::cout<<"Mass Matrix: \n"<< mass <<std::endl;
}

Softbody::~Softbody() {
    delete mesh;
}
void Softbody::computeEnergy() {
    int numTets = mesh->getNumTets();
    energy = 0;
    for(int tet =0; tet<numTets; tet++)
    {
        const double sigmam12Sum =    (Sigmav[tet] - Eigen::Vector3d::Ones()).squaredNorm();
        const double sigmaProdm1 = Sigmav[tet].prod() - 1.0;
        double e = mu * sigmam12Sum + lambda / 2.0 * sigmaProdm1 * sigmaProdm1;
        double dt = world->getDt();
        e *= mesh->getRestVolume()[tet] * dt * dt;
        energy += e;
    }
    //std::cout<<"Sum Energy: "<< energy <<std::endl;
}

void Softbody::computeInternalForce() {
    int numTets = mesh->getNumTets();
    const Eigen::MatrixXd &positions = mesh->getX();
    const Eigen::MatrixXi &tets = mesh->getTets();
    gradient.setZero(positions.rows() * 3);
    for (int tet =0; tet < numTets; tet++)
    {
        Eigen::Matrix3d P;
        compute_dE_div_dF(F[tet],Uv[tet],Vv[tet],Sigmav[tet],P);
        double dt = world->getDt();
        double w = dt * dt * mesh->getRestVolume()[tet];
        P *= w;
        Eigen::Matrix3d H = P * B[tet].transpose();
        Eigen::Vector3d f1 = H.col(0); int sf1 = tets(tet,1) * 3; gradient.segment(sf1,3) += f1;
        Eigen::Vector3d f2 = H.col(1); int sf2 = tets(tet,2) * 3; gradient.segment(sf2,3) += f2;
        Eigen::Vector3d f3 = H.col(2); int sf3 = tets(tet,3) * 3; gradient.segment(sf3,3) += f3;
        Eigen::Vector3d f0 = -f1 -f2 -f3; int sf0 = tets(tet,0)* 3; gradient.segment(sf0,3) += f0;
    }
    //std::cout<<"internal force: \n" << gradient <<std::endl;
}

void Softbody::computeStiffnessMatrixByEle(Eigen::Matrix<double, 12, 12> &hessian, const Eigen::Matrix3d &U,
                                           const Eigen::Matrix3d &V, const Eigen::Vector3d &sigma,
                                           const Eigen::Matrix3d &A, const double rVolume)
{
    Eigen::Vector3d dE_div_dsigma;
    compute_dE_div_dsigma(sigma,dE_div_dsigma);
    Eigen::Matrix3d d2E_divdsigma2;
    compute_d2E_div_dsigma2(sigma, d2E_divdsigma2);
    makePD(d2E_divdsigma2);

    Eigen::Vector3d BLeftCoef;
    compute_BLeftCoef(sigma, BLeftCoef);

    Eigen::Matrix2d B[3];
    for (int cI = 0; cI < 3; cI++) {
        int cI_post = (cI + 1) % 3;

        double rightCoef = dE_div_dsigma[cI] + dE_div_dsigma[cI_post];
        double sum_sigma = sigma[cI] + sigma[cI_post];
        const double eps = 1.0e-6;
        if (sum_sigma < eps) {
            rightCoef /= 2.0 * eps;
        }
        else {
            rightCoef /= 2.0 * sum_sigma;
        }
        const double& leftCoef = BLeftCoef[cI];
        B[cI](0, 0) = B[cI](1, 1) = leftCoef + rightCoef;
        B[cI](0, 1) = B[cI](1, 0) = leftCoef - rightCoef;
        makePD2d(B[cI]);
    }
    double w = world->getDt() * world->getDt() * rVolume;
    Eigen::Matrix<double,9,9> M;
    M.setZero();
    M(0,0) = w * d2E_divdsigma2(0,0);
    M(0,4) = w * d2E_divdsigma2(0,1);
    M(0,8) = w * d2E_divdsigma2(0,2);
    M(4,0) = w * d2E_divdsigma2(1,0);
    M(4,4) = w * d2E_divdsigma2(1,1);
    M(4,8) = w * d2E_divdsigma2(1,2);
    M(8,0) = w * d2E_divdsigma2(2,0);
    M(8,4) = w * d2E_divdsigma2(2,1);
    M(8,8) = w * d2E_divdsigma2(2,2);
    // B01
    M(1, 1) = w * B[0](0, 0);
    M(1, 3) = w * B[0](0, 1);
    M(3, 1) = w * B[0](1, 0);
    M(3, 3) = w * B[0](1, 1);
    // B12
    M(5, 5) = w * B[1](0, 0);
    M(5, 7) = w * B[1](0, 1);
    M(7, 5) = w * B[1](1, 0);
    M(7, 7) = w * B[1](1, 1);
    // B20
    M(2, 2) = w * B[2](1, 1);
    M(2, 6) = w * B[2](1, 0);
    M(6, 2) = w * B[2](0, 1);
    M(6, 6) = w * B[2](0, 0);

    Eigen::Matrix<double, 9, 9> wdP_div_dF;
    for (int i = 0; i < 3; i++) {
        int _dim_i = i * 3;
        for (int j = 0; j < 3; j++) {
            int ij = _dim_i + j;
            for (int r = 0; r < 3; r++) {
                int _dim_r = r * 3;
                for (int s = 0; s < 3; s++) {
                    int rs = _dim_r + s;
                    if (ij > rs) {
                        // bottom left, same as upper right
                        continue;
                    }
                    wdP_div_dF(ij, rs) = M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0) + M(0, 4) * U(i, 0) * V(j, 0) * U(r, 1) * V(s, 1) + M(0, 8) * U(i, 0) * V(j, 0) * U(r, 2) * V(s, 2) + M(4, 0) * U(i, 1) * V(j, 1) * U(r, 0) * V(s, 0) + M(4, 4) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1) + M(4, 8) * U(i, 1) * V(j, 1) * U(r, 2) * V(s, 2) + M(8, 0) * U(i, 2) * V(j, 2) * U(r, 0) * V(s, 0) + M(8, 4) * U(i, 2) * V(j, 2) * U(r, 1) * V(s, 1) + M(8, 8) * U(i, 2) * V(j, 2) * U(r, 2) * V(s, 2) + M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1) + M(1, 3) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0) + M(3, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1) + M(3, 3) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0) + M(5, 5) * U(i, 1) * V(j, 2) * U(r, 1) * V(s, 2) + M(5, 7) * U(i, 1) * V(j, 2) * U(r, 2) * V(s, 1) + M(7, 5) * U(i, 2) * V(j, 1) * U(r, 1) * V(s, 2) + M(7, 7) * U(i, 2) * V(j, 1) * U(r, 2) * V(s, 1) + M(2, 2) * U(i, 0) * V(j, 2) * U(r, 0) * V(s, 2) + M(2, 6) * U(i, 0) * V(j, 2) * U(r, 2) * V(s, 0) + M(6, 2) * U(i, 2) * V(j, 0) * U(r, 0) * V(s, 2) + M(6, 6) * U(i, 2) * V(j, 0) * U(r, 2) * V(s, 0);
                    if (ij < rs) {
                        wdP_div_dF(rs, ij) = wdP_div_dF(ij, rs);
                    }
                }
            }
        }
    }

    Eigen::Matrix<double, 12, 9> wdP_div_dx;
    dF_div_dx_mult<9>(wdP_div_dF.transpose(), A, wdP_div_dx);
    dF_div_dx_mult<12>(wdP_div_dx.transpose(), A, hessian);
}

void Softbody::computeStifnessMatrix()
{
    int numTets = mesh->getNumTets();
    const Eigen::MatrixXd &positions = mesh->getX();
    const Eigen::MatrixXi &tets = mesh->getTets();
    int numVert = positions.rows();
    hessian.resize(numVert*3, numVert*3);
    for(int tet=0; tet< numTets; tet++)
    {
        auto idx = tets.row(tet);
        Eigen::Matrix<double, 12,12> hessianEle;
        computeStiffnessMatrixByEle(hessianEle, Uv[tet], Vv[tet], Sigmav[tet],B[tet],mesh->getRestVolume()[tet]);
        for(int i =0; i<4;i++)
            for(int j = 0; j<4;j++)
            {
                int startRow = i * 3;
                int startCol = j * 3;
                int desRow = idx[i] * 3;
                int desCol = idx[j] * 3;
                hessian.block(desRow,desCol,3,3) = hessianEle.block(startRow,startCol,3,3);
            }
    }
}
void Softbody::assambleGlobalMatrix() {
    auto &A = hessian;
    int numVert = mesh->getX().rows();
    for(int vert = 0; vert < numVert; vert++)
    {
        int start = vert * 3;
        A(start + 0, start + 0) += mass[vert];
        A(start + 1, start + 1) += mass[vert];
        A(start + 2, start + 2) += mass[vert];
    }
}

void Softbody::computeDeformationGradient() {
    int numTets = mesh->getNumTets();
    const Eigen::MatrixXd &positions = mesh->getX();
    const Eigen::MatrixXi &tets = mesh->getTets();
    F.resize(numTets);
    Uv.resize(numTets);
    Vv.resize(numTets);
    Sigmav.resize(numTets);
    for(int i = 0; i<numTets; i++)
    {
        auto tetIdx = tets.row(i);
        const Eigen::Vector3d& p0 = positions.row(tetIdx[0]);
        const Eigen::Vector3d& p1 = positions.row(tetIdx[1]);
        const Eigen::Vector3d& p2 = positions.row(tetIdx[2]);
        const Eigen::Vector3d& p3 = positions.row(tetIdx[3]);
        Eigen::Matrix3d D;
        D.col(0) = p1 - p0;
        D.col(1) = p2 - p0;
        D.col(2) = p3 - p0;
        Eigen::Matrix3d dg;
        dg = D * B[i];   // deformation gradient
        F[i] = dg;
        Eigen::Vector3d S;
        Eigen::Matrix3d U;
        Eigen::Matrix3d V;
        JIXIE::singularValueDecomposition(dg,U,S,V);
        Uv[i] = U;
        Vv[i] = V;
        Sigmav[i] = S;
    }
}

void Softbody::compute_dE_div_dF(const Eigen::Matrix3d &F, const Eigen::Matrix3d& U, const Eigen::Matrix3d V, const Eigen::Vector3d &sigma, Eigen::Matrix3d &P) {
    P.setZero();
    Eigen::Matrix3d JFInvT;
    JFInvT(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
    JFInvT(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
    JFInvT(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
    JFInvT(1, 0) = F(0, 2) * F(2, 1) - F(0, 1) * F(2, 2);
    JFInvT(1, 1) = F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0);
    JFInvT(1, 2) = F(0, 1) * F(2, 0) - F(0, 0) * F(2, 1);
    JFInvT(2, 0) = F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1);
    JFInvT(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
    JFInvT(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
    P = mu * 2 * (F - U * V.transpose()) + lambda * (sigma.prod() - 1) * JFInvT;
}

void Softbody::update() {
    int maxIte = 1;
    Eigen::MatrixXd &positions = mesh->getX();
    Eigen::MatrixXd tmpPos = positions;
    int numVert = positions.rows();
    Eigen::MatrixXd xTileta;
    computeXTilta(xTileta);
    for(int ite = 0; ite < maxIte; ite++)
    {
        // 1. initX  warm start searchdir
        Eigen::VectorXd searchDir;
        searchDir.setZero(numVert * 3);
        int stepsize = 1.0;
        // 2. step forward
        for(int i = 0; i < numVert; i++)
            positions.row(i) += stepsize * searchDir.segment<3>(i * 3).transpose();
        // 3. computeEnergyVal
        computeDeformationGradient();
        computeEnergy();
        double initEnergy = energy;
        std::cout<<"init energy: "<< initEnergy;
        // 4. compute gradient
        computeInternalForce();
        // gradient += m * v;
        for(int i = 0; i < numVert; i++)
            gradient.segment<3>(i * 3) += mass[i] * (positions.row(i) - xTileta.row(i)).transpose();
        std::cout<<" ||g||^2 = " << gradient.squaredNorm() <<std::endl;
        // 5. implicit solve iteration
        int iterCap = 1000, currentIter = 0;
        do {
            computeStifnessMatrix(); // K
            assambleGlobalMatrix(); // H = K + M
            searchDir = hessian.ldlt().solve(-gradient); //solve: H * searchdir = -gradient
            std::cout << "stepLen = " << (stepsize * searchDir).squaredNorm() << std::endl;
            // TO: Backtrace Line Search
            // step forward
            for(int i = 0; i < numVert; i++)
                positions.row(i) += stepsize * searchDir.segment<3>(i * 3).transpose();
            // Re-compute Energy and Gradient
            computeDeformationGradient();
            computeEnergy();
            double currentEnergy = energy;
            std::cout<<"E_cur_smooth: "<< currentEnergy <<std::endl;

            // gradient += m * v;
            computeInternalForce();
            for(int i = 0; i < numVert; i++)
                gradient.segment<3>(i * 3) += mass[i] * (positions.row(i) - xTileta.row(i)).transpose();
            std::cout<<"||gradient||^2 = " << gradient.squaredNorm() << std::endl;

            if(++currentIter > iterCap)
                break;
        }while(gradient.squaredNorm() > 0.000016);

        //update positions and velocities
        for(int i = 0; i < numVert; i++)
            vels.segment<3>(i * 3) = (positions.row(i) - tmpPos.row(i)).transpose()/world->getDt();
        tmpPos = positions;
    }
}

void Softbody::compute_dE_div_dsigma(const Eigen::Vector3d &sigma, Eigen::Vector3d &dE_div_dsigma) {
    double sigmaProdm1lambda  = lambda * (sigma.prod() - 1.0);
    double sigmaProd_noI_0 = sigma[1] * sigma[2];
    double sigmaProd_noI_1 = sigma[2] * sigma[0];
    double sigmaProd_noI_2 = sigma[0] * sigma[1];

    dE_div_dsigma.setZero();
    dE_div_dsigma[0] = 2 * mu * (sigma[0] - 1.0) + sigmaProd_noI_0 * sigmaProdm1lambda;
    dE_div_dsigma[1] = 2 * mu * (sigma[1] - 1.0) + sigmaProd_noI_1 * sigmaProdm1lambda;
    dE_div_dsigma[2] = 2 * mu * (sigma[2] - 1.0) + sigmaProd_noI_2 * sigmaProdm1lambda;
}

void Softbody::compute_d2E_div_dsigma2(const Eigen::Vector3d &sigma, Eigen::Matrix3d &d2E_div_dsigma2) {
    d2E_div_dsigma2.setZero();
    double sigmaProd_noI_0 = sigma[1] * sigma[2];
    double sigmaProd_noI_1 = sigma[2] * sigma[0];
    double sigmaProd_noI_2 = sigma[0] * sigma[1];
    d2E_div_dsigma2(0,0) = 2 * mu + lambda * sigmaProd_noI_0 * sigmaProd_noI_0;
    d2E_div_dsigma2(1,1) = 2 * mu + lambda * sigmaProd_noI_1 * sigmaProd_noI_1;
    d2E_div_dsigma2(2,2) = 2 * mu + lambda * sigmaProd_noI_2 * sigmaProd_noI_2;

    d2E_div_dsigma2(0,1) = d2E_div_dsigma2(1,0) = lambda * (sigma[2] * (sigma.prod() - 1.0) + sigmaProd_noI_0 * sigmaProd_noI_1);
    d2E_div_dsigma2(0,2) = d2E_div_dsigma2(2,0) = lambda * (sigma[1] * (sigma.prod() - 1.0) + sigmaProd_noI_0 * sigmaProd_noI_2);
    d2E_div_dsigma2(1,2) = d2E_div_dsigma2(2,1) = lambda * (sigma[0] * (sigma.prod() - 1.0) + sigmaProd_noI_2 * sigmaProd_noI_1);
}

void Softbody::makePD(Eigen::Matrix3d &A) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(A);
    if (eigenSolver.eigenvalues()[0] >= 0.0) {
        return;
    }
    Eigen::DiagonalMatrix<double, 3> D(eigenSolver.eigenvalues());
    for (int i = 0; i < 3; i++) {
        if (D.diagonal()[i] < 0.0) {
            D.diagonal()[i] = 0.0;
        }
        else {
            break;
        }
    }
    A = eigenSolver.eigenvectors() * D * eigenSolver.eigenvectors().transpose();
}

void Softbody::compute_BLeftCoef(const Eigen::Vector3d &sigma, Eigen::Vector3d &BLeftCoef) {
    BLeftCoef.setZero();
    double halfLambda = lambda * 0.5;
    BLeftCoef(0) = mu - halfLambda * sigma[2] * (sigma.prod() - 1.0);
    BLeftCoef(1) = mu - halfLambda * sigma[1] * (sigma.prod() - 1.0);
    BLeftCoef(2) = mu - halfLambda * sigma[0] * (sigma.prod() - 1.0);
}

void Softbody::makePD2d(Eigen::Matrix2d &symMtr) {
    const double a = symMtr(0, 0);
    const double b = (symMtr(0, 1) + symMtr(1, 0)) / 2.0;
    const double d = symMtr(1, 1);

    double b2 = b * b;
    const double D = a * d - b2;
    const double T_div_2 = (a + d) / 2.0;
    const double sqrtTT4D = std::sqrt(T_div_2 * T_div_2 - D);
    const double L2 = T_div_2 - sqrtTT4D;
    if (L2 < 0.0) {
        const double L1 = T_div_2 + sqrtTT4D;
        if (L1 <= 0.0) {
            symMtr.setZero();
        }
        else {
            if (b2 == 0.0) {
                symMtr << L1, 0.0, 0.0, 0.0;
            }
            else {
                const double L1md = L1 - d;
                const double L1md_div_L1 = L1md / L1;
                symMtr(0, 0) = L1md_div_L1 * L1md;
                symMtr(0, 1) = symMtr(1, 0) = b * L1md_div_L1;
                symMtr(1, 1) = b2 / L1;
            }
        }
    }
}


template <int colSize>
void Softbody::dF_div_dx_mult(const Eigen::Matrix<double, 9 , colSize>& right,
                           const Eigen::Matrix<double, 3, 3>& A,
                           Eigen::Matrix<double, 12, colSize>& result)
{
    if (colSize == Eigen::Dynamic) {
        assert(right.cols() > 0);
    }
    else {
        assert(colSize > 0);
    }
    for (int colI = 0; colI < right.cols(); colI++) {
        result(3, colI)  = (A.row(0) * right.block(0, colI, 3, 1))[0];
        result(4, colI)  = (A.row(0) * right.block(3, colI, 3, 1))[0];
        result(5, colI)  = (A.row(0) * right.block(6, colI, 3, 1))[0];
        result(6, colI)  = (A.row(1) * right.block(0, colI, 3, 1))[0];
        result(7, colI)  = (A.row(1) * right.block(3, colI, 3, 1))[0];
        result(8, colI)  = (A.row(1) * right.block(6, colI, 3, 1))[0];
        result(9, colI)  = (A.row(2) * right.block(0, colI, 3, 1))[0];
        result(10, colI) = (A.row(2) * right.block(3, colI, 3, 1))[0];
        result(11, colI) = (A.row(2) * right.block(6, colI, 3, 1))[0];
        result(0, colI) = -result(3, colI) - result(6, colI) - result(9, colI);
        result(1, colI) = -result(4, colI) - result(7, colI) - result(10, colI);
        result(2, colI) = -result(5, colI) - result(8, colI) - result(11, colI);
    }
}

void Softbody::computeXTilta(Eigen::MatrixXd &xTileta) {
    Eigen::MatrixXd& positions = mesh->getX();
    Eigen::MatrixXd& pre_pos = mesh->getPreX();
    xTileta.resize(positions.rows(),positions.cols());
    int numVert = positions.rows();
    double dt = world->getDt();
    Eigen::Vector3d gravity = world->getGravity();
    for(int i =0 ; i<numVert; i++) {
        Eigen::Vector3d vel = vels.segment(i * 3, 3);
        xTileta.row(i) = positions.row(i) + dt * vel.transpose() + dt * dt * gravity.transpose();
        vels.segment(i * 3, 3)  = (positions.row(i) - pre_pos.row(i)).transpose()/dt;
    }
}