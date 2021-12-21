#include "Cloth.h"
#include <iostream>
#include <fstream>

Cloth::Cloth(int N):N(N){
    num_pos = (N+1)*(N+1);
    num_tris = 2 * N * N;
    num_edges = 2 * N * (N+1) + N * N;
    init_pos();
    init_tri();
    init_edge();

    init_inv_mass();
    init_vel();
    init_rest_len();
}

void Cloth::init_pos(){
    x.resize(num_pos, 3);
    for(int i = 0; i < N+1; i++)
        for(int j = 0; j < N+1; j++){
            int idx = i * (N+1) + j;
            x.row(idx) << double(i) / N, 0.0, double(j) / N;
        }
    x_pre = x;
    x_rest = x;
    x_next = x;
}

void Cloth::init_tri(){
    faces.resize(num_tris, 3);
    for(int i = 0; i < N; i++)
        for(int j =0; j < N; j++){
            int tri_idx = 2 * (i * N + j);
            int pos_idx = i * (N + 1) + j;
            if((i+j) % 2 == 0){
                faces.row(tri_idx) << pos_idx, pos_idx + N + 2, pos_idx + 1;
                faces.row(tri_idx + 1) << pos_idx, pos_idx + N + 1, pos_idx + N+2;
            }
            else{
                faces.row(tri_idx) << pos_idx, pos_idx + N + 1, pos_idx + 1;
                faces.row(tri_idx + 1) << pos_idx+1, pos_idx + N + 1, pos_idx + N + 2;
            }
        }
}

void Cloth::init_edge(){
    edges.resize(num_edges, 2);
    for(int i = 0; i < N+1; i++)
        for(int j =0; j < N; j++){
            int edge_idx = i * N + j;
            int pos_idx = i *(N+1)+j;
            edges.row(edge_idx) <<  pos_idx, pos_idx+1;
        }
    int start = N * (N+1);
    for(int i = 0; i < N; i++)
        for(int j =0; j < N+1; j++){
            int edge_idx = start +  j * N + i;
            int pos_idx = i *(N+1)+j;
            edges.row(edge_idx) << pos_idx, pos_idx + N + 1;
        }
    start = 2 * N * (N+1);
    for(int i = 0; i < N; i++)
        for(int j =0; j < N; j++){
            int edge_idx = start + i * N + j;
            int pos_idx = i * (N+1) +j;
            if((i+j) % 2 == 0){
                edges.row(edge_idx) << pos_idx, pos_idx + N + 2;
            }else{
                edges.row(edge_idx) << pos_idx + 1, pos_idx + N + 1;
            }
        }
}

void Cloth::print_x() {
    for(int i = 0; i < num_pos; i++)
        std::cout<< x.row(i) <<std::endl;
}

void Cloth::init_inv_mass() {
    inv_mass.resize(num_pos);
    for(int i = 0; i < num_pos; i++)
        inv_mass(i) = 1.0;
    inv_mass(0) = 0.0;
    inv_mass(N) = 0.0;
}

void Cloth::init_vel() {
    vel.resize(num_pos, 3);
    vel.setZero();
}

void Cloth::init_rest_len() {
    rest_len.resize(num_edges);
    for(int i =0; i < num_edges; i++){
        const Eigen::Vector2i &edge = edges.row(i);
        const Eigen::Vector3d& p1 = x.row(edge(0));
        const Eigen::Vector3d& p2 = x.row(edge(1));
        rest_len(i) = (p1-p2).norm();
    }
}

void Cloth::step(double h, int num_ite) {
    semi_euler(h);
    for(int i = 0; i < num_ite; i++)
    {
        solve_constraints();

        x_pre = x;
    }
    update_vel(h);
}

void Cloth::chebyshev_step(double h, int num_ite)
{
    //char* output_file = "../../statistics/pbd_primal_chebyshev_dual_residual.txt";
    char* output_file = "../../statistics/pbd_dual_residual.txt";
    std::fstream residual_file(output_file, std::ios_base::app);
    if (!residual_file.is_open()) {
        printf("Failed to open residual files\n");
    }
    // semi-euler integration
    for (int i = 0; i < num_pos; i++) {
        if (inv_mass(i) == 0.0) continue;
        vel.row(i) += h * Eigen::Vector3d(0.0, -9.8, 0.0);
        x.row(i) += h * vel.row(i);
    }
    x_old = x;
    for (int ite = 0; ite < num_ite; ite++) {
        x_next = x;
        double dual_residual = 0;
        // solve constraints
        for (int i = 0; i < num_edges; i++) {
            const Eigen::Vector2i& edge = edges.row(i);
            int idx0 = edge(0);
            int idx1 = edge(1);
            const Eigen::Vector3d& p0 = x_next.row(idx0);
            const Eigen::Vector3d& p1 = x_next.row(idx1);
            double inv_M0 = inv_mass(idx0);
            double inv_M1 = inv_mass(idx1);
            Eigen::Vector3d dis = p0 - p1;
            Eigen::Vector3d gradient = dis.normalized();
            double constraint = dis.norm() - rest_len(i);
            dual_residual += constraint;
            double l = -constraint / (inv_M0 + inv_M1);
            if (inv_M0 != 0.0)
                x_next.row(idx0) += inv_M0 * l * gradient;
            if (inv_M1 != 0.0)
                x_next.row(idx1) -= inv_M1 * l * gradient;
        }
        residual_file << dual_residual << "\n";

        double omega = 1.0, rho = 0.9, under_relax = 1.0;
        bool primal_chebyshev = false;
        if (primal_chebyshev)
        {
            if (ite <= 10) omega = 1.0;
            else if (ite == 11) omega = 2 / (2 - rho * rho);
            else omega = 4 / (4 - rho * rho);
            for (int i = 0; i < num_pos; i++)
            {
                //x_next.row(i) = (x_next.row(i) - x.row(i)) * under_relax;
                x_next.row(i) = omega * (x_next.row(i) - x_pre.row(i)) + x_pre.row(i);
            }
        }
        x_pre = x;
        x = x_next;
    }
    // update velocities
    for (int i = 0; i < num_pos; i++) {
        if (inv_mass(i) == 0.0) continue;
        vel.row(i) += (x.row(i) - x_old.row(i)) / h;
    }
    residual_file.close();
}


void Cloth::semi_euler(double h){
    x_old = x;
    for(int i = 0; i < num_pos; i++){
        if(inv_mass(i) == 0.0) continue;
        vel.row(i) +=  h * Eigen::Vector3d(0.0, -9.8, 0.0);
        x.row(i) += h * vel.row(i);
    }
}
void Cloth::solve_constraints() {
    for (int i = 0; i < num_edges; i++) {
        const Eigen::Vector2i &edge = edges.row(i);
        int idx0 = edge(0);
        int idx1 = edge(1);
        Eigen::Vector3d p0 = x.row(idx0);
        Eigen::Vector3d p1 = x.row(idx1);
        double inv_M0 = inv_mass(idx0);
        double inv_M1 = inv_mass(idx1);
        Eigen::Vector3d dis = p0 - p1;
        Eigen::Vector3d gradient = dis.normalized();
        double constraint = dis.norm() - rest_len(i);
        double l = -constraint / (inv_M0 + inv_M1);
        if(inv_M0 != 0.0)
            x.row(idx0) += inv_M0 * l * gradient;
        if(inv_M1 != 0.0)
            x.row(idx1) -= inv_M1 * l * gradient;
    }
}
void Cloth::update_vel(double h){
    for(int i = 0; i < num_pos; i++){
        if(inv_mass(i) == 0.0) continue;
        vel.row(i) =  (x.row(i) - x_old.row(i))/h;
    }
}