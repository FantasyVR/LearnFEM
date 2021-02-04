#include "Mesh.h"

#include <igl/readOBJ.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>

Mesh::Mesh(std::string filepath)
{
    bool success = loadMsh(filepath);
    if (success)
        std::cout<<"Load MSH file successfully"<<std::endl;
    else
        std::cout<< "Failed to load MSH file" <<std::endl;
    computeRestVolume();
    x_rest = x;
    x_pre = x;
}

bool Mesh::loadMsh(const std::string &filepath)
{
    x.resize(0,3);
    tets.resize(0,4);
    faces.resize(0,3);
    //std::string filepath = std::string(RESOURCE) + "tet.msh";
    std::filesystem::path p = filepath;
    if (p.extension() != ".msh")
    {
        std::cout<<"Not .msh file"<<std::endl;
        return false;
    }
    std::fstream f1(filepath, std::ios::in);
    if (!f1.is_open())
    {
        std::cout<<"faild to open the file" <<std::endl;
        return false;
    }
    std::string line;
    int vAmt = 0;
    while (!f1.eof() && std::getline(f1,line))
    {
        if (line == "$Nodes")
        {
            std::getline(f1, line);
            std::stringstream s(line);
            int a;
            s >> a >> vAmt;
            std::getline(f1, line);
            break;
        }
    }
    x.resize(vAmt,3);
    x_rest.resize(vAmt,3);
    x_pre.resize(vAmt,3);
    for(int i = 0; i < vAmt; i++)
    {
        std::getline(f1, line);
        std::stringstream s(line);
        int idx;
        double x0, x1, x2;
        s >> idx >> x0 >> x1 >> x2;
        x(i,0) = x0;
        x(i,1) = x1;
        x(i,2) = x2;
    }
    std::cout<<"Postions: \n"<< x <<std::endl;
    while (!f1.eof() && std::getline(f1,line))
    {
        if (line == "$Elements")
        {
            std::getline(f1, line);
            std::stringstream s(line);
            int a;
            s >> a >> vAmt;
            std::getline(f1, line);
            break;
        }
    }
    tets.resize(vAmt,4);
    for(int i = 0; i < vAmt; i++)
    {
        std::getline(f1, line);
        std::stringstream s(line);
        int idx;
        int x0, x1, x2, x3;
        s >> idx >> x0 >> x1 >> x2 >> x3;
        tets(i,0) = x0;
        tets(i,1) = x1;
        tets(i,2) = x2;
        tets(i,3) = x3;
    }
    tets.array() -= 1;
    std::cout<<"Elements: \n"<< tets <<std::endl;
    while (!f1.eof() && std::getline(f1,line))
    {
        if (line == "$Surface")
        {
            std::getline(f1, line);
            std::stringstream s(line);
            s >>  vAmt;
            break;
        }
    }
    faces.resize(vAmt,3);
    for(int i = 0; i < vAmt; i++)
    {
        std::getline(f1, line);
        std::stringstream s(line);
        int x0, x1, x2;
        s >> x0 >> x1 >> x2;
        faces(i,0) = x0;
        faces(i,1) = x1;
        faces(i,2) = x2;
    }
    faces.array() -= 1;
    std::cout<<"Faces: \n" << faces <<std::endl;
    f1.close();
    return true;
}
bool Mesh::loadObj(const std::string &filepath)
{
    return igl::readOBJ(filepath, x, faces);
}
void Mesh::computeRestVolume()
{
    int numTets = tets.rows();
    rest_volume.resize(numTets);
    for(int i = 0; i< numTets; i++)
    {
        auto tetIdx = tets.row(i);
        const Eigen::Vector3d& p0 = x.row(tetIdx[0]);
        const Eigen::Vector3d& p1 = x.row(tetIdx[1]);
        const Eigen::Vector3d& p2 = x.row(tetIdx[2]);
        const Eigen::Vector3d& p3 = x.row(tetIdx[3]);
        Eigen::Matrix3d D;
        D.col(0) = p1 - p0;
        D.col(1) = p2 - p0;
        D.col(2) = p3 - p1;
        rest_volume[i] =(double)(1.0/6) *  fabs(D.determinant());
    }
}

Mesh::Mesh(std::string filepath, Eigen::Vector3d &translation, Eigen::Vector3d &scale, Eigen::Quaterniond &orentation) {
    std::string filetype = filepath.substr(filepath.find_last_of('.'));
    if (filetype == ".msh")
    {
        bool success = loadMsh(filepath);
        if (success)
            std::cout<<"Load MSH file successfully"<<std::endl;
        else
            std::cout<< "Failed to load MSH file" <<std::endl;
    }
    else if (filetype == ".obj")
    {
        bool success = loadObj(filepath);
        if (success)
            std::cout<<"Load OBJ file successfully"<<std::endl;
        else
            std::cout<< "Failed to load OBJ file" <<std::endl;
    }

    int numVert = x.rows();
    // Apply scale, rotation, translation
    for(int i = 0; i < numVert; i++)
    {
        x.row(i).cwiseProduct(scale.transpose());
        x.row(i) = orentation.toRotationMatrix() * x.row(i).transpose();
        x.row(i) += translation;
    }

    x_rest = x;
    x_pre = x;
    computeRestVolume();
}
