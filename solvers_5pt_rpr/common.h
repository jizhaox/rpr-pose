#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

// Reference:
// [1] Ji Zhao, Laurent Kneip, Yijia He, and Jiayi Ma.
//     Minimal Case Relative Pose Computation using Ray-Point-Ray Features.
//     IEEE Transactions on Pattern Analysis and Machine Intelligence,
//     42(5): 1176 - 1190, 2020.
// Author: Ji Zhao
// Email: zhaoji84@gmail.com

using namespace std;
using namespace Eigen;

#define NEAR_ZERO_THRESHOLD 1e-14

void format_covert_5pt(double*input, std::vector<Eigen::Vector3d>& P1, std::vector<Eigen::Vector3d>& P2)
{
    Eigen::Vector3d X1, X2;
    P1.clear();
    P1.clear();
    for (int k = 0; k < 5; k++)
    {
        for (int i = 0; i < 2; i++)
        {
            X1(i) = input[k*2+i];
        }
        X1(2) = 1.0;
        P1.push_back(X1);
    }
    for (int k = 5; k < 10; k++)
    {
        for (int i = 0; i < 2; i++)
        {
            X2(i) = input[k*2+i];
        }
        X2(2) = 1.0;
        P2.push_back(X2);
    }
    return;
}

void format_covert_2rpr90_1pt(double*input, std::vector<Eigen::Vector3d>& P1, std::vector<Eigen::Vector3d>& P2)
{
    Eigen::Vector3d X1, X2;
    P1.clear();
    P1.clear();

    // P1
    X1(0) = input[0];
    X1(1) = input[1];
    X1(2) = 1.0;
    P1.push_back(X1);

    X1(0) = input[6];
    X1(1) = input[7];
    P1.push_back(X1);

    X1(0) = input[12];
    X1(1) = input[13];
    P1.push_back(X1);
   
    // P2
    X2(0) = input[14];
    X2(1) = input[15];
    X2(2) = 1.0;
    P2.push_back(X2);

    X2(0) = input[20];
    X2(1) = input[21];
    P2.push_back(X2);

    X2(0) = input[26];
    X2(1) = input[27];
    P2.push_back(X2);
    return;
}

void format_covert_1rpr90_3pt(double*input, std::vector<Eigen::Vector3d>& P1, std::vector<Eigen::Vector3d>& P2)
{
    Eigen::Vector3d X1, X2;
    P1.clear();
    P1.clear();

    // P1
    X1(0) = input[0];
    X1(1) = input[1];
    X1(2) = 1.0;
    P1.push_back(X1);

    X1(0) = input[6];
    X1(1) = input[7];
    P1.push_back(X1);

    X1(0) = input[8];
    X1(1) = input[9];
    P1.push_back(X1);

    X1(0) = input[10];
    X1(1) = input[11];
    P1.push_back(X1);
   
    // P2
    X2(0) = input[12];
    X2(1) = input[13];
    X2(2) = 1.0;
    P2.push_back(X2);

    X2(0) = input[18];
    X2(1) = input[19];
    P2.push_back(X2);

    X2(0) = input[20];
    X2(1) = input[21];
    P2.push_back(X2);

    X2(0) = input[22];
    X2(1) = input[23];
    P2.push_back(X2);
    return;
}

void cayley2rotm(Eigen::Matrix<double,3,3>& rotm, Eigen::Matrix<double,3,1>& q)
{
    double qx = q(0);
    double qy = q(1);
    double qz = q(2);
    double qx2 = qx*qx;
    double qy2 = qy*qy;
    double qz2 = qz*qz;
    double s = 1.0/(1+qx2+qy2+qz2);
    rotm << 1+qx2-qy2-qz2, 2*qx*qy-2*qz, 2*qy+2*qx*qz, 
        2*qx*qy+2*qz, 1-qx2+qy2-qz2, 2*qy*qz-2*qx, 
        2*qx*qz-2*qy, 2*qx+2*qy*qz, 1-qx2-qy2+qz2;
    rotm = rotm*s;

    return;
}

void cayley2rotm(std::vector<Eigen::Matrix<double,3,3>>& rotm, std::vector<Eigen::Matrix<double,3,1>>& q_arr)
{
    rotm.clear();
    Eigen::Matrix<double,3,3> M;
    for (int i = 0; i < q_arr.size(); i++)
    {
        Eigen::Matrix<double,3,1> q;
        q = q_arr[i];
        double qx = q(0);
        double qy = q(1);
        double qz = q(2);
        double qx2 = qx*qx;
        double qy2 = qy*qy;
        double qz2 = qz*qz;
        double s = 1.0/(1+qx2+qy2+qz2);
        M << 1+qx2-qy2-qz2, 2*qx*qy-2*qz, 2*qy+2*qx*qz, 
            2*qx*qy+2*qz, 1-qx2+qy2-qz2, 2*qy*qz-2*qx, 
            2*qx*qz-2*qy, 2*qx+2*qy*qz, 1-qx2-qy2+qz2;
        M = M*s;
        rotm.push_back(M);
    }
    return;
}

void calculate_cross_dot_matrix(Eigen::MatrixXd& M, 
    std::vector<Eigen::Vector3d>& P1, std::vector<Eigen::Vector3d>& P2, 
    Eigen::Matrix3d& R)
{
    if (P1.size()<3 || P2.size()<3 || P1.size()!=P2.size())
    {
        std::cout << "size of input has an error!" << std::endl;
        return;
    }
    for(int i = 0; i < P1.size(); i++)
    {
        Eigen::Vector3d pt1 = P1[i];
        Eigen::Vector3d pt2 = P2[i];
        M.row(i) = pt2.cross(R*pt1);
    }
    return;
}

void calculate_translation(
    Eigen::MatrixXcd sols, std::vector<Eigen::Vector3d>& P1, std::vector<Eigen::Vector3d>& P2, 
    std::vector<Eigen::Matrix<double,3,3>>& rotm_arr, 
    std::vector<Eigen::Matrix<double,3,1>>& q_arr, std::vector<Eigen::Matrix<double,3,1>>& t_arr)
{
    q_arr.clear();
    t_arr.clear();
    if (sols.rows()!=3)
    {
        std::cout << "size of solution has an error!" << std::endl;
        return;
    }

    for (int j = 0; j < sols.cols(); j++)
    {
        std::complex<double> cx, cy, cz;
        double x, y, z;
        cx = sols(0, j);
        cy = sols(1, j);
        cz = sols(2, j);

        if (abs(cx.imag()) > NEAR_ZERO_THRESHOLD || abs(cy.imag()) > NEAR_ZERO_THRESHOLD || abs(cz.imag()) > NEAR_ZERO_THRESHOLD)
            continue;

        x = cx.real();
        y = cy.real();
        z = cz.real();

        Eigen::Vector3d q;
        q << x, y, z;
        q_arr.push_back(q);

        Eigen::Matrix<double,3,3> rotm;
        cayley2rotm(rotm, q);
        rotm_arr.push_back(rotm);

        Eigen::MatrixXd M(P1.size(), 3);
        calculate_cross_dot_matrix(M, P1, P2, rotm);

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXd V = svd.matrixV();
        Eigen::MatrixXd v_real = V.rightCols(1);

        Eigen::Matrix<double, 3, 1> t_vec;
        t_vec << v_real(0), v_real(1), v_real(2);
        t_arr.push_back(t_vec);
    }
    return;
}




