#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

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

enum ROTATION_AXIS{
    X_AXIS,
    Y_AXIS,
    Z_AXIS
};

void format_covert_1rpr90_1pt(double*input, std::vector<Eigen::Vector3d>& P1, std::vector<Eigen::Vector3d>& P2)
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
   
    // P2
    X2(0) = input[8];
    X2(1) = input[9];
    X2(2) = 1.0;
    P2.push_back(X2);

    X2(0) = input[14];
    X2(1) = input[15];
    P2.push_back(X2);

    return;
}

void calculate_cross_dot_matrix(Eigen::MatrixXd& M, 
    std::vector<Eigen::Vector3d>& P1, std::vector<Eigen::Vector3d>& P2, 
    Eigen::Matrix3d& R)
{
    if (P1.size()!=2 || P2.size()!=2)
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


Eigen::VectorXcd solver_1rpr90_1pt(double *input, ROTATION_AXIS axis)
{
    double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12;
    double q0, q1, q2, q3, q4;
    Eigen::Vector3d P, D1, D2, Rslt;

    // RPR90 #1
    P << input[0], input[1], 1.0;
    D1 << input[2], input[3], 0;
    D2 << input[4], input[5], 0;

    Rslt = P.cross(D1);
    x1 = Rslt(0);
    x2 = Rslt(1);
    x3 = Rslt(2);

    Rslt = P.cross(D2);
    x4 = Rslt(0);
    x5 = Rslt(1);
    x6 = Rslt(2);

    // RPR90 #2
    P << input[8], input[9], 1.0;
    D1 << input[10], input[11], 0;
    D2 << input[12], input[13], 0;

    Rslt = P.cross(D1);
    x7 = Rslt(0);
    x8 = Rslt(1);
    x9 = Rslt(2);

    Rslt = P.cross(D2);
    x10 = Rslt(0);
    x11 = Rslt(1);
    x12 = Rslt(2);

    switch(axis)
    {
        case ROTATION_AXIS::X_AXIS:
        {
            q4 = x1*x4*x8*x11 + x1*x4*x9*x12 + x1*x5*x8*x10 + x1*x6*x9*x10 + x2*x4*x7*x11 + x2*x5*x7*x10 + x2*x5*x9*x12 - x2*x6*x9*x11 + x3*x4*x7*x12 - x3*x5*x8*x12 + x3*x6*x7*x10 + x3*x6*x8*x11;
            q3 = -2*x1*x5*x9*x10 + 2*x1*x6*x8*x10 - 2*x2*x4*x7*x12 + 2*x2*x5*x8*x12 + 2*x2*x5*x9*x11 - 2*x2*x6*x8*x11 + 2*x2*x6*x9*x12 + 2*x3*x4*x7*x11 - 2*x3*x5*x8*x11 + 2*x3*x5*x9*x12 - 2*x3*x6*x8*x12 - 2*x3*x6*x9*x11;
            q2 = 2*x1*x4*x8*x11 + 2*x1*x4*x9*x12 + 2*x2*x5*x7*x10 + 4*x2*x5*x8*x11 - 2*x2*x5*x9*x12 + 4*x2*x6*x8*x12 + 2*x2*x6*x9*x11 + 2*x3*x5*x8*x12 + 4*x3*x5*x9*x11 + 2*x3*x6*x7*x10 - 2*x3*x6*x8*x11 + 4*x3*x6*x9*x12;
            q1 = -2*x1*x5*x9*x10 + 2*x1*x6*x8*x10 - 2*x2*x4*x7*x12 - 2*x2*x5*x8*x12 - 2*x2*x5*x9*x11 + 2*x2*x6*x8*x11 - 2*x2*x6*x9*x12 + 2*x3*x4*x7*x11 + 2*x3*x5*x8*x11 - 2*x3*x5*x9*x12 + 2*x3*x6*x8*x12 + 2*x3*x6*x9*x11;
            q0 = x1*x4*x8*x11 + x1*x4*x9*x12 - x1*x5*x8*x10 - x1*x6*x9*x10 - x2*x4*x7*x11 + x2*x5*x7*x10 + x2*x5*x9*x12 - x2*x6*x9*x11 - x3*x4*x7*x12 - x3*x5*x8*x12 + x3*x6*x7*x10 + x3*x6*x8*x11;
            break;
        }
        case ROTATION_AXIS::Y_AXIS:
        {
            q4 = x1*x4*x8*x11 + x1*x4*x9*x12 + x1*x5*x8*x10 - x1*x6*x9*x10 + x2*x4*x7*x11 + x2*x5*x7*x10 + x2*x5*x9*x12 + x2*x6*x9*x11 - x3*x4*x7*x12 + x3*x5*x8*x12 + x3*x6*x7*x10 + x3*x6*x8*x11;
            q3 = -2*x1*x4*x7*x12 - 2*x1*x4*x9*x10 + 2*x1*x5*x8*x12 + 2*x1*x6*x7*x10 - 2*x1*x6*x9*x12 + 2*x2*x4*x9*x11 - 2*x2*x6*x7*x11 + 2*x3*x4*x7*x10 - 2*x3*x4*x9*x12 - 2*x3*x5*x8*x10 + 2*x3*x6*x7*x12 + 2*x3*x6*x9*x10;
            q2 = 4*x1*x4*x7*x10 + 2*x1*x4*x8*x11 - 2*x1*x4*x9*x12 + 4*x1*x6*x7*x12 + 2*x1*x6*x9*x10 + 2*x2*x5*x7*x10 + 2*x2*x5*x9*x12 + 2*x3*x4*x7*x12 + 4*x3*x4*x9*x10 - 2*x3*x6*x7*x10 + 2*x3*x6*x8*x11 + 4*x3*x6*x9*x12;
            q1 = 2*x1*x4*x7*x12 + 2*x1*x4*x9*x10 + 2*x1*x5*x8*x12 - 2*x1*x6*x7*x10 + 2*x1*x6*x9*x12 + 2*x2*x4*x9*x11 - 2*x2*x6*x7*x11 - 2*x3*x4*x7*x10 + 2*x3*x4*x9*x12 - 2*x3*x5*x8*x10 - 2*x3*x6*x7*x12 - 2*x3*x6*x9*x10;
            q0 = x1*x4*x8*x11 + x1*x4*x9*x12 - x1*x5*x8*x10 - x1*x6*x9*x10 - x2*x4*x7*x11 + x2*x5*x7*x10 + x2*x5*x9*x12 - x2*x6*x9*x11 - x3*x4*x7*x12 - x3*x5*x8*x12 + x3*x6*x7*x10 + x3*x6*x8*x11;
            break;

        }
        case ROTATION_AXIS::Z_AXIS:
        {
            q4 = x1*x4*x8*x11 + x1*x4*x9*x12 - x1*x5*x8*x10 + x1*x6*x9*x10 - x2*x4*x7*x11 + x2*x5*x7*x10 + x2*x5*x9*x12 + x2*x6*x9*x11 + x3*x4*x7*x12 + x3*x5*x8*x12 + x3*x6*x7*x10 + x3*x6*x8*x11;
            q3 = 2*x1*x4*x7*x11 + 2*x1*x4*x8*x10 - 2*x1*x5*x7*x10 + 2*x1*x5*x8*x11 - 2*x1*x6*x9*x11 - 2*x2*x4*x7*x10 + 2*x2*x4*x8*x11 - 2*x2*x5*x7*x11 - 2*x2*x5*x8*x10 + 2*x2*x6*x9*x10 - 2*x3*x4*x8*x12 + 2*x3*x5*x7*x12;
            q2 = 4*x1*x4*x7*x10 - 2*x1*x4*x8*x11 + 2*x1*x4*x9*x12 + 4*x1*x5*x7*x11 + 2*x1*x5*x8*x10 + 2*x2*x4*x7*x11 + 4*x2*x4*x8*x10 - 2*x2*x5*x7*x10 + 4*x2*x5*x8*x11 + 2*x2*x5*x9*x12 + 2*x3*x6*x7*x10 + 2*x3*x6*x8*x11;
            q1 = -2*x1*x4*x7*x11 - 2*x1*x4*x8*x10 + 2*x1*x5*x7*x10 - 2*x1*x5*x8*x11 - 2*x1*x6*x9*x11 + 2*x2*x4*x7*x10 - 2*x2*x4*x8*x11 + 2*x2*x5*x7*x11 + 2*x2*x5*x8*x10 + 2*x2*x6*x9*x10 - 2*x3*x4*x8*x12 + 2*x3*x5*x7*x12;
            q0 = x1*x4*x8*x11 + x1*x4*x9*x12 - x1*x5*x8*x10 - x1*x6*x9*x10 - x2*x4*x7*x11 + x2*x5*x7*x10 + x2*x5*x9*x12 - x2*x6*x9*x11 - x3*x4*x7*x12 - x3*x5*x8*x12 + x3*x6*x7*x10 + x3*x6*x8*x11;
            break;
        }
        default:
        {
            std::cout << "axis type 2 has an error!" << std::endl;
        }
    }
    Eigen::Matrix<double, 4, 4> A;
    A(0, 0) = -q3/q4;
    A(0, 1) = -q2/q4;
    A(0, 2) = -q1/q4;
    A(0, 3) = -q0/q4;
    A(1, 0) = 1.0;
    A(2, 1) = 1.0;
    A(3, 2) = 1.0;

    Eigen::EigenSolver<Eigen::Matrix<double,4,4>> es(A);
    Eigen::VectorXcd ev = es.eigenvalues();

    return ev;
}

void calculate_translation_common_direction(
    Eigen::VectorXcd sols, std::vector<Eigen::Vector3d>& P1, std::vector<Eigen::Vector3d>& P2, 
    std::vector<Eigen::Matrix<double,3,3>>& rotm_arr, 
    std::vector<double>& q_arr, std::vector<Eigen::Matrix<double,3,1>>& t_arr,
    ROTATION_AXIS axis)
{
    q_arr.clear();
    t_arr.clear();
//    if (sols.rows()!=3)
//    {
//        std::cout << "size of solution has an error!" << std::endl;
//        return;
//    }

    for (int j = 0; j < 4; j++)
    {
        if (abs(sols(j).imag()) > NEAR_ZERO_THRESHOLD)
            continue;
        double q = sols(j).real();
        q_arr.push_back(q);

        Eigen::Matrix<double,3,3> rotm;
        rotm.setZero();
        double cos_theta = (1-q*q)/(1+q*q);
        double sin_theta = 2*q/(1+q*q);

        switch(axis)
        {
            case ROTATION_AXIS::X_AXIS:
            {
                rotm(0, 0) = 1.0;
                rotm(1, 1) = cos_theta;
                rotm(1, 2) = -sin_theta;
                rotm(2, 1) = sin_theta;
                rotm(2, 2) = cos_theta;
                break;
            }
            case ROTATION_AXIS::Y_AXIS:
            {
                rotm(0, 0) = cos_theta;
                rotm(0, 2) = sin_theta;
                rotm(1, 1) = 1.0;
                rotm(2, 0) = -sin_theta;
                rotm(2, 2) = cos_theta;
                break;
            }
            case ROTATION_AXIS::Z_AXIS:
            {
                rotm(0, 0) = cos_theta;
                rotm(0, 1) = -sin_theta;
                rotm(1, 0) = sin_theta;
                rotm(1, 1) = cos_theta;
                rotm(2, 2) = 1.0;
                break;
            }
            default:
            {
                std::cout << "axis type 1 has an error!" << std::endl;
            }
        }        
        rotm_arr.push_back(rotm);

        Eigen::MatrixXd M(P1.size(), 3);
        calculate_cross_dot_matrix(M, P1, P2, rotm);
        
        Eigen::Vector3d V1 = M.row(0);
        Eigen::Vector3d V2 = M.row(1);

        Eigen::Vector3d t_vec = V1.cross(V2);
        t_vec.normalize();
        t_arr.push_back(t_vec);
    }
    return;
}
