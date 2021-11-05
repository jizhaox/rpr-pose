#include "mex.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "common.h"
#include "solver_1rpr90_3pt_cayley.h"

// Reference:
// [1] Ji Zhao, Laurent Kneip, Yijia He, and Jiayi Ma.
//     Minimal Case Relative Pose Computation using Ray-Point-Ray Features.
//     IEEE Transactions on Pattern Analysis and Machine Intelligence,
//     42(5): 1176 - 1190, 2020.
// Author: Ji Zhao
// Email: zhaoji84@gmail.com

// Matlab mex command
// eigen_dir = '/usr/include/eigen3'; 
// mex(['-I"' eigen_dir '"'],'-O','solver_1rpr90_3pt.cpp')

using namespace std;
using namespace Eigen;

Eigen::MatrixXcd solve_equation_1rpr90_3pt_cayley(double *input)
{
    const VectorXd data = Map<const VectorXd>(input, 24);
    Eigen::MatrixXcd sols = solver_1rpr90_3pt_cayley(data);
    return sols;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    if (nrhs < 1)
        mexErrMsgTxt("1 input is required!");
    if (nlhs < 1)
        mexErrMsgTxt("at least 1 output is required!");
    if (nlhs > 3)
        mexErrMsgTxt("at most 3 outputs!");
    
    if (mxIsEmpty(prhs[0]))
        mexErrMsgTxt("input parameter should not be an empty array!");
    
    if (mxGetM(prhs[0])*mxGetN(prhs[0])!=24)
        mexErrMsgTxt("the 1st parameter should be 24 by 1!");

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    {
        mexErrMsgIdAndTxt("2rpr90_1pt:notDouble", "Input data must be type double.");
    }
    
    double *input = (double *)mxGetData(prhs[0]);
    
    Eigen::MatrixXcd sols;
    sols = solve_equation_1rpr90_3pt_cayley(input);    

    std::vector<Eigen::Matrix<double,3,1>> q_arr, t_arr;
    std::vector<Eigen::Matrix<double,3,3>> rotm;
    std::vector<Eigen::Vector3d> P1, P2;
    format_covert_1rpr90_3pt(input, P1, P2);
    calculate_translation(sols, P1, P2, rotm, q_arr, t_arr);

    int n_real_sol = q_arr.size();
    double *q_real_sols, *t_real_sols, *rotm_real_sols;
    if (n_real_sol > 0)
    {
        mwSize dims[3] = {3,3,n_real_sol};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(3, n_real_sol, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(3, n_real_sol, mxREAL);

        rotm_real_sols = mxGetPr(plhs[0]);
        t_real_sols = mxGetPr(plhs[1]);
        q_real_sols = mxGetPr(plhs[2]);
        
        for (int i = 0; i < n_real_sol; i++)
        {
            Eigen::Matrix<double,3,1> q = q_arr[i];
            Eigen::Matrix<double,3,1> t = t_arr[i];
            Eigen::Matrix<double,3,3> r = rotm[i];
            q_real_sols[i*3] = q(0);
            q_real_sols[i*3+1] = q(1);
            q_real_sols[i*3+2] = q(2);
            t_real_sols[i*3] = t(0);
            t_real_sols[i*3+1] = t(1);
            t_real_sols[i*3+2] = t(2);

            rotm_real_sols[i*9] = r(0,0);
            rotm_real_sols[i*9+1] = r(1,0);
            rotm_real_sols[i*9+2] = r(2,0);
            rotm_real_sols[i*9+3] = r(0,1);
            rotm_real_sols[i*9+4] = r(1,1);
            rotm_real_sols[i*9+5] = r(2,1);
            rotm_real_sols[i*9+6] = r(0,2);
            rotm_real_sols[i*9+7] = r(1,2);
            rotm_real_sols[i*9+8] = r(2,2);
        }
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }

    return;
}
