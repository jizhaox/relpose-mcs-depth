#include "mex.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "common_stew.h"
#include "solver_depth_inter_2ac_ka_core_36.h"
#include "solver_depth_inter_2ac_ka_core_42.h"

// 
// Ji Zhao
// 09/26/2020

// Matlab mex command
// eigen_dir = '/usr/include/eigen3'; 
// mex(['-I"' eigen_dir '"'],'-O','solver_depth_inter_2ac_ka.cpp')

using namespace std;
using namespace Eigen;

Eigen::MatrixXcd solve_equation_inter_cam_stew_ka_36(double *input, double* zr, double* zi)
{
    const VectorXd data = Map<const VectorXd>(input, 781); // 781 = 84*8 + 4 + 35*3
//    std::cout << data << std::endl;
    Eigen::MatrixXcd sols = solver_depth_inter_2ac_ka_core_36(data);
    for (Index i = 0; i < sols.size(); i++) {
        zr[i] = sols(i).real();
        zi[i] = sols(i).imag();
    }
    return sols;
}

Eigen::MatrixXcd solve_equation_inter_cam_stew_ka_42(double *input, double* zr, double* zi)
{
    const VectorXd data = Map<const VectorXd>(input, 676); // 676 = 84*8 + 4
//    std::cout << data << std::endl;
    Eigen::MatrixXcd sols = solver_depth_inter_2ac_ka_core_42(data);
    for (Index i = 0; i < sols.size(); i++) {
        zr[i] = sols(i).real();
        zi[i] = sols(i).imag();
    }
    return sols;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    if (nrhs < 6)
        mexErrMsgTxt("at least 6 inputs are required!");
    if (nlhs < 1)
        mexErrMsgTxt("at least 1 output is required!");
    
    if (mxIsEmpty(prhs[0]) || mxIsEmpty(prhs[1]) || mxIsEmpty(prhs[2]) || mxIsEmpty(prhs[3]) || mxIsEmpty(prhs[4]) || mxIsEmpty(prhs[5]))
        mexErrMsgTxt("input parameter should not be an empty array!");
    
//    if (mxGetM(prhs[0])!=2 || mxGetN(prhs[0])!=2)
//        mexErrMsgTxt("the 1st parameter should be 2 by 2!");
    
//    if (mxGetM(prhs[1])!=2 || mxGetN(prhs[1])!=2)
//        mexErrMsgTxt("the 2nd parameter should be 2 by 2!");

 //   if (mxGetM(prhs[2])!=2 || mxGetN(prhs[2])!=2)
 //       mexErrMsgTxt("the 3rd parameter should be 2*2*2!");

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) 
        || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
        || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
        || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])
        || !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])
        || !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]))
    {
        mexErrMsgIdAndTxt("JPL:dp_inter_2ac_ka:notDouble", "Input data must be type double.");
    }
    
    bool is_known_angle = true;
    AC_TYPE actype = AC_TYPE::INTER_CAM_CONSTRAINT_FULL;
    if (nrhs >= 7)
    {
        double method = mxGetScalar(prhs[6]);
        if (method < -NEAR_ZERO_THRESHOLD || method > NEAR_ZERO_THRESHOLD)
        {
            actype = AC_TYPE::INTER_CAM_CONSTRAINT_PARTIAL;
        }
    }

    double *Image_1;
    double *Image_2;
    double *affine_tran;
    double *extrinsic_R_camera;
    double *extrinsic_T_camera;
    double theta;

    Image_1 = (double *)mxGetData(prhs[0]);
    Image_2 = (double *)mxGetData(prhs[1]);
    affine_tran = (double *)mxGetData(prhs[2]);
    extrinsic_R_camera = (double *)mxGetData(prhs[3]);
    extrinsic_T_camera = (double *)mxGetData(prhs[4]);
    theta = mxGetScalar(prhs[5]);
    double tmp = tan(theta*0.5);
    double sq_tan_half_theta = tmp*tmp;

    double *coeffs, *input;
    if (actype == AC_TYPE::INTER_CAM_CONSTRAINT_PARTIAL)
    {
        plhs[3] = mxCreateDoubleMatrix(3, 42, mxCOMPLEX);
        plhs[4] = mxCreateDoubleMatrix(20, 84, mxREAL);
        input = new double[676]; // 676 = 84*8 + 4
    }
    else
    {
        plhs[3] = mxCreateDoubleMatrix(3, 36, mxCOMPLEX);
        plhs[4] = mxCreateDoubleMatrix(26, 84, mxREAL);
        input = new double[781]; // 781 = 84*8 + 4 + 35*3
    }
    input[672] = 1.0;
    input[673] = 1.0;
    input[674] = 1.0;
    input[675] = -sq_tan_half_theta;

    double* zr = mxGetPr(plhs[3]);
    double* zi = mxGetPi(plhs[3]);
    coeffs = mxGetPr(plhs[4]);
    std::vector<std::vector<Eigen::Matrix<double,1,10>>> M;
    std::vector<Eigen::Matrix<double,6,1>> Line_i_all, Line_j_all;
    create_coeffs_stew(coeffs, input, M, Line_i_all, Line_j_all, Image_1, Image_2, affine_tran, extrinsic_R_camera, extrinsic_T_camera, actype, is_known_angle);
    Eigen::MatrixXcd sols;
    if (actype == AC_TYPE::INTER_CAM_CONSTRAINT_PARTIAL)
    {
        sols = solve_equation_inter_cam_stew_ka_42(input, zr, zi);
    }
    else
    {
        sols = solve_equation_inter_cam_stew_ka_36(input, zr, zi);
    }

    std::vector<Eigen::Matrix<double,3,1>> q_arr, t_arr;
    std::vector<Eigen::Matrix<double,3,3>> rotm;
    calculate_translation_stew(sols, M, Line_i_all, Line_j_all, q_arr, t_arr, rotm, is_known_angle);

    int n_real_sol = q_arr.size();
    double *q_real_sols, *t_real_sols, *rotm_real_sols;
    if (n_real_sol > 0)
    {
        plhs[0] = mxCreateDoubleMatrix(3, n_real_sol, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(3, n_real_sol, mxREAL);

        mwSize dims[3] = {3,3,n_real_sol};
        plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);

        q_real_sols = mxGetPr(plhs[0]);
        t_real_sols = mxGetPr(plhs[1]);
        rotm_real_sols = mxGetPr(plhs[2]);
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

    delete [] input;
    return;
}
