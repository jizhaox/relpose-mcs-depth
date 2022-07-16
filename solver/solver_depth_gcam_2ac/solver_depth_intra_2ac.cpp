#include "mex.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "common_stew.h"
#include "solver_depth_intra_2ac_core_48.h"

// Matlab mex command
// eigen_dir = '/usr/include/eigen3'; 
// mex(['-I"' eigen_dir '"'],'-O','solver_depth_intra_2ac.cpp')

using namespace std;
using namespace Eigen;

Eigen::MatrixXcd solve_equation_intra_cam_stew(double *input, double* zr, double* zi)
{
    const VectorXd data = Map<const VectorXd>(input, 1870); // 1870 = 83*20+35*6
//    std::cout << data << std::endl;
    Eigen::MatrixXcd sols = solver_depth_intra_2ac_core_48(data);
    for (Index i = 0; i < sols.size(); i++) {
        zr[i] = sols(i).real();
        zi[i] = sols(i).imag();
    }
    return sols;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    if (nrhs < 5)
        mexErrMsgTxt("at least 5 inputs are required!");
    if (nlhs < 1)
        mexErrMsgTxt("at least 1 output is required!");
    
    if (mxIsEmpty(prhs[0]) || mxIsEmpty(prhs[1]) || mxIsEmpty(prhs[2]) || mxIsEmpty(prhs[3]) || mxIsEmpty(prhs[4]))
        mexErrMsgTxt("input parameter should not be an empty array!");

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) 
        || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
        || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
        || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])
        || !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]))
    {
        mexErrMsgIdAndTxt("JPL:dp_intra_2ac:notDouble", "Input data must be type double.");
    }
    
    bool is_known_angle = false;
    AC_TYPE actype = AC_TYPE::INTRA_CAM_CONSTRAINT_FULL;

    double *Image_1;
    double *Image_2;
    double *affine_tran;
    double *extrinsic_R_camera;
    double *extrinsic_T_camera;

    Image_1 = (double *)mxGetData(prhs[0]);
    Image_2 = (double *)mxGetData(prhs[1]);
    affine_tran = (double *)mxGetData(prhs[2]);
    extrinsic_R_camera = (double *)mxGetData(prhs[3]);
    extrinsic_T_camera = (double *)mxGetData(prhs[4]);

    plhs[3] = mxCreateDoubleMatrix(3, 48, mxCOMPLEX);
    double* zr = mxGetPr(plhs[3]);
    double* zi = mxGetPi(plhs[3]);

    double *coeffs, *input;
    plhs[4] = mxCreateDoubleMatrix(26, 84, mxREAL);
    coeffs = mxGetPr(plhs[4]);

    input = new double[1870]; // 1870 = 83*20+35*6
    std::vector<std::vector<Eigen::Matrix<double,1,10>>> M;
    std::vector<Eigen::Matrix<double,6,1>> Line_i_all, Line_j_all;
    create_coeffs_stew(coeffs, input, M, Line_i_all, Line_j_all, Image_1, Image_2, affine_tran, extrinsic_R_camera, extrinsic_T_camera, actype, is_known_angle);
    Eigen::MatrixXcd sols = solve_equation_intra_cam_stew(input, zr, zi);

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
}
