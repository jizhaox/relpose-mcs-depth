#include "mex.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "common_stew_mono.h"
#include "solver_depth_mono_2ac_core.h"

// Matlab mex command
// eigen_dir = '/usr/include/eigen3'; 
// mex(['-I"' eigen_dir '"'],'-O','solver_depth_mono_2ac.cpp')

using namespace std;
using namespace Eigen;


Eigen::MatrixXcd solver_depth_mono_2ac(double *input, double* zr, double* zi)
{
    const VectorXd data = Map<const VectorXd>(input, 420); // 35*12
    Eigen::MatrixXcd sols = solver_depth_mono_2ac_core(data);
    for (Index i = 0; i < sols.size(); i++) {
        zr[i] = sols(i).real();
        zi[i] = sols(i).imag();
    }
    return sols;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    if (nrhs != 3)
        mexErrMsgTxt("3 inputs are required!");
    if (nlhs < 1)
        mexErrMsgTxt("at least 1 output is required!");
    
    if (mxIsEmpty(prhs[0]) || mxIsEmpty(prhs[1]) || mxIsEmpty(prhs[2]))
        mexErrMsgTxt("input parameter should not be an empty array!");
    
    if (mxGetM(prhs[0])!=2 || mxGetN(prhs[0])!=2)
        mexErrMsgTxt("the 1st parameter should be 2 by 2!");
    
    if (mxGetM(prhs[1])!=2 || mxGetN(prhs[1])!=2)
        mexErrMsgTxt("the 2nd parameter should be 2 by 2!");

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) 
        || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
        || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
    {
        mexErrMsgIdAndTxt("JPL:dp_mono_2ac:notDouble", "Input data must be type double.");
    }

    double *Image_1;
    double *Image_2;
    double *affine_tran;

    Image_1 = (double *)mxGetData(prhs[0]);
    Image_2 = (double *)mxGetData(prhs[1]);
    affine_tran = (double *)mxGetData(prhs[2]);

    plhs[3] = mxCreateDoubleMatrix(3, 20, mxCOMPLEX);
    double* zr = mxGetPr(plhs[3]);
    double* zi = mxGetPi(plhs[3]);

    double *coeffs, *input;
    int r = 12;
    int c = 35;
    plhs[4] = mxCreateDoubleMatrix(r, c, mxREAL);
    coeffs = mxGetPr(plhs[4]);

    input = new double[420]; // 35*12

    std::vector<std::vector<Eigen::Matrix<double,1,10>>> M;
    create_coeffs_stew_single(coeffs, input, M, Image_1, Image_2, affine_tran);
    Eigen::MatrixXcd sols = solver_depth_mono_2ac(input, zr, zi);

    std::vector<Eigen::Matrix<double,3,1>> q_arr, t_arr;
    std::vector<Eigen::Matrix<double,3,3>> rotm;
    bool is_known_angle = false;
    calculate_translation_stew_single(sols, M, Image_1, Image_2, q_arr, t_arr, rotm, is_known_angle);

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
