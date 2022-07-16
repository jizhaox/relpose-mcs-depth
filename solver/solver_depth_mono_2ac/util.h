#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define NEAR_ZERO_THRESHOLD 1e-14

void format_convert(
    double* input_Image_1, double* input_Image_2, double* input_affine_tran, 
    std::vector<Eigen::Vector2d>& Image1, std::vector<Eigen::Vector2d>& Image2, std::vector<Eigen::Matrix2d>& Ac)
{
    // data format convertion
    Eigen::Matrix2d A;
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                // Matlab uses column-first order
                A(i, j) = input_affine_tran[k*4+j*2+i];
            }
        }
        Ac.push_back(A);
    }

    Eigen::Vector2d X1, X2;
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < 2; i++)
        {
            X1(i) = input_Image_1[k*2+i];
            X2(i) = input_Image_2[k*2+i];
        }
        Image1.push_back(X1);
        Image2.push_back(X2);
    }
    
    return;
}

void set_value_column_major(double* array, int num_r, int num_c, int r, int c, double v)
{
    int idx = num_r*c + r;
    array[idx] = v;
    return;
}

void quad2rotm(Eigen::Matrix<double,3,3>& rotm, Eigen::Matrix<double,3,1>& q)
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

void quad2rotm(std::vector<Eigen::Matrix<double,3,3>>& rotm, std::vector<Eigen::Matrix<double,3,1>>& q_arr)
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

