#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define NEAR_ZERO_THRESHOLD 1e-14

enum AC_TYPE{
    INTER_CAM_CONSTRAINT_FULL,
    INTRA_CAM_CONSTRAINT_FULL,
    INTER_CAM_CONSTRAINT_PARTIAL,
    INTRA_CAM_CONSTRAINT_PARTIAL
};

void format_convert(
    double* input_Image_1, double* input_Image_2, double* input_affine_tran, 
    double* extrinsic_R_camera, double* extrinsic_T_camera, 
    std::vector<Eigen::Vector3d>& Image1, std::vector<Eigen::Vector3d>& Image2, std::vector<Eigen::Matrix3d>& Ac,
    std::vector<Eigen::Matrix3d>& R_camera, std::vector<Eigen::Vector3d>& T_camera)
{
    // data format convertion
    Eigen::Matrix3d R;
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                // Matlab uses column-first order
                R(i, j) = extrinsic_R_camera[k*9+j*3+i];
            }
        }
        R_camera.push_back(R);
    }
    Eigen::Matrix3d A;
    for (int k = 0; k < 2; k++)
    {
        A.setZero();
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

    Eigen::Vector3d T;
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            T(i) = extrinsic_T_camera[k*3+i];
        }
        T_camera.push_back(T);
    }
    Eigen::Vector3d X1, X2;
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            X1(i) = input_Image_1[k*3+i];
            X2(i) = input_Image_2[k*3+i];
        }
        Image1.push_back(X1);
        Image2.push_back(X2);
    }

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

void var3_order2_two_multiplication(Eigen::Matrix<double,1,10>& a_arr, Eigen::Matrix<double,1,10>& b_arr, Eigen::Matrix<double,1,35>& c_arr)
{
    // x^2, x*y, x*z, x, y^2, y*z, y, z^2, z, 1
    double a_x2 = a_arr(0);
    double a_xy = a_arr(1);
    double a_xz = a_arr(2);
    double a_x  = a_arr(3);
    double a_y2 = a_arr(4);
    double a_yz = a_arr(5);
    double a_y  = a_arr(6);
    double a_z2 = a_arr(7);
    double a_z  = a_arr(8);
    double a_1  = a_arr(9);

    double b_x2 = b_arr(0);
    double b_xy = b_arr(1);
    double b_xz = b_arr(2);
    double b_x  = b_arr(3);
    double b_y2 = b_arr(4);
    double b_yz = b_arr(5);
    double b_y  = b_arr(6);
    double b_z2 = b_arr(7);
    double b_z  = b_arr(8);
    double b_1  = b_arr(9);

    // results
    // x^4, x^3*y, x^3*z, x^3, x^2*y^2, x^2*y*z, x^2*y, x^2*z^2, x^2*z, x^2, 
    // x*y^3, x*y^2*z, x*y^2, x*y*z^2, x*y*z, x*y, x*z^3, x*z^2, x*z, x, 
    // y^4, y^3*z, y^3, y^2*z^2, y^2*z, y^2, y*z^3, y*z^2, y*z, y, 
    // z^4, z^3, z^2, z, 1
    c_arr(0) = a_x2*b_x2;
    c_arr(1) = a_x2*b_xy + a_xy*b_x2;
    c_arr(2) = a_x2*b_xz + a_xz*b_x2;
    c_arr(3) = a_x*b_x2 + a_x2*b_x;
    c_arr(4) = a_xy*b_xy + a_x2*b_y2 + a_y2*b_x2;
    c_arr(5) = a_xy*b_xz + a_xz*b_xy + a_x2*b_yz + a_yz*b_x2;
    c_arr(6) = a_x*b_xy + a_xy*b_x + a_x2*b_y + a_y*b_x2;
    c_arr(7) = a_xz*b_xz + a_x2*b_z2 + a_z2*b_x2;
    c_arr(8) = a_x*b_xz + a_xz*b_x + a_x2*b_z + a_z*b_x2;
    c_arr(9) = a_1*b_x2 + a_x2*b_1 + a_x*b_x;
    c_arr(10) = a_xy*b_y2 + a_y2*b_xy;
    c_arr(11) = a_xz*b_y2 + a_y2*b_xz + a_xy*b_yz + a_yz*b_xy;
    c_arr(12) = a_x*b_y2 + a_y2*b_x + a_xy*b_y + a_y*b_xy;
    c_arr(13) = a_xz*b_yz + a_yz*b_xz + a_xy*b_z2 + a_z2*b_xy;
    c_arr(14) = a_x*b_yz + a_xz*b_y + a_y*b_xz + a_yz*b_x + a_xy*b_z + a_z*b_xy;
    c_arr(15) = a_1*b_xy + a_xy*b_1 + a_x*b_y + a_y*b_x;
    c_arr(16) = a_xz*b_z2 + a_z2*b_xz;
    c_arr(17) = a_x*b_z2 + a_z2*b_x + a_xz*b_z + a_z*b_xz;
    c_arr(18) = a_1*b_xz + a_xz*b_1 + a_x*b_z + a_z*b_x;
    c_arr(19) = a_1*b_x + a_x*b_1;
    c_arr(20) = a_y2*b_y2;
    c_arr(21) = a_y2*b_yz + a_yz*b_y2;
    c_arr(22) = a_y*b_y2 + a_y2*b_y;
    c_arr(23) = a_yz*b_yz + a_y2*b_z2 + a_z2*b_y2;
    c_arr(24) = a_y*b_yz + a_yz*b_y + a_y2*b_z + a_z*b_y2;
    c_arr(25) = a_1*b_y2 + a_y2*b_1 + a_y*b_y;
    c_arr(26) = a_yz*b_z2 + a_z2*b_yz;
    c_arr(27) = a_y*b_z2 + a_z2*b_y + a_yz*b_z + a_z*b_yz;
    c_arr(28) = a_1*b_yz + a_yz*b_1 + a_y*b_z + a_z*b_y;
    c_arr(29) = a_1*b_y + a_y*b_1;
    c_arr(30) = a_z2*b_z2;
    c_arr(31) = a_z*b_z2 + a_z2*b_z;
    c_arr(32) = a_1*b_z2 + a_z2*b_1 + a_z*b_z;
    c_arr(33) = a_1*b_z + a_z*b_1;
    c_arr(34) = a_1*b_1;

    return;
}

void var3_order4_multiplication(Eigen::Matrix<double,1,35>& a_arr, Eigen::Matrix<double,1,35>& b_arr, Eigen::Matrix<double,1,165>& c_arr)
{
    // x^4, x^3*y, x^3*z, x^3, x^2*y^2, x^2*y*z, x^2*y, x^2*z^2, x^2*z, x^2, 
    // x*y^3, x*y^2*z, x*y^2, x*y*z^2, x*y*z, x*y, x*z^3, x*z^2, x*z, x, 
    // y^4, y^3*z, y^3, y^2*z^2, y^2*z, y^2, y*z^3, y*z^2, y*z, y, 
    // z^4, z^3, z^2, z, 1
    double a_x4   = a_arr(0);
    double a_x3y  = a_arr(1);
    double a_x3z  = a_arr(2);
    double a_x3   = a_arr(3);
    double a_x2y2 = a_arr(4);
    double a_x2yz = a_arr(5);
    double a_x2y  = a_arr(6);
    double a_x2z2 = a_arr(7);
    double a_x2z  = a_arr(8);
    double a_x2   = a_arr(9);
    double a_xy3  = a_arr(10);
    double a_xy2z = a_arr(11);
    double a_xy2  = a_arr(12);
    double a_xyz2 = a_arr(13);
    double a_xyz  = a_arr(14);
    double a_xy   = a_arr(15);
    double a_xz3  = a_arr(16);
    double a_xz2  = a_arr(17);
    double a_xz   = a_arr(18);
    double a_x    = a_arr(19);
    double a_y4   = a_arr(20);
    double a_y3z  = a_arr(21);
    double a_y3   = a_arr(22);
    double a_y2z2 = a_arr(23);
    double a_y2z  = a_arr(24);
    double a_y2   = a_arr(25);
    double a_yz3  = a_arr(26);
    double a_yz2  = a_arr(27);
    double a_yz   = a_arr(28);
    double a_y    = a_arr(29);
    double a_z4   = a_arr(30);
    double a_z3   = a_arr(31);
    double a_z2   = a_arr(32);
    double a_z    = a_arr(33);
    double a_1    = a_arr(34);

    double b_x4   = b_arr(0);
    double b_x3y  = b_arr(1);
    double b_x3z  = b_arr(2);
    double b_x3   = b_arr(3);
    double b_x2y2 = b_arr(4);
    double b_x2yz = b_arr(5);
    double b_x2y  = b_arr(6);
    double b_x2z2 = b_arr(7);
    double b_x2z  = b_arr(8);
    double b_x2   = b_arr(9);
    double b_xy3  = b_arr(10);
    double b_xy2z = b_arr(11);
    double b_xy2  = b_arr(12);
    double b_xyz2 = b_arr(13);
    double b_xyz  = b_arr(14);
    double b_xy   = b_arr(15);
    double b_xz3  = b_arr(16);
    double b_xz2  = b_arr(17);
    double b_xz   = b_arr(18);
    double b_x    = b_arr(19);
    double b_y4   = b_arr(20);
    double b_y3z  = b_arr(21);
    double b_y3   = b_arr(22);
    double b_y2z2 = b_arr(23);
    double b_y2z  = b_arr(24);
    double b_y2   = b_arr(25);
    double b_yz3  = b_arr(26);
    double b_yz2  = b_arr(27);
    double b_yz   = b_arr(28);
    double b_y    = b_arr(29);
    double b_z4   = b_arr(30);
    double b_z3   = b_arr(31);
    double b_z2   = b_arr(32);
    double b_z    = b_arr(33);
    double b_1    = b_arr(34);

    c_arr(0) = a_x4*b_x4;
    c_arr(1) = a_x4*b_x3y + a_x3y*b_x4;
    c_arr(2) = a_x4*b_x3z + a_x3z*b_x4;
    c_arr(3) = a_x3*b_x4 + a_x4*b_x3;
    c_arr(4) = a_x4*b_x2y2 + a_x2y2*b_x4 + a_x3y*b_x3y;
    c_arr(5) = a_x4*b_x2yz + a_x3y*b_x3z + a_x3z*b_x3y + a_x2yz*b_x4;
    c_arr(6) = a_x3*b_x3y + a_x4*b_x2y + a_x2y*b_x4 + a_x3y*b_x3;
    c_arr(7) = a_x4*b_x2z2 + a_x2z2*b_x4 + a_x3z*b_x3z;
    c_arr(8) = a_x3*b_x3z + a_x4*b_x2z + a_x2z*b_x4 + a_x3z*b_x3;
    c_arr(9) = a_x2*b_x4 + a_x3*b_x3 + a_x4*b_x2;
    c_arr(10) = a_x4*b_xy3 + a_xy3*b_x4 + a_x3y*b_x2y2 + a_x2y2*b_x3y;
    c_arr(11) = a_x4*b_xy2z + a_xy2z*b_x4 + a_x3z*b_x2y2 + a_x2y2*b_x3z + a_x3y*b_x2yz + a_x2yz*b_x3y;
    c_arr(12) = a_x4*b_xy2 + a_xy2*b_x4 + a_x3*b_x2y2 + a_x2y2*b_x3 + a_x2y*b_x3y + a_x3y*b_x2y;
    c_arr(13) = a_x4*b_xyz2 + a_xyz2*b_x4 + a_x3y*b_x2z2 + a_x2z2*b_x3y + a_x3z*b_x2yz + a_x2yz*b_x3z;
    c_arr(14) = a_x4*b_xyz + a_xyz*b_x4 + a_x3*b_x2yz + a_x2y*b_x3z + a_x2z*b_x3y + a_x3y*b_x2z + a_x3z*b_x2y + a_x2yz*b_x3;
    c_arr(15) = a_x4*b_xy + a_xy*b_x4 + a_x2*b_x3y + a_x3*b_x2y + a_x2y*b_x3 + a_x3y*b_x2;
    c_arr(16) = a_x4*b_xz3 + a_xz3*b_x4 + a_x3z*b_x2z2 + a_x2z2*b_x3z;
    c_arr(17) = a_x4*b_xz2 + a_xz2*b_x4 + a_x3*b_x2z2 + a_x2z2*b_x3 + a_x2z*b_x3z + a_x3z*b_x2z;
    c_arr(18) = a_x4*b_xz + a_xz*b_x4 + a_x2*b_x3z + a_x3*b_x2z + a_x2z*b_x3 + a_x3z*b_x2;
    c_arr(19) = a_x*b_x4 + a_x4*b_x + a_x2*b_x3 + a_x3*b_x2;
    c_arr(20) = a_x3y*b_xy3 + a_xy3*b_x3y + a_x4*b_y4 + a_y4*b_x4 + a_x2y2*b_x2y2;
    c_arr(21) = a_xy3*b_x3z + a_x3z*b_xy3 + a_x3y*b_xy2z + a_xy2z*b_x3y + a_x4*b_y3z + a_y3z*b_x4 + a_x2y2*b_x2yz + a_x2yz*b_x2y2;
    c_arr(22) = a_x3*b_xy3 + a_xy3*b_x3 + a_xy2*b_x3y + a_x3y*b_xy2 + a_x4*b_y3 + a_y3*b_x4 + a_x2y*b_x2y2 + a_x2y2*b_x2y;
    c_arr(23) = a_x2y2*b_x2z2 + a_x2z2*b_x2y2 + a_x3y*b_xyz2 + a_xyz2*b_x3y + a_x3z*b_xy2z + a_xy2z*b_x3z + a_x4*b_y2z2 + a_y2z2*b_x4 + a_x2yz*b_x2yz;
    c_arr(24) = a_x3*b_xy2z + a_xy2*b_x3z + a_x3z*b_xy2 + a_xy2z*b_x3 + a_x2z*b_x2y2 + a_x2y2*b_x2z + a_x3y*b_xyz + a_xyz*b_x3y + a_x2y*b_x2yz + a_x2yz*b_x2y + a_x4*b_y2z + a_y2z*b_x4;
    c_arr(25) = a_x3*b_xy2 + a_xy2*b_x3 + a_x2*b_x2y2 + a_x2y2*b_x2 + a_xy*b_x3y + a_x3y*b_xy + a_x2y*b_x2y + a_x4*b_y2 + a_y2*b_x4;
    c_arr(26) = a_x3y*b_xz3 + a_xz3*b_x3y + a_x3z*b_xyz2 + a_xyz2*b_x3z + a_x4*b_yz3 + a_yz3*b_x4 + a_x2z2*b_x2yz + a_x2yz*b_x2z2;
    c_arr(27) = a_x3*b_xyz2 + a_x3y*b_xz2 + a_xz2*b_x3y + a_xyz2*b_x3 + a_x2y*b_x2z2 + a_x2z2*b_x2y + a_x3z*b_xyz + a_xyz*b_x3z + a_x2z*b_x2yz + a_x2yz*b_x2z + a_x4*b_yz2 + a_yz2*b_x4;
    c_arr(28) = a_x3*b_xyz + a_xy*b_x3z + a_xz*b_x3y + a_x3y*b_xz + a_x3z*b_xy + a_xyz*b_x3 + a_x2*b_x2yz + a_x2y*b_x2z + a_x2z*b_x2y + a_x2yz*b_x2 + a_x4*b_yz + a_yz*b_x4;
    c_arr(29) = a_x*b_x3y + a_x3*b_xy + a_xy*b_x3 + a_x3y*b_x + a_x2*b_x2y + a_x2y*b_x2 + a_x4*b_y + a_y*b_x4;
    c_arr(30) = a_x3z*b_xz3 + a_xz3*b_x3z + a_x2z2*b_x2z2 + a_x4*b_z4 + a_z4*b_x4;
    c_arr(31) = a_x3*b_xz3 + a_xz3*b_x3 + a_xz2*b_x3z + a_x3z*b_xz2 + a_x2z*b_x2z2 + a_x2z2*b_x2z + a_x4*b_z3 + a_z3*b_x4;
    c_arr(32) = a_x3*b_xz2 + a_xz2*b_x3 + a_x2*b_x2z2 + a_x2z2*b_x2 + a_xz*b_x3z + a_x3z*b_xz + a_x2z*b_x2z + a_x4*b_z2 + a_z2*b_x4;
    c_arr(33) = a_x*b_x3z + a_x3*b_xz + a_xz*b_x3 + a_x3z*b_x + a_x2*b_x2z + a_x2z*b_x2 + a_x4*b_z + a_z*b_x4;
    c_arr(34) = a_1*b_x4 + a_x4*b_1 + a_x*b_x3 + a_x3*b_x + a_x2*b_x2;
    c_arr(35) = a_xy3*b_x2y2 + a_x2y2*b_xy3 + a_x3y*b_y4 + a_y4*b_x3y;
    c_arr(36) = a_xy3*b_x2yz + a_x2yz*b_xy3 + a_x3z*b_y4 + a_y4*b_x3z + a_x2y2*b_xy2z + a_xy2z*b_x2y2 + a_x3y*b_y3z + a_y3z*b_x3y;
    c_arr(37) = a_x2y*b_xy3 + a_xy3*b_x2y + a_x3*b_y4 + a_y4*b_x3 + a_xy2*b_x2y2 + a_x2y2*b_xy2 + a_x3y*b_y3 + a_y3*b_x3y;
    c_arr(38) = a_xy3*b_x2z2 + a_x2z2*b_xy3 + a_x2y2*b_xyz2 + a_xyz2*b_x2y2 + a_x2yz*b_xy2z + a_xy2z*b_x2yz + a_x3z*b_y3z + a_y3z*b_x3z + a_x3y*b_y2z2 + a_y2z2*b_x3y;
    c_arr(39) = a_x2z*b_xy3 + a_xy3*b_x2z + a_x2y*b_xy2z + a_xy2*b_x2yz + a_x2yz*b_xy2 + a_xy2z*b_x2y + a_x3*b_y3z + a_x3z*b_y3 + a_y3*b_x3z + a_y3z*b_x3 + a_x2y2*b_xyz + a_xyz*b_x2y2 + a_x3y*b_y2z + a_y2z*b_x3y;
    c_arr(40) = a_x2*b_xy3 + a_xy3*b_x2 + a_x2y*b_xy2 + a_xy2*b_x2y + a_x3*b_y3 + a_y3*b_x3 + a_xy*b_x2y2 + a_x2y2*b_xy + a_x3y*b_y2 + a_y2*b_x3y;
    c_arr(41) = a_xz3*b_x2y2 + a_x2y2*b_xz3 + a_x2z2*b_xy2z + a_xy2z*b_x2z2 + a_x2yz*b_xyz2 + a_xyz2*b_x2yz + a_x3y*b_yz3 + a_yz3*b_x3y + a_x3z*b_y2z2 + a_y2z2*b_x3z;
    c_arr(42) = a_xy2*b_x2z2 + a_xz2*b_x2y2 + a_x2y2*b_xz2 + a_x2z2*b_xy2 + a_x2y*b_xyz2 + a_xyz2*b_x2y + a_x2z*b_xy2z + a_xy2z*b_x2z + a_x3*b_y2z2 + a_y2z2*b_x3 + a_xyz*b_x2yz + a_x2yz*b_xyz + a_x3y*b_yz2 + a_yz2*b_x3y + a_x3z*b_y2z + a_y2z*b_x3z;
    c_arr(43) = a_x2*b_xy2z + a_xy2*b_x2z + a_x2z*b_xy2 + a_xy2z*b_x2 + a_xz*b_x2y2 + a_x2y2*b_xz + a_xy*b_x2yz + a_x2y*b_xyz + a_xyz*b_x2y + a_x2yz*b_xy + a_x3*b_y2z + a_x3z*b_y2 + a_y2*b_x3z + a_y2z*b_x3 + a_x3y*b_yz + a_yz*b_x3y;
    c_arr(44) = a_x2*b_xy2 + a_xy2*b_x2 + a_x*b_x2y2 + a_x2y2*b_x + a_xy*b_x2y + a_x2y*b_xy + a_x3*b_y2 + a_y2*b_x3 + a_x3y*b_y + a_y*b_x3y;
    c_arr(45) = a_xz3*b_x2yz + a_x2yz*b_xz3 + a_x2z2*b_xyz2 + a_xyz2*b_x2z2 + a_x3z*b_yz3 + a_yz3*b_x3z + a_x3y*b_z4 + a_z4*b_x3y;
    c_arr(46) = a_x2y*b_xz3 + a_xz3*b_x2y + a_x2z*b_xyz2 + a_xz2*b_x2yz + a_x2yz*b_xz2 + a_xyz2*b_x2z + a_x3*b_yz3 + a_yz3*b_x3 + a_x2z2*b_xyz + a_xyz*b_x2z2 + a_x3z*b_yz2 + a_yz2*b_x3z + a_x3y*b_z3 + a_z3*b_x3y;
    c_arr(47) = a_x2*b_xyz2 + a_x2y*b_xz2 + a_xz2*b_x2y + a_xyz2*b_x2 + a_xy*b_x2z2 + a_x2z2*b_xy + a_xz*b_x2yz + a_x2z*b_xyz + a_xyz*b_x2z + a_x2yz*b_xz + a_x3*b_yz2 + a_yz2*b_x3 + a_x3z*b_yz + a_yz*b_x3z + a_x3y*b_z2 + a_z2*b_x3y;
    c_arr(48) = a_x*b_x2yz + a_x2*b_xyz + a_xy*b_x2z + a_xz*b_x2y + a_x2y*b_xz + a_x2z*b_xy + a_xyz*b_x2 + a_x2yz*b_x + a_x3*b_yz + a_x3z*b_y + a_y*b_x3z + a_yz*b_x3 + a_x3y*b_z + a_z*b_x3y;
    c_arr(49) = a_1*b_x3y + a_x3y*b_1 + a_x*b_x2y + a_x2*b_xy + a_xy*b_x2 + a_x2y*b_x + a_x3*b_y + a_y*b_x3;
    c_arr(50) = a_xz3*b_x2z2 + a_x2z2*b_xz3 + a_x3z*b_z4 + a_z4*b_x3z;
    c_arr(51) = a_x2z*b_xz3 + a_xz3*b_x2z + a_xz2*b_x2z2 + a_x2z2*b_xz2 + a_x3*b_z4 + a_z4*b_x3 + a_x3z*b_z3 + a_z3*b_x3z;
    c_arr(52) = a_x2*b_xz3 + a_xz3*b_x2 + a_x2z*b_xz2 + a_xz2*b_x2z + a_xz*b_x2z2 + a_x2z2*b_xz + a_x3*b_z3 + a_z3*b_x3 + a_x3z*b_z2 + a_z2*b_x3z;
    c_arr(53) = a_x2*b_xz2 + a_xz2*b_x2 + a_x*b_x2z2 + a_x2z2*b_x + a_xz*b_x2z + a_x2z*b_xz + a_x3*b_z2 + a_z2*b_x3 + a_x3z*b_z + a_z*b_x3z;
    c_arr(54) = a_1*b_x3z + a_x3z*b_1 + a_x*b_x2z + a_x2*b_xz + a_xz*b_x2 + a_x2z*b_x + a_x3*b_z + a_z*b_x3;
    c_arr(55) = a_1*b_x3 + a_x3*b_1 + a_x*b_x2 + a_x2*b_x;
    c_arr(56) = a_xy3*b_xy3 + a_x2y2*b_y4 + a_y4*b_x2y2;
    c_arr(57) = a_xy3*b_xy2z + a_xy2z*b_xy3 + a_x2yz*b_y4 + a_y4*b_x2yz + a_x2y2*b_y3z + a_y3z*b_x2y2;
    c_arr(58) = a_xy2*b_xy3 + a_xy3*b_xy2 + a_x2y*b_y4 + a_y4*b_x2y + a_x2y2*b_y3 + a_y3*b_x2y2;
    c_arr(59) = a_xy3*b_xyz2 + a_xyz2*b_xy3 + a_x2z2*b_y4 + a_y4*b_x2z2 + a_xy2z*b_xy2z + a_x2y2*b_y2z2 + a_y2z2*b_x2y2 + a_x2yz*b_y3z + a_y3z*b_x2yz;
    c_arr(60) = a_xy3*b_xyz + a_xyz*b_xy3 + a_xy2*b_xy2z + a_xy2z*b_xy2 + a_x2z*b_y4 + a_y4*b_x2z + a_x2y*b_y3z + a_x2yz*b_y3 + a_y3*b_x2yz + a_y3z*b_x2y + a_x2y2*b_y2z + a_y2z*b_x2y2;
    c_arr(61) = a_xy*b_xy3 + a_xy3*b_xy + a_xy2*b_xy2 + a_x2*b_y4 + a_y4*b_x2 + a_x2y*b_y3 + a_y3*b_x2y + a_x2y2*b_y2 + a_y2*b_x2y2;
    c_arr(62) = a_xy3*b_xz3 + a_xz3*b_xy3 + a_xy2z*b_xyz2 + a_xyz2*b_xy2z + a_x2y2*b_yz3 + a_yz3*b_x2y2 + a_x2z2*b_y3z + a_y3z*b_x2z2 + a_x2yz*b_y2z2 + a_y2z2*b_x2yz;
    c_arr(63) = a_xy3*b_xz2 + a_xz2*b_xy3 + a_xy2*b_xyz2 + a_xyz2*b_xy2 + a_x2z2*b_y3 + a_y3*b_x2z2 + a_xyz*b_xy2z + a_xy2z*b_xyz + a_x2z*b_y3z + a_y3z*b_x2z + a_x2y*b_y2z2 + a_x2y2*b_yz2 + a_yz2*b_x2y2 + a_y2z2*b_x2y + a_x2yz*b_y2z + a_y2z*b_x2yz;
    c_arr(64) = a_xz*b_xy3 + a_xy3*b_xz + a_xy*b_xy2z + a_xy2*b_xyz + a_xyz*b_xy2 + a_xy2z*b_xy + a_x2*b_y3z + a_x2z*b_y3 + a_y3*b_x2z + a_y3z*b_x2 + a_x2y*b_y2z + a_x2yz*b_y2 + a_y2*b_x2yz + a_y2z*b_x2y + a_x2y2*b_yz + a_yz*b_x2y2;
    c_arr(65) = a_x*b_xy3 + a_xy3*b_x + a_xy*b_xy2 + a_xy2*b_xy + a_x2*b_y3 + a_y3*b_x2 + a_x2y*b_y2 + a_y2*b_x2y + a_x2y2*b_y + a_y*b_x2y2;
    c_arr(66) = a_xz3*b_xy2z + a_xy2z*b_xz3 + a_xyz2*b_xyz2 + a_x2z2*b_y2z2 + a_y2z2*b_x2z2 + a_x2yz*b_yz3 + a_yz3*b_x2yz + a_x2y2*b_z4 + a_z4*b_x2y2;
    c_arr(67) = a_xy2*b_xz3 + a_xz3*b_xy2 + a_xz2*b_xy2z + a_xy2z*b_xz2 + a_xyz*b_xyz2 + a_xyz2*b_xyz + a_x2y*b_yz3 + a_yz3*b_x2y + a_x2z*b_y2z2 + a_x2z2*b_y2z + a_y2z*b_x2z2 + a_y2z2*b_x2z + a_x2yz*b_yz2 + a_yz2*b_x2yz + a_x2y2*b_z3 + a_z3*b_x2y2;
    c_arr(68) = a_xy2*b_xz2 + a_xz2*b_xy2 + a_xy*b_xyz2 + a_xyz2*b_xy + a_xz*b_xy2z + a_xy2z*b_xz + a_x2*b_y2z2 + a_x2z2*b_y2 + a_y2*b_x2z2 + a_y2z2*b_x2 + a_xyz*b_xyz + a_x2y*b_yz2 + a_yz2*b_x2y + a_x2z*b_y2z + a_y2z*b_x2z + a_x2yz*b_yz + a_yz*b_x2yz + a_x2y2*b_z2 + a_z2*b_x2y2;
    c_arr(69) = a_x*b_xy2z + a_xz*b_xy2 + a_xy2*b_xz + a_xy2z*b_x + a_xy*b_xyz + a_xyz*b_xy + a_x2*b_y2z + a_x2z*b_y2 + a_y2*b_x2z + a_y2z*b_x2 + a_x2y*b_yz + a_x2yz*b_y + a_y*b_x2yz + a_yz*b_x2y + a_x2y2*b_z + a_z*b_x2y2;
    c_arr(70) = a_1*b_x2y2 + a_x2y2*b_1 + a_x*b_xy2 + a_xy2*b_x + a_xy*b_xy + a_x2*b_y2 + a_y2*b_x2 + a_x2y*b_y + a_y*b_x2y;
    c_arr(71) = a_xz3*b_xyz2 + a_xyz2*b_xz3 + a_x2z2*b_yz3 + a_yz3*b_x2z2 + a_x2yz*b_z4 + a_z4*b_x2yz;
    c_arr(72) = a_xz3*b_xyz + a_xyz*b_xz3 + a_xz2*b_xyz2 + a_xyz2*b_xz2 + a_x2z*b_yz3 + a_yz3*b_x2z + a_x2z2*b_yz2 + a_yz2*b_x2z2 + a_x2y*b_z4 + a_z4*b_x2y + a_x2yz*b_z3 + a_z3*b_x2yz;
    c_arr(73) = a_xy*b_xz3 + a_xz3*b_xy + a_xz*b_xyz2 + a_xz2*b_xyz + a_xyz*b_xz2 + a_xyz2*b_xz + a_x2*b_yz3 + a_yz3*b_x2 + a_x2z*b_yz2 + a_yz2*b_x2z + a_x2z2*b_yz + a_yz*b_x2z2 + a_x2y*b_z3 + a_z3*b_x2y + a_x2yz*b_z2 + a_z2*b_x2yz;
    c_arr(74) = a_x*b_xyz2 + a_xy*b_xz2 + a_xz2*b_xy + a_xyz2*b_x + a_xz*b_xyz + a_xyz*b_xz + a_x2*b_yz2 + a_yz2*b_x2 + a_x2z2*b_y + a_y*b_x2z2 + a_x2z*b_yz + a_yz*b_x2z + a_x2y*b_z2 + a_z2*b_x2y + a_x2yz*b_z + a_z*b_x2yz;
    c_arr(75) = a_1*b_x2yz + a_x2yz*b_1 + a_x*b_xyz + a_xy*b_xz + a_xz*b_xy + a_xyz*b_x + a_x2*b_yz + a_x2z*b_y + a_y*b_x2z + a_yz*b_x2 + a_x2y*b_z + a_z*b_x2y;
    c_arr(76) = a_1*b_x2y + a_x2y*b_1 + a_x*b_xy + a_xy*b_x + a_x2*b_y + a_y*b_x2;
    c_arr(77) = a_xz3*b_xz3 + a_x2z2*b_z4 + a_z4*b_x2z2;
    c_arr(78) = a_xz2*b_xz3 + a_xz3*b_xz2 + a_x2z*b_z4 + a_z4*b_x2z + a_x2z2*b_z3 + a_z3*b_x2z2;
    c_arr(79) = a_xz*b_xz3 + a_xz3*b_xz + a_xz2*b_xz2 + a_x2*b_z4 + a_z4*b_x2 + a_x2z*b_z3 + a_z3*b_x2z + a_x2z2*b_z2 + a_z2*b_x2z2;
    c_arr(80) = a_x*b_xz3 + a_xz3*b_x + a_xz*b_xz2 + a_xz2*b_xz + a_x2*b_z3 + a_z3*b_x2 + a_x2z*b_z2 + a_z2*b_x2z + a_x2z2*b_z + a_z*b_x2z2;
    c_arr(81) = a_1*b_x2z2 + a_x2z2*b_1 + a_x*b_xz2 + a_xz2*b_x + a_xz*b_xz + a_x2*b_z2 + a_z2*b_x2 + a_x2z*b_z + a_z*b_x2z;
    c_arr(82) = a_1*b_x2z + a_x2z*b_1 + a_x*b_xz + a_xz*b_x + a_x2*b_z + a_z*b_x2;
    c_arr(83) = a_1*b_x2 + a_x2*b_1 + a_x*b_x;
    c_arr(84) = a_xy3*b_y4 + a_y4*b_xy3;
    c_arr(85) = a_xy3*b_y3z + a_xy2z*b_y4 + a_y4*b_xy2z + a_y3z*b_xy3;
    c_arr(86) = a_xy2*b_y4 + a_xy3*b_y3 + a_y3*b_xy3 + a_y4*b_xy2;
    c_arr(87) = a_xyz2*b_y4 + a_y4*b_xyz2 + a_xy3*b_y2z2 + a_y2z2*b_xy3 + a_xy2z*b_y3z + a_y3z*b_xy2z;
    c_arr(88) = a_xyz*b_y4 + a_y4*b_xyz + a_xy2*b_y3z + a_xy3*b_y2z + a_xy2z*b_y3 + a_y3*b_xy2z + a_y2z*b_xy3 + a_y3z*b_xy2;
    c_arr(89) = a_xy*b_y4 + a_y4*b_xy + a_xy2*b_y3 + a_xy3*b_y2 + a_y2*b_xy3 + a_y3*b_xy2;
    c_arr(90) = a_xz3*b_y4 + a_y4*b_xz3 + a_xy3*b_yz3 + a_yz3*b_xy3 + a_xyz2*b_y3z + a_y3z*b_xyz2 + a_xy2z*b_y2z2 + a_y2z2*b_xy2z;
    c_arr(91) = a_xz2*b_y4 + a_y4*b_xz2 + a_xy3*b_yz2 + a_xyz2*b_y3 + a_y3*b_xyz2 + a_yz2*b_xy3 + a_xy2*b_y2z2 + a_y2z2*b_xy2 + a_xyz*b_y3z + a_y3z*b_xyz + a_xy2z*b_y2z + a_y2z*b_xy2z;
    c_arr(92) = a_xz*b_y4 + a_y4*b_xz + a_xy*b_y3z + a_xy3*b_yz + a_xyz*b_y3 + a_y3*b_xyz + a_yz*b_xy3 + a_y3z*b_xy + a_xy2*b_y2z + a_xy2z*b_y2 + a_y2*b_xy2z + a_y2z*b_xy2;
    c_arr(93) = a_x*b_y4 + a_y4*b_x + a_xy*b_y3 + a_xy3*b_y + a_y*b_xy3 + a_y3*b_xy + a_xy2*b_y2 + a_y2*b_xy2;
    c_arr(94) = a_xz3*b_y3z + a_y3z*b_xz3 + a_xy2z*b_yz3 + a_yz3*b_xy2z + a_xy3*b_z4 + a_z4*b_xy3 + a_xyz2*b_y2z2 + a_y2z2*b_xyz2;
    c_arr(95) = a_xz3*b_y3 + a_y3*b_xz3 + a_xy2*b_yz3 + a_yz3*b_xy2 + a_xz2*b_y3z + a_y3z*b_xz2 + a_xy2z*b_yz2 + a_xyz2*b_y2z + a_y2z*b_xyz2 + a_yz2*b_xy2z + a_xy3*b_z3 + a_z3*b_xy3 + a_xyz*b_y2z2 + a_y2z2*b_xyz;
    c_arr(96) = a_xz2*b_y3 + a_y3*b_xz2 + a_xz*b_y3z + a_xy2*b_yz2 + a_xyz2*b_y2 + a_y2*b_xyz2 + a_yz2*b_xy2 + a_y3z*b_xz + a_xy*b_y2z2 + a_y2z2*b_xy + a_xyz*b_y2z + a_xy2z*b_yz + a_yz*b_xy2z + a_y2z*b_xyz + a_xy3*b_z2 + a_z2*b_xy3;
    c_arr(97) = a_x*b_y3z + a_xz*b_y3 + a_y3*b_xz + a_y3z*b_x + a_xy*b_y2z + a_xy2*b_yz + a_xyz*b_y2 + a_xy2z*b_y + a_y*b_xy2z + a_y2*b_xyz + a_yz*b_xy2 + a_y2z*b_xy + a_xy3*b_z + a_z*b_xy3;
    c_arr(98) = a_1*b_xy3 + a_xy3*b_1 + a_x*b_y3 + a_y3*b_x + a_xy*b_y2 + a_xy2*b_y + a_y*b_xy2 + a_y2*b_xy;
    c_arr(99) = a_xz3*b_y2z2 + a_y2z2*b_xz3 + a_xyz2*b_yz3 + a_yz3*b_xyz2 + a_xy2z*b_z4 + a_z4*b_xy2z;
    c_arr(100) = a_xz3*b_y2z + a_y2z*b_xz3 + a_xz2*b_y2z2 + a_y2z2*b_xz2 + a_xyz*b_yz3 + a_yz3*b_xyz + a_xyz2*b_yz2 + a_yz2*b_xyz2 + a_xy2*b_z4 + a_z4*b_xy2 + a_xy2z*b_z3 + a_z3*b_xy2z;
    c_arr(101) = a_xz3*b_y2 + a_y2*b_xz3 + a_xy*b_yz3 + a_yz3*b_xy + a_xz2*b_y2z + a_y2z*b_xz2 + a_xz*b_y2z2 + a_y2z2*b_xz + a_xyz*b_yz2 + a_xyz2*b_yz + a_yz*b_xyz2 + a_yz2*b_xyz + a_xy2*b_z3 + a_z3*b_xy2 + a_xy2z*b_z2 + a_z2*b_xy2z;
    c_arr(102) = a_xz2*b_y2 + a_y2*b_xz2 + a_x*b_y2z2 + a_y2z2*b_x + a_xy*b_yz2 + a_xyz2*b_y + a_y*b_xyz2 + a_yz2*b_xy + a_xz*b_y2z + a_y2z*b_xz + a_xyz*b_yz + a_yz*b_xyz + a_xy2*b_z2 + a_z2*b_xy2 + a_xy2z*b_z + a_z*b_xy2z;
    c_arr(103) = a_1*b_xy2z + a_xy2z*b_1 + a_x*b_y2z + a_xz*b_y2 + a_y2*b_xz + a_y2z*b_x + a_xy*b_yz + a_xyz*b_y + a_y*b_xyz + a_yz*b_xy + a_xy2*b_z + a_z*b_xy2;
    c_arr(104) = a_1*b_xy2 + a_xy2*b_1 + a_x*b_y2 + a_y2*b_x + a_xy*b_y + a_y*b_xy;
    c_arr(105) = a_xz3*b_yz3 + a_yz3*b_xz3 + a_xyz2*b_z4 + a_z4*b_xyz2;
    c_arr(106) = a_xz2*b_yz3 + a_xz3*b_yz2 + a_yz2*b_xz3 + a_yz3*b_xz2 + a_xyz*b_z4 + a_z4*b_xyz + a_xyz2*b_z3 + a_z3*b_xyz2;
    c_arr(107) = a_xz*b_yz3 + a_xz3*b_yz + a_yz*b_xz3 + a_yz3*b_xz + a_xz2*b_yz2 + a_yz2*b_xz2 + a_xy*b_z4 + a_z4*b_xy + a_xyz*b_z3 + a_z3*b_xyz + a_xyz2*b_z2 + a_z2*b_xyz2;
    c_arr(108) = a_x*b_yz3 + a_xz3*b_y + a_y*b_xz3 + a_yz3*b_x + a_xz*b_yz2 + a_xz2*b_yz + a_yz*b_xz2 + a_yz2*b_xz + a_xy*b_z3 + a_z3*b_xy + a_xyz*b_z2 + a_xyz2*b_z + a_z*b_xyz2 + a_z2*b_xyz;
    c_arr(109) = a_1*b_xyz2 + a_xyz2*b_1 + a_x*b_yz2 + a_xz2*b_y + a_y*b_xz2 + a_yz2*b_x + a_xz*b_yz + a_yz*b_xz + a_xy*b_z2 + a_z2*b_xy + a_xyz*b_z + a_z*b_xyz;
    c_arr(110) = a_1*b_xyz + a_xyz*b_1 + a_x*b_yz + a_xz*b_y + a_y*b_xz + a_yz*b_x + a_xy*b_z + a_z*b_xy;
    c_arr(111) = a_1*b_xy + a_xy*b_1 + a_x*b_y + a_y*b_x;
    c_arr(112) = a_xz3*b_z4 + a_z4*b_xz3;
    c_arr(113) = a_xz2*b_z4 + a_xz3*b_z3 + a_z3*b_xz3 + a_z4*b_xz2;
    c_arr(114) = a_xz*b_z4 + a_z4*b_xz + a_xz2*b_z3 + a_xz3*b_z2 + a_z2*b_xz3 + a_z3*b_xz2;
    c_arr(115) = a_x*b_z4 + a_z4*b_x + a_xz*b_z3 + a_xz3*b_z + a_z*b_xz3 + a_z3*b_xz + a_xz2*b_z2 + a_z2*b_xz2;
    c_arr(116) = a_1*b_xz3 + a_xz3*b_1 + a_x*b_z3 + a_z3*b_x + a_xz*b_z2 + a_xz2*b_z + a_z*b_xz2 + a_z2*b_xz;
    c_arr(117) = a_1*b_xz2 + a_xz2*b_1 + a_x*b_z2 + a_z2*b_x + a_xz*b_z + a_z*b_xz;
    c_arr(118) = a_1*b_xz + a_xz*b_1 + a_x*b_z + a_z*b_x;
    c_arr(119) = a_1*b_x + a_x*b_1;
    c_arr(120) = a_y4*b_y4;
    c_arr(121) = a_y4*b_y3z + a_y3z*b_y4;
    c_arr(122) = a_y3*b_y4 + a_y4*b_y3;
    c_arr(123) = a_y4*b_y2z2 + a_y2z2*b_y4 + a_y3z*b_y3z;
    c_arr(124) = a_y3*b_y3z + a_y4*b_y2z + a_y2z*b_y4 + a_y3z*b_y3;
    c_arr(125) = a_y2*b_y4 + a_y3*b_y3 + a_y4*b_y2;
    c_arr(126) = a_y4*b_yz3 + a_yz3*b_y4 + a_y3z*b_y2z2 + a_y2z2*b_y3z;
    c_arr(127) = a_y4*b_yz2 + a_yz2*b_y4 + a_y3*b_y2z2 + a_y2z2*b_y3 + a_y2z*b_y3z + a_y3z*b_y2z;
    c_arr(128) = a_y4*b_yz + a_yz*b_y4 + a_y2*b_y3z + a_y3*b_y2z + a_y2z*b_y3 + a_y3z*b_y2;
    c_arr(129) = a_y*b_y4 + a_y4*b_y + a_y2*b_y3 + a_y3*b_y2;
    c_arr(130) = a_y3z*b_yz3 + a_yz3*b_y3z + a_y4*b_z4 + a_z4*b_y4 + a_y2z2*b_y2z2;
    c_arr(131) = a_y3*b_yz3 + a_yz3*b_y3 + a_yz2*b_y3z + a_y3z*b_yz2 + a_y4*b_z3 + a_z3*b_y4 + a_y2z*b_y2z2 + a_y2z2*b_y2z;
    c_arr(132) = a_y3*b_yz2 + a_yz2*b_y3 + a_y2*b_y2z2 + a_y2z2*b_y2 + a_yz*b_y3z + a_y3z*b_yz + a_y2z*b_y2z + a_y4*b_z2 + a_z2*b_y4;
    c_arr(133) = a_y*b_y3z + a_y3*b_yz + a_yz*b_y3 + a_y3z*b_y + a_y2*b_y2z + a_y2z*b_y2 + a_y4*b_z + a_z*b_y4;
    c_arr(134) = a_1*b_y4 + a_y4*b_1 + a_y*b_y3 + a_y3*b_y + a_y2*b_y2;
    c_arr(135) = a_yz3*b_y2z2 + a_y2z2*b_yz3 + a_y3z*b_z4 + a_z4*b_y3z;
    c_arr(136) = a_y2z*b_yz3 + a_yz3*b_y2z + a_y3*b_z4 + a_z4*b_y3 + a_yz2*b_y2z2 + a_y2z2*b_yz2 + a_y3z*b_z3 + a_z3*b_y3z;
    c_arr(137) = a_y2*b_yz3 + a_yz3*b_y2 + a_y2z*b_yz2 + a_yz2*b_y2z + a_y3*b_z3 + a_z3*b_y3 + a_yz*b_y2z2 + a_y2z2*b_yz + a_y3z*b_z2 + a_z2*b_y3z;
    c_arr(138) = a_y2*b_yz2 + a_yz2*b_y2 + a_y*b_y2z2 + a_y2z2*b_y + a_yz*b_y2z + a_y2z*b_yz + a_y3*b_z2 + a_z2*b_y3 + a_y3z*b_z + a_z*b_y3z;
    c_arr(139) = a_1*b_y3z + a_y3z*b_1 + a_y*b_y2z + a_y2*b_yz + a_yz*b_y2 + a_y2z*b_y + a_y3*b_z + a_z*b_y3;
    c_arr(140) = a_1*b_y3 + a_y3*b_1 + a_y*b_y2 + a_y2*b_y;
    c_arr(141) = a_yz3*b_yz3 + a_y2z2*b_z4 + a_z4*b_y2z2;
    c_arr(142) = a_yz2*b_yz3 + a_yz3*b_yz2 + a_y2z*b_z4 + a_z4*b_y2z + a_y2z2*b_z3 + a_z3*b_y2z2;
    c_arr(143) = a_yz*b_yz3 + a_yz3*b_yz + a_yz2*b_yz2 + a_y2*b_z4 + a_z4*b_y2 + a_y2z*b_z3 + a_z3*b_y2z + a_y2z2*b_z2 + a_z2*b_y2z2;
    c_arr(144) = a_y*b_yz3 + a_yz3*b_y + a_yz*b_yz2 + a_yz2*b_yz + a_y2*b_z3 + a_z3*b_y2 + a_y2z*b_z2 + a_z2*b_y2z + a_y2z2*b_z + a_z*b_y2z2;
    c_arr(145) = a_1*b_y2z2 + a_y2z2*b_1 + a_y*b_yz2 + a_yz2*b_y + a_yz*b_yz + a_y2*b_z2 + a_z2*b_y2 + a_y2z*b_z + a_z*b_y2z;
    c_arr(146) = a_1*b_y2z + a_y2z*b_1 + a_y*b_yz + a_yz*b_y + a_y2*b_z + a_z*b_y2;
    c_arr(147) = a_1*b_y2 + a_y2*b_1 + a_y*b_y;
    c_arr(148) = a_yz3*b_z4 + a_z4*b_yz3;
    c_arr(149) = a_yz2*b_z4 + a_yz3*b_z3 + a_z3*b_yz3 + a_z4*b_yz2;
    c_arr(150) = a_yz*b_z4 + a_z4*b_yz + a_yz2*b_z3 + a_yz3*b_z2 + a_z2*b_yz3 + a_z3*b_yz2;
    c_arr(151) = a_y*b_z4 + a_z4*b_y + a_yz*b_z3 + a_yz3*b_z + a_z*b_yz3 + a_z3*b_yz + a_yz2*b_z2 + a_z2*b_yz2;
    c_arr(152) = a_1*b_yz3 + a_yz3*b_1 + a_y*b_z3 + a_z3*b_y + a_yz*b_z2 + a_yz2*b_z + a_z*b_yz2 + a_z2*b_yz;
    c_arr(153) = a_1*b_yz2 + a_yz2*b_1 + a_y*b_z2 + a_z2*b_y + a_yz*b_z + a_z*b_yz;
    c_arr(154) = a_1*b_yz + a_yz*b_1 + a_y*b_z + a_z*b_y;
    c_arr(155) = a_1*b_y + a_y*b_1;
    c_arr(156) = a_z4*b_z4;
    c_arr(157) = a_z3*b_z4 + a_z4*b_z3;
    c_arr(158) = a_z2*b_z4 + a_z3*b_z3 + a_z4*b_z2;
    c_arr(159) = a_z*b_z4 + a_z4*b_z + a_z2*b_z3 + a_z3*b_z2;
    c_arr(160) = a_1*b_z4 + a_z4*b_1 + a_z*b_z3 + a_z3*b_z + a_z2*b_z2;
    c_arr(161) = a_1*b_z3 + a_z3*b_1 + a_z*b_z2 + a_z2*b_z;
    c_arr(162) = a_1*b_z2 + a_z2*b_1 + a_z*b_z;
    c_arr(163) = a_1*b_z + a_z*b_1;
    c_arr(164) = a_1*b_1;

    return;
}

void var3_order4_var3_order2_multiplication(Eigen::Matrix<double,1,35>& a_arr, Eigen::Matrix<double,1,10>& b_arr, Eigen::Matrix<double,1,84>& c_arr)
{
    // x^4, x^3*y, x^3*z, x^3, x^2*y^2, x^2*y*z, x^2*y, x^2*z^2, x^2*z, x^2, 
    // x*y^3, x*y^2*z, x*y^2, x*y*z^2, x*y*z, x*y, x*z^3, x*z^2, x*z, x, 
    // y^4, y^3*z, y^3, y^2*z^2, y^2*z, y^2, y*z^3, y*z^2, y*z, y, 
    // z^4, z^3, z^2, z, 1
    double a_x4   = a_arr(0);
    double a_x3y  = a_arr(1);
    double a_x3z  = a_arr(2);
    double a_x3   = a_arr(3);
    double a_x2y2 = a_arr(4);
    double a_x2yz = a_arr(5);
    double a_x2y  = a_arr(6);
    double a_x2z2 = a_arr(7);
    double a_x2z  = a_arr(8);
    double a_x2   = a_arr(9);
    double a_xy3  = a_arr(10);
    double a_xy2z = a_arr(11);
    double a_xy2  = a_arr(12);
    double a_xyz2 = a_arr(13);
    double a_xyz  = a_arr(14);
    double a_xy   = a_arr(15);
    double a_xz3  = a_arr(16);
    double a_xz2  = a_arr(17);
    double a_xz   = a_arr(18);
    double a_x    = a_arr(19);
    double a_y4   = a_arr(20);
    double a_y3z  = a_arr(21);
    double a_y3   = a_arr(22);
    double a_y2z2 = a_arr(23);
    double a_y2z  = a_arr(24);
    double a_y2   = a_arr(25);
    double a_yz3  = a_arr(26);
    double a_yz2  = a_arr(27);
    double a_yz   = a_arr(28);
    double a_y    = a_arr(29);
    double a_z4   = a_arr(30);
    double a_z3   = a_arr(31);
    double a_z2   = a_arr(32);
    double a_z    = a_arr(33);
    double a_1    = a_arr(34);

    // x^2, x*y, x*z, x, y^2, y*z, y, z^2, z, 1
    double b_x2 = b_arr(0);
    double b_xy = b_arr(1);
    double b_xz = b_arr(2);
    double b_x  = b_arr(3);
    double b_y2 = b_arr(4);
    double b_yz = b_arr(5);
    double b_y  = b_arr(6);
    double b_z2 = b_arr(7);
    double b_z  = b_arr(8);
    double b_1  = b_arr(9);

    // results
    // x^6, x^5*y, x^5*z, x^5, x^4*y^2, x^4*y*z, x^4*y, x^4*z^2, x^4*z, x^4, 
    // x^3*y^3, x^3*y^2*z, x^3*y^2, x^3*y*z^2, x^3*y*z, x^3*y, x^3*z^3, x^3*z^2, x^3*z, x^3, 
    // x^2*y^4, x^2*y^3*z, x^2*y^3, x^2*y^2*z^2, x^2*y^2*z, x^2*y^2, x^2*y*z^3, x^2*y*z^2, x^2*y*z, x^2*y, 
    // x^2*z^4, x^2*z^3, x^2*z^2, x^2*z, x^2, x*y^5, x*y^4*z, x*y^4, x*y^3*z^2, x*y^3*z, 
    // x*y^3, x*y^2*z^3, x*y^2*z^2, x*y^2*z, x*y^2, x*y*z^4, x*y*z^3, x*y*z^2, x*y*z, x*y, 
    // x*z^5, x*z^4, x*z^3, x*z^2, x*z, x, y^6, y^5*z, y^5, y^4*z^2, 
    // y^4*z, y^4, y^3*z^3, y^3*z^2, y^3*z, y^3, y^2*z^4, y^2*z^3, y^2*z^2, y^2*z, 
    // y^2, y*z^5, y*z^4, y*z^3, y*z^2, y*z, y, z^6, z^5, z^4, 
    // z^3, z^2, z, 1
 
    c_arr(0) = a_x4*b_x2;
    c_arr(1) = a_x4*b_xy + a_x3y*b_x2;
    c_arr(2) = a_x2y2*b_x2 + a_x3y*b_xy + a_x4*b_y2;
    c_arr(3) = a_xy3*b_x2 + a_x2y2*b_xy + a_x3y*b_y2;
    c_arr(4) = a_xy3*b_xy + a_y4*b_x2 + a_x2y2*b_y2;
    c_arr(5) = a_y4*b_xy + a_xy3*b_y2;
    c_arr(6) = a_y4*b_y2;
    c_arr(7) = a_x4*b_xz + a_x3z*b_x2;
    c_arr(8) = a_x3y*b_xz + a_x3z*b_xy + a_x2yz*b_x2 + a_x4*b_yz;
    c_arr(9) = a_xy2z*b_x2 + a_x2y2*b_xz + a_x2yz*b_xy + a_x3z*b_y2 + a_x3y*b_yz;
    c_arr(10) = a_xy3*b_xz + a_xy2z*b_xy + a_y3z*b_x2 + a_x2yz*b_y2 + a_x2y2*b_yz;
    c_arr(11) = a_y4*b_xz + a_xy3*b_yz + a_y3z*b_xy + a_xy2z*b_y2;
    c_arr(12) = a_y4*b_yz + a_y3z*b_y2;
    c_arr(13) = a_x2z2*b_x2 + a_x3z*b_xz + a_x4*b_z2;
    c_arr(14) = a_xyz2*b_x2 + a_x2z2*b_xy + a_x2yz*b_xz + a_x3z*b_yz + a_x3y*b_z2;
    c_arr(15) = a_xyz2*b_xy + a_xy2z*b_xz + a_x2z2*b_y2 + a_y2z2*b_x2 + a_x2yz*b_yz + a_x2y2*b_z2;
    c_arr(16) = a_xyz2*b_y2 + a_y3z*b_xz + a_y2z2*b_xy + a_xy2z*b_yz + a_xy3*b_z2;
    c_arr(17) = a_y2z2*b_y2 + a_y3z*b_yz + a_y4*b_z2;
    c_arr(18) = a_xz3*b_x2 + a_x2z2*b_xz + a_x3z*b_z2;
    c_arr(19) = a_xz3*b_xy + a_xyz2*b_xz + a_yz3*b_x2 + a_x2z2*b_yz + a_x2yz*b_z2;
    c_arr(20) = a_xz3*b_y2 + a_yz3*b_xy + a_y2z2*b_xz + a_xyz2*b_yz + a_xy2z*b_z2;
    c_arr(21) = a_yz3*b_y2 + a_y2z2*b_yz + a_y3z*b_z2;
    c_arr(22) = a_xz3*b_xz + a_z4*b_x2 + a_x2z2*b_z2;
    c_arr(23) = a_xz3*b_yz + a_yz3*b_xz + a_z4*b_xy + a_xyz2*b_z2;
    c_arr(24) = a_yz3*b_yz + a_z4*b_y2 + a_y2z2*b_z2;
    c_arr(25) = a_z4*b_xz + a_xz3*b_z2;
    c_arr(26) = a_z4*b_yz + a_yz3*b_z2;
    c_arr(27) = a_z4*b_z2;
    c_arr(28) = a_x4*b_x + a_x3*b_x2;
    c_arr(29) = a_x3*b_xy + a_x3y*b_x + a_x2y*b_x2 + a_x4*b_y;
    c_arr(30) = a_xy2*b_x2 + a_x2y2*b_x + a_x2y*b_xy + a_x3*b_y2 + a_x3y*b_y;
    c_arr(31) = a_xy3*b_x + a_xy2*b_xy + a_y3*b_x2 + a_x2y*b_y2 + a_x2y2*b_y;
    c_arr(32) = a_y4*b_x + a_xy3*b_y + a_y3*b_xy + a_xy2*b_y2;
    c_arr(33) = a_y4*b_y + a_y3*b_y2;
    c_arr(34) = a_x3*b_xz + a_x3z*b_x + a_x2z*b_x2 + a_x4*b_z;
    c_arr(35) = a_x2y*b_xz + a_x2z*b_xy + a_xyz*b_x2 + a_x2yz*b_x + a_x3*b_yz + a_x3z*b_y + a_x3y*b_z;
    c_arr(36) = a_xy2*b_xz + a_xy2z*b_x + a_xyz*b_xy + a_x2z*b_y2 + a_y2z*b_x2 + a_x2y*b_yz + a_x2yz*b_y + a_x2y2*b_z;
    c_arr(37) = a_y3*b_xz + a_y3z*b_x + a_xy2*b_yz + a_xyz*b_y2 + a_xy2z*b_y + a_y2z*b_xy + a_xy3*b_z;
    c_arr(38) = a_y3*b_yz + a_y3z*b_y + a_y2z*b_y2 + a_y4*b_z;
    c_arr(39) = a_xz2*b_x2 + a_x2z2*b_x + a_x2z*b_xz + a_x3*b_z2 + a_x3z*b_z;
    c_arr(40) = a_xz2*b_xy + a_xyz2*b_x + a_xyz*b_xz + a_yz2*b_x2 + a_x2z2*b_y + a_x2z*b_yz + a_x2y*b_z2 + a_x2yz*b_z;
    c_arr(41) = a_xz2*b_y2 + a_y2z2*b_x + a_xyz2*b_y + a_yz2*b_xy + a_y2z*b_xz + a_xyz*b_yz + a_xy2*b_z2 + a_xy2z*b_z;
    c_arr(42) = a_yz2*b_y2 + a_y2z2*b_y + a_y2z*b_yz + a_y3*b_z2 + a_y3z*b_z;
    c_arr(43) = a_xz3*b_x + a_xz2*b_xz + a_z3*b_x2 + a_x2z*b_z2 + a_x2z2*b_z;
    c_arr(44) = a_xz3*b_y + a_yz3*b_x + a_xz2*b_yz + a_yz2*b_xz + a_z3*b_xy + a_xyz*b_z2 + a_xyz2*b_z;
    c_arr(45) = a_yz3*b_y + a_yz2*b_yz + a_z3*b_y2 + a_y2z*b_z2 + a_y2z2*b_z;
    c_arr(46) = a_z4*b_x + a_xz3*b_z + a_z3*b_xz + a_xz2*b_z2;
    c_arr(47) = a_z4*b_y + a_yz3*b_z + a_z3*b_yz + a_yz2*b_z2;
    c_arr(48) = a_z4*b_z + a_z3*b_z2;
    c_arr(49) = a_x4*b_1 + a_x3*b_x + a_x2*b_x2;
    c_arr(50) = a_x3y*b_1 + a_x2*b_xy + a_xy*b_x2 + a_x2y*b_x + a_x3*b_y;
    c_arr(51) = a_x2y2*b_1 + a_xy2*b_x + a_xy*b_xy + a_x2*b_y2 + a_y2*b_x2 + a_x2y*b_y;
    c_arr(52) = a_xy3*b_1 + a_y3*b_x + a_xy*b_y2 + a_xy2*b_y + a_y2*b_xy;
    c_arr(53) = a_y4*b_1 + a_y3*b_y + a_y2*b_y2;
    c_arr(54) = a_x3z*b_1 + a_x2*b_xz + a_xz*b_x2 + a_x2z*b_x + a_x3*b_z;
    c_arr(55) = a_x2yz*b_1 + a_xy*b_xz + a_xz*b_xy + a_xyz*b_x + a_x2*b_yz + a_x2z*b_y + a_yz*b_x2 + a_x2y*b_z;
    c_arr(56) = a_xy2z*b_1 + a_xz*b_y2 + a_y2*b_xz + a_y2z*b_x + a_xy*b_yz + a_xyz*b_y + a_yz*b_xy + a_xy2*b_z;
    c_arr(57) = a_y3z*b_1 + a_y2*b_yz + a_yz*b_y2 + a_y2z*b_y + a_y3*b_z;
    c_arr(58) = a_x2z2*b_1 + a_xz2*b_x + a_xz*b_xz + a_x2*b_z2 + a_z2*b_x2 + a_x2z*b_z;
    c_arr(59) = a_xyz2*b_1 + a_xz2*b_y + a_yz2*b_x + a_xz*b_yz + a_yz*b_xz + a_xy*b_z2 + a_z2*b_xy + a_xyz*b_z;
    c_arr(60) = a_y2z2*b_1 + a_yz2*b_y + a_yz*b_yz + a_y2*b_z2 + a_z2*b_y2 + a_y2z*b_z;
    c_arr(61) = a_xz3*b_1 + a_z3*b_x + a_xz*b_z2 + a_xz2*b_z + a_z2*b_xz;
    c_arr(62) = a_yz3*b_1 + a_z3*b_y + a_yz*b_z2 + a_yz2*b_z + a_z2*b_yz;
    c_arr(63) = a_z4*b_1 + a_z3*b_z + a_z2*b_z2;
    c_arr(64) = a_x3*b_1 + a_x*b_x2 + a_x2*b_x;
    c_arr(65) = a_x2y*b_1 + a_x*b_xy + a_xy*b_x + a_x2*b_y + a_y*b_x2;
    c_arr(66) = a_xy2*b_1 + a_x*b_y2 + a_y2*b_x + a_xy*b_y + a_y*b_xy;
    c_arr(67) = a_y3*b_1 + a_y*b_y2 + a_y2*b_y;
    c_arr(68) = a_x2z*b_1 + a_x*b_xz + a_xz*b_x + a_x2*b_z + a_z*b_x2;
    c_arr(69) = a_xyz*b_1 + a_x*b_yz + a_xz*b_y + a_y*b_xz + a_yz*b_x + a_xy*b_z + a_z*b_xy;
    c_arr(70) = a_y2z*b_1 + a_y*b_yz + a_yz*b_y + a_y2*b_z + a_z*b_y2;
    c_arr(71) = a_xz2*b_1 + a_x*b_z2 + a_z2*b_x + a_xz*b_z + a_z*b_xz;
    c_arr(72) = a_yz2*b_1 + a_y*b_z2 + a_z2*b_y + a_yz*b_z + a_z*b_yz;
    c_arr(73) = a_z3*b_1 + a_z*b_z2 + a_z2*b_z;
    c_arr(74) = a_1*b_x2 + a_x2*b_1 + a_x*b_x;
    c_arr(75) = a_1*b_xy + a_xy*b_1 + a_x*b_y + a_y*b_x;
    c_arr(76) = a_1*b_y2 + a_y2*b_1 + a_y*b_y;
    c_arr(77) = a_1*b_xz + a_xz*b_1 + a_x*b_z + a_z*b_x;
    c_arr(78) = a_1*b_yz + a_yz*b_1 + a_y*b_z + a_z*b_y;
    c_arr(79) = a_1*b_z2 + a_z2*b_1 + a_z*b_z;
    c_arr(80) = a_1*b_x + a_x*b_1;
    c_arr(81) = a_1*b_y + a_y*b_1;
    c_arr(82) = a_1*b_z + a_z*b_1;
    c_arr(83) = a_1*b_1;

    return;
}

void var3_order2_three_multiplication(
    Eigen::Matrix<double,1,10>& a_arr, Eigen::Matrix<double,1,10>& b_arr, 
    Eigen::Matrix<double,1,10>& c_arr, Eigen::Matrix<double,1,84>& d_arr)
{
    Eigen::Matrix<double,1,35> d1;
    var3_order2_two_multiplication(a_arr, b_arr, d1);
    var3_order4_var3_order2_multiplication(d1, c_arr, d_arr);
    return;
}

void var3_order2_four_multiplication(
    Eigen::Matrix<double,1,10>& a_arr, Eigen::Matrix<double,1,10>& b_arr, 
    Eigen::Matrix<double,1,10>& c_arr, Eigen::Matrix<double,1,10>& d_arr,
    Eigen::Matrix<double,1,165>& e_arr)
{
    Eigen::Matrix<double,1,35> e1, e2;
    var3_order2_two_multiplication(a_arr, b_arr, e1);
    var3_order2_two_multiplication(c_arr, d_arr, e2);
    var3_order4_multiplication(e1, e2, e_arr);
    return;
}
