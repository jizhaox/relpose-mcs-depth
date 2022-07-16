#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "util.h"

using namespace std;
using namespace Eigen;

void build_coefficient_matrix_depth_parameter(
    std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M,
    std::vector<Eigen::Vector2d>& Image1, std::vector<Eigen::Vector2d>& Image2, 
    std::vector<Eigen::Matrix2d>& Ac, int i)
{
    std::vector<Eigen::Matrix<double,1,10>> r1, r2, r3, r4, r5;
    // reference is built on i-th AC
    int j = -1;
    if (i == 0)
        j = 1;
    else
        j = 0;
    
    // 1st AC
    // from world reference to frame 1
    double u11 = Image1[i](0);
    double u12 = Image1[i](1);
    double u21 = Image2[i](0);
    double u22 = Image2[i](1);

    double a1 = Ac[i](0, 0);
    double a2 = Ac[i](0, 1);
    double a3 = Ac[i](1, 0);
    double a4 = Ac[i](1, 1);

    Eigen::Matrix<double,1,10> f1_New_C11, f1_New_C12, f2_New_C21, f2_New_C22, f3_New_C31, f3_New_C32;
    f2_New_C21 << u22 - u12, -2*u21, 2*u12*u21, - 2*u12*u22 - 2, - u12 - u22, 2*u12*u22 - 2, 2*u12*u21, u12 + u22, 2*u21, u12 - u22;
    f2_New_C22 << a1*u12 + a3*u11 - a1*u22 + a3*u21, 2*a3*u12 - 2*a1*u11, 2*a3 + 2*a1*u11*u22 - 2*a3*u11*u21, 2*a1 + 2*a1*u12*u22 - 2*a3*u12*u21, a3*u21 - a3*u11 - a1*u22 - a1*u12, 2*a1*u12*u22 - 2*a1 - 2*a3*u12*u21, 2*a3 - 2*a1*u11*u22 + 2*a3*u11*u21, a1*u12 - a3*u11 + a1*u22 - a3*u21, - 2*a1*u11 - 2*a3*u12, a3*u11 - a1*u12 + a1*u22 - a3*u21;
    f3_New_C31 << u11 + u21, 2*u22, 2 - 2*u11*u21, 2*u11*u22, u11 - u21, -2*u11*u22, - 2*u11*u21 - 2, - u11 - u21, 2*u22, u21 - u11;
    f3_New_C32 << a2*u12 + a4*u11 - a2*u22 + a4*u21, 2*a4*u12 - 2*a2*u11, 2*a4 + 2*a2*u11*u22 - 2*a4*u11*u21, 2*a2 + 2*a2*u12*u22 - 2*a4*u12*u21, a4*u21 - a4*u11 - a2*u22 - a2*u12, 2*a2*u12*u22 - 2*a2 - 2*a4*u12*u21, 2*a4 - 2*a2*u11*u22 + 2*a4*u11*u21, a2*u12 - a4*u11 + a2*u22 - a4*u21, - 2*a2*u11 - 2*a4*u12, a4*u11 - a2*u12 + a2*u22 - a4*u21;

    r1.push_back(f2_New_C21);
    r1.push_back(f2_New_C22);
    r2.push_back(f3_New_C31);
    r2.push_back(f3_New_C32);

    // 2nd AC
    // from world reference to frame 1
    double v11 = Image1[j](0);
    double v12 = Image1[j](1);
    double v21 = Image2[j](0);
    double v22 = Image2[j](1);

    a1 = Ac[j](0, 0);
    a2 = Ac[j](0, 1);
    a3 = Ac[j](1, 0);
    a4 = Ac[j](1, 1);

    f1_New_C11 << u11*v12 - u12*v11 - u11*v22 - u12*v21 + v11*v22 + v12*v21, 2*u11*v21 - 2*u12*v22 - 2*v11*v21 + 2*v12*v22, 2*v12 - 2*u12 - 2*u11*v12*v21 + 2*u12*v11*v21, 2*u11 - 2*v11 + 2*u11*v12*v22 - 2*u12*v11*v22, u11*v12 - u12*v11 + u11*v22 + u12*v21 - v11*v22 - v12*v21, 2*u11 - 2*v11 - 2*u11*v12*v22 + 2*u12*v11*v22, 2*u12 - 2*v12 - 2*u11*v12*v21 + 2*u12*v11*v21, u12*v11 - u11*v12 - u11*v22 + u12*v21 + v11*v22 - v12*v21, 2*v11*v21 - 2*u12*v22 - 2*u11*v21 + 2*v12*v22, u12*v11 - u11*v12 + u11*v22 - u12*v21 - v11*v22 + v12*v21;
    f1_New_C12 << u21*v22 - u22*v11 - u21*v12 - u22*v21 + v11*v22 + v12*v21, 2*u21*v11 - 2*u22*v12 - 2*v11*v21 + 2*v12*v22, 2*v22 - 2*u22 - 2*u21*v11*v22 + 2*u22*v11*v21, 2*v21 - 2*u21 - 2*u21*v12*v22 + 2*u22*v12*v21, u21*v12 + u22*v11 + u21*v22 - u22*v21 - v11*v22 - v12*v21, 2*u21 - 2*v21 - 2*u21*v12*v22 + 2*u22*v12*v21, 2*v22 - 2*u22 + 2*u21*v11*v22 - 2*u22*v11*v21, u22*v11 - u21*v12 - u21*v22 + u22*v21 - v11*v22 + v12*v21, 2*u21*v11 + 2*u22*v12 - 2*v11*v21 - 2*v12*v22, u21*v12 - u22*v11 - u21*v22 + u22*v21 + v11*v22 - v12*v21;
    f2_New_C21 << v22 - u12 - a1*u12 - a3*u11 + a1*v12 + a3*v11, 2*a1*u11 - 2*v21 - 2*a3*u12 - 2*a1*v11 + 2*a3*v12, 2*u12*v21 - 2*a1*u11*v12 + 2*a1*u12*v11, 2*a3*u11*v12 - 2*u12*v22 - 2*a3*u12*v11 - 2, a1*u12 - v22 - u12 + a3*u11 - a1*v12 - a3*v11, 2*u12*v22 - 2*a3*u11*v12 + 2*a3*u12*v11 - 2, 2*u12*v21 - 2*a1*u11*v12 + 2*a1*u12*v11, u12 + v22 + a1*u12 - a3*u11 - a1*v12 + a3*v11, 2*v21 - 2*a1*u11 - 2*a3*u12 + 2*a1*v11 + 2*a3*v12, u12 - v22 - a1*u12 + a3*u11 + a1*v12 - a3*v11;
    f2_New_C22 << v22 - u22 - a1*u22 + a3*u21 + a1*v12 + a3*v11, 2*u21 - 2*v21 - 2*a1*v11 + 2*a3*v12, 2*a3 - 2*u21*v22 + 2*u22*v21 + 2*a1*u22*v11 - 2*a3*u21*v11, 2*a1 + 2*a1*u22*v12 - 2*a3*u21*v12, u22 - v22 - a1*u22 + a3*u21 - a1*v12 - a3*v11, 2*a1*u22*v12 - 2*a1 - 2*a3*u21*v12, 2*a3 + 2*u21*v22 - 2*u22*v21 - 2*a1*u22*v11 + 2*a3*u21*v11, u22 - v22 + a1*u22 - a3*u21 + a1*v12 - a3*v11, 2*u21 - 2*v21 - 2*a1*v11 - 2*a3*v12, v22 - u22 + a1*u22 - a3*u21 - a1*v12 + a3*v11;
    f3_New_C31 << u11 + v21 - a2*u12 - a4*u11 + a2*v12 + a4*v11, 2*v22 + 2*a2*u11 - 2*a4*u12 - 2*a2*v11 + 2*a4*v12, 2*a2*u12*v11 - 2*a2*u11*v12 - 2*u11*v21 + 2, 2*u11*v22 + 2*a4*u11*v12 - 2*a4*u12*v11, u11 - v21 + a2*u12 + a4*u11 - a2*v12 - a4*v11, 2*a4*u12*v11 - 2*a4*u11*v12 - 2*u11*v22, 2*a2*u12*v11 - 2*a2*u11*v12 - 2*u11*v21 - 2, a2*u12 - v21 - u11 - a4*u11 - a2*v12 + a4*v11, 2*v22 - 2*a2*u11 - 2*a4*u12 + 2*a2*v11 + 2*a4*v12, v21 - u11 - a2*u12 + a4*u11 + a2*v12 - a4*v11;
    f3_New_C32 << v21 - u21 - a2*u22 + a4*u21 + a2*v12 + a4*v11, 2*v22 - 2*u22 - 2*a2*v11 + 2*a4*v12, 2*a4 + 2*a2*u22*v11 - 2*a4*u21*v11, 2*a2 - 2*u21*v22 + 2*u22*v21 + 2*a2*u22*v12 - 2*a4*u21*v12, u21 - v21 - a2*u22 + a4*u21 - a2*v12 - a4*v11, 2*u22*v21 - 2*u21*v22 - 2*a2 + 2*a2*u22*v12 - 2*a4*u21*v12, 2*a4 - 2*a2*u22*v11 + 2*a4*u21*v11, v21 - u21 + a2*u22 - a4*u21 + a2*v12 - a4*v11, 2*u22 - 2*v22 - 2*a2*v11 - 2*a4*v12, u21 - v21 + a2*u22 - a4*u21 - a2*v12 + a4*v11;

    r3.push_back(f1_New_C11);
    r3.push_back(f1_New_C12);
    r4.push_back(f2_New_C21);
    r4.push_back(f2_New_C22);
    r5.push_back(f3_New_C31);
    r5.push_back(f3_New_C32);

    M.push_back(r3);
    if (i==0)
    {
        M.push_back(r1);
        M.push_back(r2);
        M.push_back(r4);
        M.push_back(r5);
    }
    else
    {
        M.push_back(r4);
        M.push_back(r5);
        M.push_back(r1);
        M.push_back(r2);
    }

    return;
}

void build_coefficient_matrix_qxqyqz(Eigen::MatrixXd& C, 
    std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M, Eigen::MatrixXi Idx_all, int idx_ac)
{
    if (Idx_all.rows()==0 || Idx_all.cols()!=2)
    {
        std::cerr << "size of index array is wrong!" << std::endl;
    }
    // all the 2*2 sub-determinants of M(r) must equal zero. 
    // This gives C_3^2 = 3 equations which only involve the rotation parameters.
    int row_offset_M = 5*idx_ac; // each AC introduces 5 equations, corresponding 5 rows in M
    int n_det_equ_per_ac = Idx_all.rows();

    for (int i = 0; i < n_det_equ_per_ac; i++)
    {
        int i0 = Idx_all(i,0);
        int i1 = Idx_all(i,1);

        Eigen::Matrix<double,1,10> M11, M12, M21, M22;
        M11 = M[row_offset_M+i0][0];
        M12 = M[row_offset_M+i0][1];
        M21 = M[row_offset_M+i1][0];
        M22 = M[row_offset_M+i1][1];
        double b1 = M11(0);
        double b2 = M11(1);
        double b3 = M11(2);
        double b4 = M11(3);
        double b5 = M11(4);
        double b6 = M11(5);
        double b7 = M11(6);
        double b8 = M11(7);
        double b9 = M11(8);
        double b10 = M11(9);

        double c1 = M12(0);
        double c2 = M12(1);
        double c3 = M12(2);
        double c4 = M12(3);
        double c5 = M12(4);
        double c6 = M12(5);
        double c7 = M12(6);
        double c8 = M12(7);
        double c9 = M12(8);
        double c10 = M12(9);

        double d1 = M21(0);
        double d2 = M21(1);
        double d3 = M21(2);
        double d4 = M21(3);
        double d5 = M21(4);
        double d6 = M21(5);
        double d7 = M21(6);
        double d8 = M21(7);
        double d9 = M21(8);
        double d10 = M21(9);

        double f1 = M22(0);
        double f2 = M22(1);
        double f3 = M22(2);
        double f4 = M22(3);
        double f5 = M22(4);
        double f6 = M22(5);
        double f7 = M22(6);
        double f8 = M22(7);
        double f9 = M22(8);
        double f10 = M22(9);

        C.row(n_det_equ_per_ac*idx_ac+i) << b1*f1 - c1*d1, b1*f2 - c2*d1 - c1*d2 + b2*f1, b1*f3 - c3*d1 - c1*d3 + b3*f1, b1*f4 - c4*d1 - c1*d4 + b4*f1, b2*f2 - c1*d5 - c5*d1 - c2*d2 + b1*f5 + b5*f1, b2*f3 - c3*d2 - c1*d6 - c6*d1 - c2*d3 + b3*f2 + b1*f6 + b6*f1, b2*f4 - c4*d2 - c1*d7 - c7*d1 - c2*d4 + b4*f2 + b1*f7 + b7*f1, b3*f3 - c1*d8 - c8*d1 - c3*d3 + b1*f8 + b8*f1, b3*f4 - c4*d3 - c1*d9 - c9*d1 - c3*d4 + b4*f3 + b1*f9 + b9*f1, b4*f4 - c1*d10 - c10*d1 - c4*d4 + b1*f10 + b10*f1, b2*f5 - c5*d2 - c2*d5 + b5*f2, b2*f6 - c3*d5 - c5*d3 - c6*d2 - c2*d6 + b3*f5 + b5*f3 + b6*f2, b2*f7 - c4*d5 - c5*d4 - c7*d2 - c2*d7 + b4*f5 + b5*f4 + b7*f2, b3*f6 - c6*d3 - c2*d8 - c8*d2 - c3*d6 + b6*f3 + b2*f8 + b8*f2, b3*f7 - c4*d6 - c6*d4 - c7*d3 - c2*d9 - c9*d2 - c3*d7 + b4*f6 + b6*f4 + b7*f3 + b2*f9 + b9*f2, b4*f7 - c7*d4 - c2*d10 - c10*d2 - c4*d7 + b7*f4 + b2*f10 + b10*f2, b3*f8 - c8*d3 - c3*d8 + b8*f3, b3*f9 - c4*d8 - c8*d4 - c9*d3 - c3*d9 + b4*f8 + b8*f4 + b9*f3, b3*f10 - c4*d9 - c9*d4 - c10*d3 - c3*d10 + b4*f9 + b9*f4 + b10*f3, b4*f10 - c10*d4 - c4*d10 + b10*f4, b5*f5 - c5*d5, b5*f6 - c6*d5 - c5*d6 + b6*f5, b5*f7 - c7*d5 - c5*d7 + b7*f5, b6*f6 - c5*d8 - c8*d5 - c6*d6 + b5*f8 + b8*f5, b6*f7 - c7*d6 - c5*d9 - c9*d5 - c6*d7 + b7*f6 + b5*f9 + b9*f5, b7*f7 - c5*d10 - c10*d5 - c7*d7 + b5*f10 + b10*f5, b6*f8 - c8*d6 - c6*d8 + b8*f6, b6*f9 - c7*d8 - c8*d7 - c9*d6 - c6*d9 + b7*f8 + b8*f7 + b9*f6, b6*f10 - c7*d9 - c9*d7 - c10*d6 - c6*d10 + b7*f9 + b9*f7 + b10*f6, b7*f10 - c10*d7 - c7*d10 + b10*f7, b8*f8 - c8*d8, b8*f9 - c9*d8 - c8*d9 + b9*f8, b8*f10 - c9*d9 - c10*d8 - c8*d10 + b9*f9 + b10*f8, b9*f10 - c10*d9 - c9*d10 + b10*f9, b10*f10 - c10*d10;
    }

    return;
}

void create_coeffs_stew_single(double* coeffs, double* input, 
    std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M, 
    double* input_Image_1, double* input_Image_2, double* input_affine_tran)
{
    Eigen::MatrixXi Idx_all(6,2);
    Idx_all << 
        0, 1, 
        0, 2,
        0, 3,
        1, 2,
        1, 3, 
        2, 3;

    M.clear();
    std::vector<Eigen::Matrix2d> Ac;
    std::vector<Eigen::Vector2d> Image1, Image2;

    format_convert(input_Image_1, input_Image_2, input_affine_tran, Image1, Image2, Ac);
    Eigen::MatrixXd C(12,35);
    Eigen::MatrixXd C2(12,35);
    for (int i = 0; i < 2; i++)
    {
        build_coefficient_matrix_depth_parameter(M, Image1, Image2, Ac, i);
        build_coefficient_matrix_qxqyqz(C2, M, Idx_all, i);
    }

    // change from plex to grevlex order
    // source oder
    // x^4, x^3*y, x^3*z, x^3, x^2*y^2, x^2*y*z, x^2*y, x^2*z^2, x^2*z, x^2, x*y^3, x*y^2*z, x*y^2, x*y*z^2, x*y*z, x*y, x*z^3, x*z^2, x*z, x, y^4, y^3*z, y^3, y^2*z^2, y^2*z, y^2, y*z^3, y*z^2, y*z, y, z^4, z^3, z^2, z, 1
    // target order
    // x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^3*z, x^2*y*z, x*y^2*z, y^3*z, x^2*z^2, x*y*z^2, y^2*z^2, x*z^3, y*z^3, z^4, x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3, x^2, x*y, y^2, x*z, y*z, z^2, x, y, z, 1
    Eigen::Matrix<int,1,35> Index;
    Index << 0, 1, 4, 10, 20, 2, 5, 11, 21, 7, 13, 23, 16, 26, 30, 3, 6, 12, 22, 8, 14, 24, 17, 27, 31, 9, 15, 25, 18, 28, 32, 19, 29, 33, 34;
    for (int i = 0; i < 35; i++)
    {
        C.col(i) = C2.col(Index(i));
    }

    // prepare data for Matlab interface
    // Matlab memory is column-major order
    int subset[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    if (true){
        for (int i = 0; i < 12; i++)
        {
            for (int j = 0; j < 35; j++)
            {
                coeffs[j*12+subset[i]] = C(subset[i], j);
            }
        }
    }
    // prepare data for the solver
    int cnt = 0;
    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 35; j++)
        {
            input[cnt] = C(subset[i], j);
            cnt++;
        }
    }
    return;
}

void create_coeffs_stew_single_ka(double* coeffs, double* input, 
    std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M, 
    double* input_Image_1, double* input_Image_2, double* input_affine_tran)
{
    Eigen::MatrixXi Idx_all(3,2);
    Idx_all << 
        0, 1, 
        0, 2, 
        1, 2;

    M.clear();
    std::vector<Eigen::Matrix2d> Ac;
    std::vector<Eigen::Vector2d> Image1, Image2;

    format_convert(input_Image_1, input_Image_2, input_affine_tran, Image1, Image2, Ac);
    Eigen::MatrixXd C(6,35);
    Eigen::MatrixXd C2(6,35);
    for (int i = 0; i < 2; i++)
    {
        build_coefficient_matrix_depth_parameter(M, Image1, Image2, Ac, i);
        build_coefficient_matrix_qxqyqz(C2, M, Idx_all, i);
    }

    // change from plex to grevlex order
    // source oder
    // x^4, x^3*y, x^3*z, x^3, x^2*y^2, x^2*y*z, x^2*y, x^2*z^2, x^2*z, x^2, x*y^3, x*y^2*z, x*y^2, x*y*z^2, x*y*z, x*y, x*z^3, x*z^2, x*z, x, y^4, y^3*z, y^3, y^2*z^2, y^2*z, y^2, y*z^3, y*z^2, y*z, y, z^4, z^3, z^2, z, 1
    // target order
    // x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^3*z, x^2*y*z, x*y^2*z, y^3*z, x^2*z^2, x*y*z^2, y^2*z^2, x*z^3, y*z^3, z^4, x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3, x^2, x*y, y^2, x*z, y*z, z^2, x, y, z, 1
    Eigen::Matrix<int,1,35> Index;
    Index << 0, 1, 4, 10, 20, 2, 5, 11, 21, 7, 13, 23, 16, 26, 30, 3, 6, 12, 22, 8, 14, 24, 17, 27, 31, 9, 15, 25, 18, 28, 32, 19, 29, 33, 34;
    for (int i = 0; i < 35; i++)
    {
        C.col(i) = C2.col(Index(i));
    }

    // prepare data for Matlab interface
    // Matlab memory is column-major order
    int subset[6] = {0, 1, 2, 3, 4, 5};
    if (true){
        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 35; j++)
            {
                coeffs[j*7+subset[i]] = C(subset[i], j);
            }
        }
    }
    // prepare data for the solver
    int cnt = 0;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 35; j++)
        {
            input[cnt] = C(subset[i], j);
            cnt++;
        }
    }
    return;
}

void calculate_translation_stew_single(
    Eigen::MatrixXcd& sols, std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M, 
    double *Image_1, double *Image_2, 
    std::vector<Eigen::Matrix<double,3,1>>& q_arr, std::vector<Eigen::Matrix<double,3,1>>& t_arr, 
    std::vector<Eigen::Matrix<double,3,3>>& rotm, bool is_known_angle)
{
    q_arr.clear();
    t_arr.clear();
    rotm.clear();
    if (sols.rows() != 3)
    {
        std::cerr << "size of solution has an error!" << std::endl;
        return;
    }

    for (int j = 0; j < sols.cols(); j++)
    {
        std::complex<double> cx, cy, cz;
        double x, y, z;
        cx = sols(0, j);
        cy = sols(1, j);
        cz = sols(2, j);

        if (abs(cx.imag()) > NEAR_ZERO_THRESHOLD || abs(cy.imag()) > NEAR_ZERO_THRESHOLD || abs(cy.imag()) > NEAR_ZERO_THRESHOLD)
            continue;

        x = cx.real();
        y = cy.real();
        z = cz.real();

        Eigen::Vector3d q;
        q << x, y, z;
        q_arr.push_back(q);

        Eigen::Matrix<double,3,3> Rwf2;
        quad2rotm(Rwf2, q);

        Eigen::Matrix<double, 1, 10> qxqyqz;
        qxqyqz << x*x, x*y, x*z, x, y*y, y*z, y, z*z, z, 1.0;

        Eigen::Matrix<double, 2, 2> V;
        if (is_known_angle)
        {
            Eigen::Matrix<double, 3, 2> C_Lam1_Lam2;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    C_Lam1_Lam2(i,j) = qxqyqz.dot(M[i][j]);
                }
            }
            JacobiSVD<Eigen::Matrix<double, 3, 2>> svd(C_Lam1_Lam2, ComputeThinU | ComputeThinV);
             V = svd.matrixV();
        }
        else
        {
            Eigen::Matrix<double, 5, 2> C_Lam1_Lam2;
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    C_Lam1_Lam2(i,j) = qxqyqz.dot(M[i][j]);
                }
            }
            JacobiSVD<Eigen::Matrix<double, 5, 2>> svd(C_Lam1_Lam2, ComputeThinU | ComputeThinV);
             V = svd.matrixV();
        }
        Eigen::Matrix<double, 2, 1> Lambda = V.rightCols(1);

        double lambda1 = Lambda(0);
        double lambda2 = Lambda(1);
        
        Eigen::Vector3d U1, U2;
        U1 << *Image_1, *(Image_1+1), 1.0;
        U2 << *Image_2, *(Image_2+1), 1.0;
        Eigen::Vector3d Twf1 = lambda1*U1;
        Eigen::Vector3d Twf2 = lambda2*U2;

        rotm.push_back(Rwf2);
        Eigen::Vector3d Tf1f2 = -Rwf2*Twf1 + Twf2;
        Tf1f2.normalize();
        t_arr.push_back(Tf1f2);
    }

    return;
}

