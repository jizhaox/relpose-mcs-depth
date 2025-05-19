#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "util.h"
#include "common_stew_core.h"

using namespace std;
using namespace Eigen;

void construct_order6_poly_stew(std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M, 
    Eigen::Matrix<double,10,84>& C)
{
    C.setZero();

    Eigen::Matrix<int,10,3> Idx_all;
    Idx_all << 
     0, 1, 2, 
     0, 1, 3, 
     0, 1, 4,
     0, 2, 3,
     0, 2, 4,
     0, 3, 4,
     1, 2, 3,
     1, 2, 4,
     1, 3, 4,
     2, 3, 4;
    
    for (int i = 0; i < 10; i++)
    {
        Eigen::Matrix<int,1,3> idx;
        idx = Idx_all.row(i);

        std::vector<Eigen::Matrix<double,1,10>> m1, m2, m3;
        Eigen::Matrix<double,1,10> 
            m11, m12, m13, 
            m21, m22, m23,
            m31, m32, m33;

        m1 = M[idx(0)];
        m2 = M[idx(1)];
        m3 = M[idx(2)];
        m11 = m1[0]; m12 = m1[1]; m13 = m1[2];
        m21 = m2[0]; m22 = m2[1]; m23 = m2[2];
        m31 = m3[0]; m32 = m3[1]; m33 = m3[2];

        Eigen::Matrix<double,1,84> d, d0;
        d.setZero();
        d0.setZero();
        var3_order2_three_multiplication(m11, m22, m33, d0); d = d + d0;
        var3_order2_three_multiplication(m11, m23, m32, d0); d = d - d0;
        var3_order2_three_multiplication(m12, m21, m33, d0); d = d - d0;
        var3_order2_three_multiplication(m12, m23, m31, d0); d = d + d0;
        var3_order2_three_multiplication(m13, m21, m32, d0); d = d + d0;
        var3_order2_three_multiplication(m13, m22, m31, d0); d = d - d0;

        C.block(i, 0, 1, 84) = d;
    }
    return;
}

void construct_order4_extra_poly_stew(std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M, 
    int idx_r_1, int idx_r_2, 
    Eigen::Matrix<double,3,84>& C)
{
    C.setZero();

    Eigen::Matrix<int,3,2> Idx_all;
    Idx_all << 
     0, 1, 
     0, 2, 
     1, 2;
    
    Eigen::Matrix<double,1,10> m_one;
    m_one.setZero(); m_one(9) = 1;
    for (int i = 0; i < 3; i++)
    {
        Eigen::Matrix<int,1,2> idx;
        idx = Idx_all.row(i);

        std::vector<Eigen::Matrix<double,1,10>> m1, m2;
        Eigen::Matrix<double,1,10> m11, m12, m21, m22;

        m1 = M[idx_r_1];
        m2 = M[idx_r_2];
        m11 = m1[idx(0)]; m12 = m1[idx(1)];
        m21 = m2[idx(0)]; m22 = m2[idx(1)];

        Eigen::Matrix<double,1,84> d, d0;
        d.setZero();
        d0.setZero();
        var3_order2_three_multiplication(m11, m22, m_one, d0); d = d + d0;
        var3_order2_three_multiplication(m12, m21, m_one, d0); d = d - d0;

        C.block(i, 0, 1, 84) = d;
    }
    return;
}


void create_coeffs_stew(double* coeffs, double* input, 
    std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M, 
    std::vector<Eigen::Matrix<double,6,1>>& Line_i_all, std::vector<Eigen::Matrix<double,6,1>>& Line_j_all,
    double* input_Image_1, double* input_Image_2, double* input_affine_tran, 
    double* extrinsic_R_camera, double* extrinsic_T_camera, AC_TYPE actype, bool is_known_angle)
{
    M.clear();
    std::vector<Eigen::Matrix3d> R_camera, Ac;
    std::vector<Eigen::Vector3d> T_camera, Image1, Image2;

    format_convert(input_Image_1, input_Image_2, input_affine_tran, extrinsic_R_camera, extrinsic_T_camera, 
        Image1, Image2, Ac, R_camera, T_camera);

    int point_num = 2;
    Eigen::Matrix<double,20,84> C;
    Eigen::Matrix<double,6,84> C_extra;
    C.setZero();
    C_extra.setZero();
    for (int k = 0; k < point_num; k++)
    {
        Eigen::Matrix<double,6,1> Line_i, Line_j;
        std::vector<Eigen::Matrix<double,1,10>> f_pt;
        std::vector<Eigen::Matrix<double,1,10>> f_af1, f_af2;
        bool pt_choose_as_origin = false;
        // 1st element is k; 2nd to last elements are other elements except for k.
        int ele_order[2];
        if (k==0)
        {
            ele_order[0] = 0;
            ele_order[1] = 1;
        }
        else if (k==1)
        {
            ele_order[0] = 1;
            ele_order[1] = 0;
        }

        for (int n = 0; n < point_num; n++)
        {
            int i = ele_order[n];

            // extract observations
            Eigen::Vector3d P1 = Image1[i];
            Eigen::Vector3d P2 = Image2[i];
            Eigen::Matrix3d Atemp = Ac[i];
            Eigen::Vector3d U1 = P1;
            U1.normalize();
            Eigen::Vector3d U2 = P2;
            U2.normalize();

            int idx1 = 0;
            int idx2 = 0;
            if (actype==INTER_CAM_CONSTRAINT_FULL || actype==INTER_CAM_CONSTRAINT_PARTIAL)
            {
                if (i==0)
                {
                    idx1 = 0;
                    idx2 = 1;
                }
                else
                {
                    idx1 = 1;
                    idx2 = 0;
                }
            }
            else if (actype==INTRA_CAM_CONSTRAINT_FULL || actype==INTRA_CAM_CONSTRAINT_PARTIAL)
            {
                if (i==0)
                {
                    idx1 = 0;
                    idx2 = 0;
                }
                else
                {
                    idx1 = 1;
                    idx2 = 1;
                }
            }
            else
            {
                std::cout << "error: unknown AC type!" << std::endl;
                return;
            }

            // extract extrinsic parameters
            Eigen::Matrix3d R1 = R_camera[idx1];
            Eigen::Vector3d T1 = T_camera[idx1];
            Eigen::Matrix3d R2 = R_camera[idx2];
            Eigen::Vector3d T2 = T_camera[idx2];

            if (n==0) // this point is choose as orign
            {
                Eigen::Vector3d V = R1*U1;
                Line_i.block(0, 0, 3, 1) = V;
                Line_i.block(3, 0, 3, 1) = T1.cross(V);
                V = R2*U2;
                Line_j.block(0, 0, 3, 1) = V;
                Line_j.block(3, 0, 3, 1) = T2.cross(V);

                Line_i_all.push_back(Line_i);
                Line_j_all.push_back(Line_j);

                pt_choose_as_origin = true;
            }
            else
            {
                pt_choose_as_origin = false;
            }

            if (i==0)
            {
                create_coefficient_core(f_pt, f_af1, Line_i, Line_j, R1, R2, T1, T2, P1, P2, Atemp, pt_choose_as_origin);
            }
            else if (i==1)
            {
                create_coefficient_core(f_pt, f_af2, Line_i, Line_j, R1, R2, T1, T2, P1, P2, Atemp, pt_choose_as_origin);
            }
        }

        std::vector<std::vector<Eigen::Matrix<double,1,10>>> M0;
        std::vector<Eigen::Matrix<double,1,10>> m1;
        // row 1
        m1.clear();
        m1.push_back(f_pt[0]);
        m1.push_back(f_pt[1]);
        m1.push_back(f_pt[2]);
        M0.push_back(m1);
        M.push_back(m1);
        // row 2
        m1.clear();
        m1.push_back(f_af1[0]);
        m1.push_back(f_af1[1]);
        m1.push_back(f_af1[2]);
        M0.push_back(m1);
        M.push_back(m1);
        // row 3
        m1.clear();
        m1.push_back(f_af1[3]);
        m1.push_back(f_af1[4]);
        m1.push_back(f_af1[5]);
        M0.push_back(m1);
        M.push_back(m1);
        // row 4
        m1.clear();
        m1.push_back(f_af2[0]);
        m1.push_back(f_af2[1]);
        m1.push_back(f_af2[2]);
        M0.push_back(m1);
        M.push_back(m1);
        // row 5
        m1.clear();
        m1.push_back(f_af2[3]);
        m1.push_back(f_af2[4]);
        m1.push_back(f_af2[5]);
        M0.push_back(m1);
        M.push_back(m1);

        Eigen::Matrix<double,10,84> C0;
        construct_order6_poly_stew(M0, C0);
        C.block(k*10, 0, 10, 84) = C0;

        if (actype == INTER_CAM_CONSTRAINT_FULL || actype == INTRA_CAM_CONSTRAINT_FULL)
        {
            Eigen::Matrix<double,3,84> C_extra0;
            C_extra0.setZero();
            int idx_r_1, idx_r_2;
            if (k==0)
            {
                idx_r_1 = 1;
                idx_r_2 = 2;
            }
            else if (k==1)
            {
                idx_r_1 = 5+3;
                idx_r_2 = 5+4;
            }
            construct_order4_extra_poly_stew(M, idx_r_1, idx_r_2, C_extra0);
            C_extra.block(k*3, 0, 3, 84) = C_extra0;
        }

    }

    // prepare data for Matlab interface
    // Matlab memory is column-major order
    int subset[8] = {0, 1, 3, 6, 10, 11, 13, 16};
    int cnt = 0;
    if (actype == INTER_CAM_CONSTRAINT_PARTIAL || actype == INTRA_CAM_CONSTRAINT_PARTIAL)
    {
        for (int j = 0; j < 84; j++)
        {
            for (int i = 0; i < 20; i++)
            {
                coeffs[cnt] = C(i, j);
                cnt++;
            }
        }
    }
    if (actype == INTER_CAM_CONSTRAINT_FULL || actype == INTRA_CAM_CONSTRAINT_FULL)
    {
        for (int j = 0; j < 84; j++)
        {
            for (int i = 0; i < 20; i++)
            {
                coeffs[cnt] = C(i, j);
                cnt++;
            }
            for (int i = 0; i < 6; i++)
            {
                coeffs[cnt] = C_extra(i, j);
                cnt++;
            }
        }
    }

    // prepare data for the solver
    cnt = 0;
    if (!is_known_angle)
    {
        if (actype == INTER_CAM_CONSTRAINT_PARTIAL)
        {
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 84; j++)
                {
                    input[cnt] = C(i, j);
                    cnt++;
                }
            }
        }
        if (actype == INTER_CAM_CONSTRAINT_FULL)
        {
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 84; j++)
                {
                    input[cnt] = C(i, j);
                    cnt++;
                }
            }
            for (int i = 0; i < 6; i++)
            {
                for (int j = 49; j < 84; j++)
                {
                    input[cnt] = C_extra(i, j);
                    cnt++;
                }
            }
        }
        if (actype == INTRA_CAM_CONSTRAINT_FULL)
        {
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 83; j++)
                {
                    input[cnt] = C(i, j);
                    cnt++;
                }
            }
            for (int i = 0; i < 6; i++)
            {
                for (int j = 49; j < 84; j++)
                {
                    input[cnt] = C_extra(i, j);
                    cnt++;
                }
            }
        }
    }
    else
    {
        if (actype == INTER_CAM_CONSTRAINT_PARTIAL){
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 84; j++)
                {
                    input[cnt] = C(subset[i], j);
                    cnt++;
                }
            }
        }
        if (actype == INTER_CAM_CONSTRAINT_FULL){
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 84; j++)
                {
                    input[cnt] = C(subset[i], j);
                    cnt++;
                }
            }
            cnt+=4;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 49; j < 84; j++)
                {
                    input[cnt] = C_extra(i, j);
                    cnt++;
                }
            }
        }
        if (actype == INTRA_CAM_CONSTRAINT_PARTIAL){
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 83; j++)
                {
                    input[cnt] = C(subset[i], j);
                    cnt++;
                }
            }
        }
        if (actype == INTRA_CAM_CONSTRAINT_FULL){
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < 83; j++)
                {
                    input[cnt] = C(subset[i], j);
                    cnt++;
                }
            }
            cnt+=4;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 49; j < 84; j++)
                {
                    input[cnt] = C_extra(i, j);
                    cnt++;
                }
            }
        }
    }
    
    return;
}

void calculate_M_5_by_3(Eigen::Matrix<double, 5, 3>& M_double, std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M, double x, double y, double z)
{
    int n_row = M.size();
    if (n_row != 10)
    {
        std::cout << "size of M has an error!" << std::endl;
        return;
    }
    for (int i = 0; i < n_row; i++)
    {
        int n_col =  M[i].size();
        if (n_col != 3)
        {
            std::cout << "size of M has an error!" << std::endl;
            return;
        }
    }
    
    Eigen::Matrix<double, 10, 1> m;
    m(0) = x*x;
    m(1) = x*y;
    m(2) = x*z;
    m(3) = x;
    m(4) = y*y;
    m(5) = y*z;
    m(6) = y;
    m(7) = z*z;
    m(8) = z;
    m(9) = 1;

    // use the first 5*3 matrix
    for (int i = 0; i < 5; i++)
    {
        std::vector<Eigen::Matrix<double,1,10>> m1 = M[i];
        for (int j = 0; j < 3; j++)
        {
            Eigen::Matrix<double,1,1> rslt= m1[j]*m;
            M_double(i,j) = rslt(0);
        }
    }

    return;
}

void calculate_translation_stew(
    Eigen::MatrixXcd sols, std::vector<std::vector<Eigen::Matrix<double,1,10>>>& M, 
    std::vector<Eigen::Matrix<double,6,1>>& Line_i_all, std::vector<Eigen::Matrix<double,6,1>>& Line_j_all, 
    std::vector<Eigen::Matrix<double,3,1>>& q_arr, std::vector<Eigen::Matrix<double,3,1>>& t_arr, 
    std::vector<Eigen::Matrix<double,3,3>>& rotm, bool is_known_angle)
{
    q_arr.clear();
    t_arr.clear();
    rotm.clear();
    if (sols.rows() != 3)
    {
        std::cout << "size of solution has an error!" << std::endl;
        return;
    }

    Eigen::Matrix<double,6,1> Line_i = Line_i_all[0];
    Eigen::Matrix<double,6,1> Line_j = Line_j_all[0];
    Eigen::Vector3d U1 = Line_i.block(0, 0, 3, 1);
    Eigen::Vector3d V1 = Line_i.block(3, 0, 3, 1);
    Eigen::Vector3d Cross_U1V1 = U1.cross(V1);

    Eigen::Vector3d U2 = Line_j.block(0, 0, 3, 1);
    Eigen::Vector3d V2 = Line_j.block(3, 0, 3, 1);
    Eigen::Vector3d Cross_U2V2 = U2.cross(V2);

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

        Eigen::Matrix<double,3,3> Rwf2;
        cayley2rotm(Rwf2, q);

        Eigen::Matrix<double, 5, 3> M_double;
        calculate_M_5_by_3(M_double, M, x, y, z);

        Eigen::Matrix<double, 2, 1> Lambda;
        if (!is_known_angle)
        {
            Eigen::Matrix<double, 5, 2> C0 = M_double.block(0, 0, 5, 2);
            Eigen::Matrix<double, 5, 1> C1 = M_double.block(0, 2, 5, 1);
            Lambda = -C0.colPivHouseholderQr().solve(C1);
        }
        else
        {
            Eigen::Matrix<double, 4, 2> C0 = M_double.block(0, 0, 4, 2);
            Eigen::Matrix<double, 4, 1> C1 = M_double.block(0, 2, 4, 1);
            Lambda = -C0.colPivHouseholderQr().solve(C1);
        }

        double lambda1 = Lambda(0);
        double lambda2 = Lambda(1);

        Eigen::Vector3d Twf1 = Cross_U1V1 + lambda1*U1;
        Eigen::Vector3d Twf2 = Cross_U2V2 + lambda2*U2;

        rotm.push_back(Rwf2);
        t_arr.push_back(-Rwf2*Twf1 + Twf2);
    }
    return;
}
