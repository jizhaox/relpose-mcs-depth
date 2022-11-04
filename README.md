# README
This package implements relative pose estimation for multi-camera systems from two affine correspondences, including two multi-camera system solvers (`solver_depth_inter_2ac`, `solver_depth_intra_2ac`) and a single camera solver (`solver_depth_mono_2ac`).

Source codes and Matlab mex files with demo code are provided in the package. The core solvers are written by C++. Matlab mex files are compiled using Ubuntu 16.04 + Matlab R2019a. Run `test_solver_AC.m` in folder "test".

# Reference

[1] Banglei Guan, and Ji Zhao. [**Relative Pose Estimation for Multi-Camera Systems from Two Affine Correspondences**](https://www.ecva.net/papers/eccv_2022/papers_ECCV/html/5358_ECCV_2022_paper.php). European Conference on Computer Vision, 2022. (Oral)

If you use this package in an academic work, please cite:

    @inproceedings{guan2022relative,
      title={Relative Pose Estimation for Multi-Camera Systems from Two Affine Correspondences},
      author={Guan, Banglei and Zhao, Ji},
      booktitle={European Conference on Computer Vision},
      year={2022}
     }


# solver_depth_inter_2ac

A minimal solver for the relative pose estimation of multi-camera systems using two inter-camera affine correspondences. Returns a maximum of 48 solutions or 56 solutions.
* **Solver**:  `solver_depth_inter_2ac.mexa64`   

* **API**: `[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_inter_2ac(Image1(1:2,:), Image2(1:2,:), At, R_cam, t_cam, 0);` `[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_inter_2ac(Image1(1:2,:), Image2(1:2,:), At, R_cam, t_cam, 1);`

* **Input data for Demo**: 

     `Image1` (3\*2 matrix): normalized homogeneous image coordinates of two inter-camera feature points expressed in view 1.

     `Image2` (3\*2 matrix): normalized homogeneous image coordinates of two inter-camera feature points expressed in view 2.

     `At` (2\*2\*2 matrix): `At(:,:,1) = At1`, `At(:,:,2) = At2`, where `At1` and `At2` are the corresponding 2\*2 local affine transformations of two inter-camera feature points.

     `R_cam` (3\*3\*2 matrix): extrinsic rotation of two cameras expressed in the reference of the multi-camera system.

     `t_cam` (3\*2 matrix): extrinsic translation of two cameras expressed in the reference of the multi-camera system.

     `Option`: 0 refers the `2AC-inter-48 solver`, 1 refers the `2AC-inter-56 solver`.

* **Output data for Demo**: 

     `cay_sols` (3\*N matrix): real number solutions of Cayley rotation parameter, N is the number of real number solutions.

     `t_sols` (3\*N matrix): real number solutions of translation.

     `R_sols` (3\*3\*N matrix): real number solutions of rotation matrix.

     `cay_sols_all` (3\*M matrix): all solutions of Cayley rotation parameter, including real number solutions and complex number solutions. M is the number of all the solutions.


# solver_depth_intra_2ac

A minimal solver for the relative pose estimation of multi-camera systems using two intra-camera affine correspondences. Returns a maximum of 48 solutions.
* **Solver**:  `solver_depth_intra_2ac.mexa64`

* **API**: `[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_intra_2ac(Image1(1:2,:), Image2(1:2,:), At, R_cam, t_cam);`

* **Input data for Demo**: 

     `Image1` (3\*2 matrix): normalized homogeneous image coordinates of two intra-camera feature points expressed in view 1.

     `Image2` (3\*2 matrix): normalized homogeneous image coordinates of two intra-camera feature points expressed in view 2.

     `At` (2\*2\*2 matrix): `At(:,:,1) = At1`, `At(:,:,2) = At2`, where `At1` and `At2` are the corresponding 2\*2 local affine transformations of two intra-camera feature points.

     `R_cam` (3\*3\*2 matrix): extrinsic rotation of two cameras expressed in the reference of the multi-camera system.

     `t_cam` (3\*2 matrix): extrinsic translation of two cameras expressed in the reference of the multi-camera system.

* **Output data for Demo**: 

     `cay_sols` (3\*N matrix): real number solutions of Cayley rotation parameter, N is the number of real number solutions.

     `t_sols` (3\*N matrix): real number solutions of translation.

     `R_sols` (3\*3\*N matrix): real number solutions of rotation matrix.

     `cay_sols_all` (3\*M matrix): all solutions of Cayley rotation parameter, including real number solutions and complex number solutions. M is the number of all the solutions.


# solver_depth_mono_2ac

A minimal solver for the relative pose estimation of a single camera using two affine correspondences. Returns a maximum of 20 solutions.
* **Solver**: `solver_depth_mono_2ac.mexa64`

* **API**: `[cay_sols, t_sols, R_sols, cay_sols_all] = solver_depth_mono_2ac(Image1(1:2,:), Image2(1:2,:), At);`

* **Input data for Demo**: 

     `Image1` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in view 1.

     `Image2` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in view 2.

     `At` (2\*2\*2 matrix): `At(:,:,1) = At1`, `At(:,:,2) = At2`, where `At1` and `At2` are the corresponding 2\*2 local affine transformations of two feature points.

* **Output data for Demo**: 

     `cay_sols` (3\*N matrix): real number solutions of Cayley rotation parameter, N is the number of real number solutions.

     `t_sols` (3\*N matrix): real number solutions of translation.

     `R_sols` (3\*3\*N matrix): real number solutions of rotation matrix.

     `cay_sols_all` (3\*M matrix): all solutions of Cayley rotation parameter, including real number solutions and complex number solutions. M is the number of all the solutions.


# Run

Compiled files using Ubuntu 16.04 + Matlab R2019a are provided. You can run the package in Matlab.

` test_solver_AC.m` is the demo which shows how to call the APIs.



