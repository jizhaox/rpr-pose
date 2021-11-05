## 1. Overview

This code is an implementation of the following paper. The core solvers are written by C++. An Matlab Mex interface is provided for easy-to-use.

```
@article{zhao2020minimal,
  title={Minimal Case Relative Pose Computation Using Ray-Point-Ray Features},
  author={Zhao, Ji and Kneip, Laurent and He, Yijia and Ma, Jiayi},
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
  year={2020},
  volume={42},
  number={5},
  pages={1176--1190}
}
```

In this paper, we proposed a few solvers for relative pose estimation from `ray-point-ray (RPR)` features.

We also proposed a few new `five-point methods` for relative pose estimation, which is a classical problem in geometric vision. 

Note: We have tried 4 different rotation representations, and found that the Cayley representation has the best performance. Hence, only the solvers with the Cayley representation are provided.



**Authors:** [Ji Zhao](https://sites.google.com/site/drjizhao)


## 2. Quick Start

Run `compile.m` in subfolders of "`solvers_5pt_rpr`" to compile the mex files. Make sure to set the proper path of Eigen library in `compile.m`. This step can be skipped once the mex files have been compiled. Compiled files using Ubuntu 16.04 and Matlab R2019a are provided.

Run `test_solver*.m` in folder "`test`"

