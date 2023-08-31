# **Data-Driven Discovery of Invariant Measures**

This repository contains MATLAB scripts to reproduce the data and figures from [Data-driven discovery of invariant measures](https://arxiv.org/abs/2308.15318) by [Jason J. Bramburger](https://hybrid.concordia.ca/jbrambur/) and [Giovanni Fantuzzi](https://dcn.nat.fau.eu/giovanni-fantuzzi/) (2023).

## **Paper Abstract**
Invariant measures encode the long-time behaviour of a dynamical system. In this work, we propose an optimization-based method to discover invariant measures directly from data gathered from a system. Our method does not require an explicit model for the dynamics and allows one to target specific invariant measures, such as physical and ergodic measures. Moreover, it applies to both deterministic and stochastic dynamics in either continuous or discrete time. We provide convergence results and illustrate the performance of our method on data from the logistic map and a stochastic double-well system, for which invariant measures can be found by other means. We then use our method to approximate the physical measure of the chaotic attractor of the Rossler system, and we extract unstable periodic orbits embedded in this attractor by identifying discrete-time periodic points of a suitably defined Poincare map. This final example is truly data-driven and shows that our method can significantly outperform previous approaches based on model identification.

## **Required MATLAB Packages**
All scripts require YALMIP and MOSEK to run. Both packages can be download for free at 
- YALMIP: https://yalmip.github.io/download/
- MOSEK: https://www.mosek.com/downloads/

In some scripts we further use the ChebFun package, which can be downloaded at: https://www.chebfun.org/

## **Repository Contents**
This repository contains MATLAB script to reproduce the results for the examples in Section 6 of the paper. Precisely, the scripts perform the following tasks:
- `logistic_physical.m` discovers the physical measure of the logistic map from data. This script accompanies Section 6.1.3.
- `logistic_ergodic.m` discovers ergodic measures of the logistic map that are approximately supported on sets of unstable periodic orbits. This script accompanies Section 6.1.2.
- `double_well.m` approximates the unique invariant measure from data of a stochastic double-well system. This script accompanies Section 6.2.
- `rossler_physical.m` discovers the physical measure of the Rossler system in a chaotic parameter regime from data. The MATLAB data `Rossler_LongQuadAvgs.mat` contains approximations of the quadratic moments from long-time averages of the system that can be compared with predictions from the discovered measure. This script accompanies Section 6.3.1.
- `rossler_poincare.m` discovers ergodic measures of the Poincare map of the Rossler system that are approximately supported on the unstable periodic orbits in the Poincare section. This script accompanies Section 6.3.2.

Additional scripts that are shared among those above are stored in the [Auxiliary Scripts](https://github.com/jbramburger/data-measures/tree/main/Auxiliary%20Scripts) folder.
