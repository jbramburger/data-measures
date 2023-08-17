# **Data-Driven Discovery of Invariant Measures**

This repository contains MATLAB scripts to reproduce the data and figures from [Data-driven discovery of invariant measures](TBD) by [Jason J. Bramburger](https://hybrid.concordia.ca/jbrambur/) and [Giovanni Fantuzzi](https://dcn.nat.fau.eu/giovanni-fantuzzi/) (2023).

## **Paper Abstract**
Coming soon.

## **Required MATLAB Packages**
All scripts require YALMIP and MOSEK to run. Both packages can be download for free at 
- YALMIP: https://yalmip.github.io/download/
- MOSEK: https://www.mosek.com/downloads/

In some scripts we further use the ChebFun package, which can be downloaded at: https://www.chebfun.org/

## **Repository Contents**
This repository contains MATLAB script to reproduce the results for the examples in Section 5 of the paper. Precisely, the scripts perform the following tasks:
- logistic_physical.m discovers the physical measure of the logistic map from data. This script accompanies Section 5.1.3.
- logistic_ergodic.m discovers ergodic measures of the logistic map that are approximately supported on sets of unstable periodic orbits. This script accompanies Section 5.1.2.
- double_well.m approximates the unique invariant measure from data of a stochastic double-well system. This script accompanies Section 5.2.
- rossler_physical.m discovers the physical measure of the Rossler system in a chaotic parameter regime from data. The MATLAB data Rossler_LongQuadAvgs.mat contains approximations of the quadratic moments from long-time averages of the system that can be compared with predictions from the discovered measure. This script accompanies Section 5.3.1.
- rossler_poincare.m discovers ergodic measures of the Poincare map of the Rossler system that are approximately supported on the unstable periodic orbits in the Poincare section. This script accompanies Section 5.3.2.

Additional scripts that are shared among those above are stored in the [Auxiliary Scripts](https://github.com/jbramburger/data-measures/tree/main/Auxiliary%20Scripts) folder.
