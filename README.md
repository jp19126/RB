# Randomized Benchmarking
This repository collects the python programmes realizing the Clifford-based Randomized Benchmarking technique in Quantum Control problems, determining the error performance of a given quantum gate.

## Version
v1. (Historical version)Rewrote the generation rule of random unitary gates from the Clifford group. No other changes were applied, which means they remained the same as the original MATLAB programme.

v2. The programme is now able to compute the fidelity after applying m Clifford group gates to the initial state. Currently, we only considered a single iteration. In the next version, the same process will be repeated to compute the average fidelity. Meanwhile, we will modify m and do curve fitting based on that.
v2.1 Function get_para simplified.

## References
[1] Barends, R., Kelly, J., Megrant, A., Veitia, A., Sank, D., Jeffrey, E., White, T.C., Mutus, J., Fowler, A.G., Campbell, B. and Chen, Y., 2014. Superconducting quantum circuits at the surface code threshold for fault tolerance. Nature, 508(7497), p.500.

[2] Knill, E., Leibfried, D., Reichle, R., Britton, J., Blakestad, R.B., Jost, J.D., Langer, C., Ozeri, R., Seidelin, S. and Wineland, D.J., 2008. Randomized benchmarking of quantum gates. Physical Review A, 77(1), p.012307.

[3] Ryan, C.A., Laforest, M. and Laflamme, R., 2009. Randomized benchmarking of single-and multi-qubit control in liquid-state NMR quantum information processing. New Journal of Physics, 11(1), p.013034.

[4] Brown, K.R., Wilson, A.C., Colombe, Y., Ospelkaus, C., Meier, A.M., Knill, E., Leibfried, D. and Wineland, D.J., 2011. Single-qubit-gate error below 10− 4 in a trapped ion. Physical Review A, 84(3), p.030303.

[5] Córcoles, A.D., Gambetta, J.M., Chow, J.M., Smolin, J.A., Ware, M., Strand, J., Plourde, B.L. and Steffen, M., 2013. Process verification of two-qubit quantum gates by randomized benchmarking. Physical Review A, 87(3), p.030301.

## Contributors
- WANG Rui
- LIU Baojie
