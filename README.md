# ADMM-Code
Matlab and C++ Code for ADMM Simulation

## Matlab Code
Original ADMM implementation in Matlab is located in the Matlab folder

## C++ Code
ADMM C++ stand alone implementation is located in the CppApp folder

## To Do
* Remove Matlab libraries from Project 3
* Rename Project 3 folder
* Test BFGS performance

## Better performance for standalone C++ program with Clang

The key idea is adding the optimization flag -Ofast or -O3

clang -fopenmp -Ofast -IC:\Users\TBSUser\repos\eigen -IC:\Users\TBSUser\repos\ADMM-Code\LBFGSpp\include -I"C:\Program Files\MATLAB\R2020b\extern\include" -L"C:\Program Files\MATLAB\R2020b\bin\win64" -L"C:\Program Files\MATLAB\R2020b\extern\lib\win64\microsoft" -llibmat -llibmx -llibmex -o myprog.exe mainmex.cpp Tet_Energy_Term.cpp Triangle_Energy_Term.cpp admm_solver.cpp 

