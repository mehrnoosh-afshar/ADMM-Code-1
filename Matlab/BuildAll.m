clearvars;

sourceFile1 = ['../CppApp/Project3/Project3/main.cpp'];
sourceFile2= ['../CppApp/Project3/Project3/Triangle_Energy_Term.cpp'];
sourceFile3= ['../CppApp/Project3/Project3/Triangle_Energy_Term.h'];
sourceFile4= ['../CppApp/Project3/Project3/admm_solver.cpp'];
sourceFile5 = ['../CppApp/Project3/Project3/admm_solver.h'];
fprintf('Building %s \n', sourceFile1)
% otherFlags = '-I/Applications/MATLAB_R2020b.app/extern/include -L/Applications/MATLAB_R2020b.app/bin/maci64 -lmat -lmx -lmex -Xpreprocessor -fopenmp -lomp -I/Users/dinula/Documents/repos/ADMM-Code/LBFGSpp/include -I/usr/local/include/eigen3/';
ipathEigen = ['-I', '/usr/local/include/eigen3/'];
ipathLBFGS = ['-I', '../LBFGSpp/include'];
iPathOpenMP = ['-I', '/usr/local/opt/libomp/include'];
lPathOpenMP = ['-L', '/usr/local/opt/libomp/lib'];
iPathMatlab = ['-I', '/Applications/MATLAB_R2020b.app/extern/include'];
lPathMatlab = ['-L', '/Applications/MATLAB_R2020b.app/bin/maci64'];
cFlags = ['CXXFLAGS=', '$CXXFLAGS -std=c++11 -lmat -lmx -lmex'];
compFlags = ['COMPFLAGS=', '$COMPFLAGS -Xpreprocessor -fopenmp -lomp'];
% ldFlags = ['LDFLAGS=', '$LDFLAGS -Xpreprocessor -fopenmp -lomp'];
ldFlags = ['LDFLAGS=', '$LDFLAGS -Xpreprocessor -fopenmp -lomp -O2'];
mex('-g','-largeArrayDims','-v',ipathEigen,iPathOpenMP,iPathMatlab,lPathMatlab,cFlags,compFlags,ldFlags,ipathLBFGS,sourceFile1, sourceFile2, sourceFile4);
fprintf('\n')

% sourceFile = ['MexTest.cpp'];
% fprintf('Building %s \n', sourceFile)
% mex('-g','-largeArrayDims',sourceFile);
% fprintf('\n')


%mex -g -largeArrayDims Source\L0_NeedleSim5DoFStep.c
%mex -g -largeArrayDims Source\L0_NeedleSim5DoFFinal.c
%mex -g -v -largeArrayDims Source\L0_NeedleSim5DoF.c
%mex -g -v -largeArrayDbeepims Source\L0_NeedleSim5DoFFinal.c
