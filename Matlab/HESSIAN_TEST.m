% Script to check the hessian of the strain energy 
% N = [1,0,0;
%    0,1,0;
%    0,0,1;
%   -1,-1,-1];  % this is for thetradron 

N_triangle = [1,0;0,1;-1,-1];

model = createpde;
R1 = [3,4,-1,1,1,-1,0.5,0.5,-0.75,-0.75]';
gm = [R1];
sf = 'R1';
ns = char('R1');
ns = ns';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);
 pdegplot(model,'EdgeLabels','on')
% axis equal
% xlim([-1.1,1.1])
mesh = generateMesh(model,'GeometricOrder','Linear');
pdeplot(model)

Elements = model.Mesh.Elements;
Nodes = model.Mesh.Nodes;

node = Elements(:,1);
X = Nodes(:,node);
x = X + rand(2,3);

Dx = X*N_triangle ; 
Dx_inv = inv(Dx);
DF_X = kron(transpose(N_triangle*Dx_inv),eye(2,2));
F = x*N_triangle* Dx_inv;

% Knn = zeros(16,16);
% Knn(1,1) = 1;
% Knn(2,5) = 1;
% Knn(3,9) = 1;
% Knn(4,13) = 1;
% 
% Knn(5,2) = 1;
% Knn(6,6) = 1;
% Knn(7,10) = 1;
% Knn(8,14) = 1;
% 
% Knn(9,3) = 1;
% Knn(10,7) = 1;
% Knn(11,11) = 1;
% Knn(12,15) = 1;
% 
% Knn(13,4) = 1;
% Knn(14,8) = 1;
% Knn(15,12) = 1;
% Knn(16,16) = 1;


Knn = [1,0,0,0;0,0,1,0;
    0,1,0,0;0,0,0,1];
i =trace(F'*F);
j = det(F);

E=3000;
rho=1000;
nu=0.49;
k = E/(3*(1-2*nu));
miu = E/(2*(1+nu));
c10=0.5*miu;
D11=2/k;

inv_F= inv(F);
vec_F = reshape(F,[],1);
vec_INVF_T = reshape(inv_F',[],1);

DS_F = ((miu/3)*(j^(-2/3))*i-k*j*(j-1))*Knn'*kron(inv_F',inv_F)+...
        (miu*(j^(-2/3))*kron(eye(2,2),eye(2,2)))+...
        (miu*(2/9)*(j^(-2/3))*i+k*j*(2*j-1))* (vec_INVF_T*vec_INVF_T')+...
        (-4/3*miu*(j^(-2/3)))* (vec_INVF_T*vec_F');

Hessian = DF_X'*DS_F*DF_X


% Calculate Analytically 
syms x1_1 x1_2  x2_1 x2_2  x3_1 x3_2 

xa =[x1_1,x2_1,x3_1;...
     x1_2,x2_2,x3_2];
 
Say = strain_energy_density_analytic (xa,N_triangle,Dx_inv,k,miu)

%d1= diff(Say,x3_2,2);
d1= diff(diff(Say,x2_2,1),x3_2,1);


s=subs(subs(subs(subs(d1,x1_1,x(1,1)),x1_2,x(2,1)),x2_1,x(1,2)),x2_2,x(2,2));
s= subs(subs(s,x3_1,x(1,3)),x3_2,x(2,3));
s= eval(s)
