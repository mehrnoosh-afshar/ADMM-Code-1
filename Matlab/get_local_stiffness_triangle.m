function K_local = get_local_stiffness_triangle (Current_Node_Value, Initial_Node_Value,E,nu)
% This equations are developed for Hyper elastic Material 
N_triangle = [1,0;0,1;-1,-1]; % This is for triangle mesh 

Knn = [1,0,0,0;0,0,1,0;
    0,1,0,0;0,0,0,1]; % This is just for triangle Mesh 

X = Initial_Node_Value;
x = Current_Node_Value; 

Dx = X*N_triangle ; 
Dx_inv = inv(Dx);
DF_X = kron(transpose(N_triangle*Dx_inv),eye(2,2));
F = x*N_triangle* Dx_inv;
i =trace(F'*F);
j = det(F);

k = E/(3*(1-2*nu));
miu = E/(2*(1+nu));

inv_F= inv(F);
vec_F = reshape(F,[],1);
vec_INVF_T = reshape(inv_F',[],1);

DS_F = ((miu/3)*(j^(-2/3))*i-k*j*(j-1))*Knn'*kron(inv_F',inv_F)+...
        (miu*(j^(-2/3))*kron(eye(2,2),eye(2,2)))+...
        (miu*(2/9)*(j^(-2/3))*i+k*j*(2*j-1))* (vec_INVF_T*vec_INVF_T')+...
        (-4/3*miu*(j^(-2/3)))* (vec_INVF_T*vec_F');

Hessian = DF_X'*DS_F*DF_X;

K_local = Hessian;


end

