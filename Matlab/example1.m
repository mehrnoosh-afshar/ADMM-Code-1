clear all 
% Cylandrical Model with Tet mesh
R = 0.1;
H = 0.5;
gm = multicylinder(R,H);
% Assign it to a generic PDE model and plot the geometry.
model = createpde;
model.Geometry = gm;
figure
pdegplot(model,'FaceLabels','on')
%Generate mesh and plot the mesh.
mesh=generateMesh(model, 'GeometricOrder','linear','Hmax',0.05,'Hmin',0.0011)
figure
pdemesh(model);

 
Elements=model.Mesh.Elements;
Nodes=model.Mesh.Nodes;

% Find Surface Nodes 
Nf2 = findNodes(mesh,'region','Face',2);
NBc = findNodes(mesh,'region','Face',1);

figure
pdemesh(model)
hold on
plot3(Nodes(1,Nf2),Nodes(2,Nf2),Nodes(3,Nf2),'ok','MarkerFaceColor','g') 
plot3(Nodes(1,NBc),Nodes(2,NBc),Nodes(3,NBc),'ok','MarkerFaceColor','r') 
%%
%----------------------------------------------------------------------
% Mechanical Properties 
E=3000;
rho=1000;
nu=0.49;
k = E/(3*(1-2*nu));
miu = E/(2*(1+nu));
c10=0.5*miu;
D11=2/k;
Dt=0.005;

% Make Extrnal force vector 
F_extr = zeros(size(Nodes,2)*3,1);
g = 9.8;
NodeNumber_whole = (1:1:size(Nodes,2));
direction_index = 3; % gravity is z direction
Node_globalIndex_direction = get_Node_globalIndex_dircetion(NodeNumber_whole,direction_index);
% F_extr(Node_globalIndex_direction) = F_extr(Node_globalIndex_direction) + (-M0*g); % Apply Gravity Force to all Nodes 
N = get_Node_globalIndex(Nf2);
Nf2_Node_globalIndex = get_Node_globalIndex_dircetion(Nf2,3);
F_extr(Nf2_Node_globalIndex,1) = F_extr (Nf2_Node_globalIndex,1)-10; % Apply -10 N on nodes on Face2 

% Make C matrix and d vector  of  dirichlet B.C 
x_n = reshape(Nodes,[],1);
NBc_Node_globalIndex = get_Node_globalIndex(NBc);
C_BC = sparse(size(NBc_Node_globalIndex,2),size(Nodes,2)*3);

for l=1:1:size(NBc_Node_globalIndex,2)
    C_BC(l,NBc_Node_globalIndex(l))=1;
end
d = x_n(NBc_Node_globalIndex);

%%
%--- Pre computations and mass assignments 
sizeX = 3*size(Nodes,2);
D = sparse(size(Elements,2)*9,sizeX);
W = sparse(size(Elements,2)*9,size(Elements,2)*9);
M0 = sparse(sizeX,sizeX);
[va,vi] = volume(mesh);
mm=va*rho/size(Nodes,2);
tic
for iter=1:1:size(Elements,2)
       NodeNumber = Elements(:,iter);
       Node_globalIndex = get_Node_globalIndex(NodeNumber');

       m = vi(iter)*rho;
       for jj=1:size(Node_globalIndex,2)
          %M0(Node_globalIndex(jj),Node_globalIndex(jj))= M0(Node_globalIndex(jj),Node_globalIndex(jj))+ m/4;
          M0(Node_globalIndex(jj),Node_globalIndex(jj))=mm;
       end

    Di = get_reduction_localMatrix(Node_globalIndex,sizeX);
    Binv_element = inv(reshape(Di*X_n,3,3));
    Binv_element_trans = Binv_element';
    B2 = sparse(blkdiag(Binv_element_trans,Binv_element_trans,Binv_element_trans));

%      F = reshape(Di*x_n,3,3)*Binv_element; % for test delete later
%      tt = reshape(DDi*x_n,3,3);
%      reshape(B2*DDi*x_n,3,3)
%      d=B2*DDi;
    DDi = get_reduction_localMatrix2(Node_globalIndex,sizeX);

      
    Wi = sparse(get_weight(vi(iter),k));
    
    D(1+(iter-1)*9:iter*9,:) = B2*DDi;
    W(1+(iter-1)*9:iter*9,1+(iter-1)*9:iter*9) = Wi;
    
end
toc

M_inv = inv(M0);

A = M0+Dt^2*D'*(W'*W)*D;
A_inv = inv(A);
C1 = Dt^2*D'*(W'*W);
Augm_Matrix = [A,C_BC';C_BC,sparse(size(NBc_Node_globalIndex,2),size(NBc_Node_globalIndex,2))];
Augm_inv = inv(Augm_Matrix);

wc = D'*(W'*W)*D;


%% Solver Global + Local 

% Calculate x_bar 
X_n = reshape(Nodes,[],1);
x_n = reshape(Nodes,[],1);
v_n = x_n*0;
x_bar = update_x_bar(x_n,v_n,Dt,M_inv,F_extr);
u = sparse(9,size(Elements,2));
z = sparse(9,size(Elements,2));

z_n =[];
u_n =[];
tspan = 0.05 ; % Simulation Duration 
time_number = length((0:Dt:tspan));
taw =5* k;  % weightin in local step

% ADMM initialization 
% temporal parameters are used in optimization 
%curr_x = x_bar;
curr_z = D*x_n;
curr_u = sparse(size(curr_z,1),1);
curr_sigma = (get_current_sigma(curr_z,size(Elements,2)))';


options = optimoptions(@fminunc,'MaxIterations',20,'Algorithm','quasi-newton');
options2 = optimoptions(@fminunc,'MaxIterations',20,'Algorithm','quasi-newton','SpecifyObjectiveGradien',true);

% Calculte initia sigma0 for the local step out side the loop 
% sigma0 is calculated based on the deriviation of element one

NodeNumber = Elements(:,1);
X_Node_globalIndex = get_Node_globalIndex(NodeNumber');
X = reshape(x_n(X_Node_globalIndex),3,4);
B = [X(:,2)-X(:,1),X(:,3)-X(:,1),X(:,4)-X(:,1)];
U = [1,0.5,0.3,1;0,0,0,0.1;2,5,3,0]*0.001;  % displacement of the nodes
z00 = B+[U(:,2)-U(:,1),U(:,3)-U(:,1),U(:,4)-U(:,1)];
F0 = z00*inv(B); 
[U,sigma0,V]=svd(F0);
sigma0 = diag(sigma0);

% parpool('local',6)
tic
for time_iter =1:1:200
    
        % Start from global loop 
%               s1=C1*(curr_z);
%               s2= M0*x_bar;
%               
              b = M0*x_bar + C1*(curr_z-curr_u);
%               vv = A_inv * b;
%               xu = vv(1:sizeX);

              
              Augm_b = vertcat(b,d);
              v = Augm_inv*Augm_b;
              vv = A_inv * b;
              curr_x = v(1:sizeX);
              norm(curr_x-x_n)
              norm(curr_x-x_bar)
              v_n = Dt^-1*(curr_x-x_n);
              x_n = curr_x;
              x_bar = update_x_bar(x_n,v_n,Dt,M_inv,F_extr);
              u= reshape(curr_u,9,size(Elements,2));
              velocity_norm(time_iter)= norm(v_n);
              
    % Local step in parallel computing 
          for iter=1:1:size(Elements,2)
              
              NodeNumber = Elements(:,iter);
              Node_globalIndex = get_Node_globalIndex(NodeNumber');
              Di = get_reduction_localMatrix(Node_globalIndex,sizeX);
              u_n_element = u(:,iter);
              
              
              
              Binv_element = inv(reshape(Di*X_n,3,3));
              Binv_element_trans = Binv_element';
              
              B2 = blkdiag(Binv_element_trans,Binv_element_trans,Binv_element_trans);
              
              DDi = get_reduction_localMatrix2(Node_globalIndex,sizeX);
              
                %dummy1 = Di*x_n;
                dummy11 = B2*DDi*x_n;  % this gives the vector elemnt of transpose(F)

%               dummy2 = reshape(Di*x_n,3,3)* Binv_element + ...  
%                  reshape(u_n_element,3,3); % if you want to work with F
              
              dummy2 = reshape(dummy11,3,3)+ ...
                  reshape(u_n_element,3,3); % transpose(F)+un
              
              F_bar = dummy2;
              [U,sigma_bar,V] = modifiedSVD(F_bar');
              %[U,sigma_bar,V]=svd(F_bar);
              sigma_bar = diag(sigma_bar);
              ve = vi(iter);
              myfun_3var = @(sigma) local_costfunction2_tetelemnt(sigma,sigma_bar,miu,k,taw,ve);

              % Do optimization 
              [Sigma_op,Fp2] = fmincon(myfun_3var,sigma0,[],[],[],[],zeros(3,1),+Inf);
              %[Sigma_op,Fp2]  = fminunc(myfun_3var,sigma0,options);
               
              S(:,iter)= Sigma_op;
              
              z_n_element = reshape((U*diag(Sigma_op)*V')',[],1);
              z(:,iter) = z_n_element;
              
              u_n_element_new = reshape(dummy2,[],1) - z_n_element;
              u(:,iter) = u_n_element_new;

          end
          
          curr_z = reshape(z,[],1);
          curr_u = reshape(u,[],1);
          residual(time_iter) = norm((D*x_n-curr_z));
          
          % global step update
% %               sigma0 = Sigma_op;
% %               z_n = reshape(z,[],1);
% %               u_n = reshape(u,[],1);
% %               b = M0*x_bar + C1*(z_n-u_n);
% %               Augm_b = vertcat(b,d);
% %               %vv = A_inv*b;
% %               v = Augm_inv*Augm_b;
% %               x_n0 = x_n;
% %               x_n = v(1:sizeX);
% %               v_n = Dt^-1*(x_n-x_n0);
% %               x_bar = update_x_bar(x_n,v_n,Dt,M_inv,F_extr);
% %               residual(time_iter) = norm((D*x_n-z_n));
end
toc

%  b = M0*x_bar + C1*(curr_z-curr_u);
%               Augm_b = vertcat(b,d);
%               v = Augm_inv*Augm_b;
%               curr_x = v(1:sizeX);
%               norm(curr_x-x_n)
%               norm(curr_x-x_bar)
%%  Show result on MESH 
model2 = createpde();
nodes = reshape(curr_x(:,1),3,size(Nodes,2));
V = reshape(v_n(:,1),3,size(Nodes,2));
XB = reshape(x_bar(:,1),3,size(Nodes,2));
for ii=1:1:size(Elements,2)
  NodeNumber = Elements(:,ii);
  vol(ii) = tetrahedronVolume(nodes(:,NodeNumber)');
end
xc=find(vol<0)
mesh2=geometryFromMesh(model2,nodes,Elements);
figure;
subplot(1,2,1)
pdemesh(model2);
subplot(1,2,2)
pdemesh(model);

figure 
pdemesh(model2);
hold on 
pdemesh(model);


NodeNumber = Elements(:,xc);
figure;
%pdemesh(model2)
plot3(nodes(1,:),nodes(2,:),nodes(3,:),'ok','MarkerFaceColor','r') 
hold on 
% plot3(nodes(1,NBc),nodes(2,NBc),nodes(3,NBc),'ok','MarkerFaceColor','y') 
hold on 
plot3(Nodes(1,:),Nodes(2,:),Nodes(3,:),'ok','MarkerFaceColor','g') 
hold on 
plot3(nodes(1,Nf2),nodes(2,Nf2),nodes(3,Nf2),'ok','MarkerFaceColor','b') 
hold on 
plot3(nodes(1,NodeNumber),nodes(2,NodeNumber),nodes(3,NodeNumber),'ok','MarkerFaceColor','y') 

bc1=[Nodes(1,Nf2)',Nodes(2,Nf2)',Nodes(3,Nf2)'];
bc2=[nodes(1,Nf2)',nodes(2,Nf2)',nodes(3,Nf2)'];
e=bc1-bc2;
ee=max(abs(e(3,:)));


figure;
plot3(nodes(1,Nf2),nodes(2,Nf2),nodes(3,Nf2),'ok','MarkerFaceColor','b') 
hold on 
plot3(Nodes(1,Nf2),Nodes(2,Nf2),Nodes(3,Nf2),'ok','MarkerFaceColor','g') 
hold on 
plot3(XB(1,Nf2),XB(2,Nf2),XB(3,Nf2),'ok','MarkerFaceColor','Y')
% hold on 
% quiver3(nodes(1,Nf2),nodes(2,Nf2),nodes(3,Nf2),V(1,Nf2),V(2,Nf2),V(3,Nf2))

figure
plot(nodes(1,Nf2(1:12)),nodes(2,Nf2(1:12)),'ok','MarkerFaceColor','b') 
hold on 
quiver(nodes(1,Nf2(1:12)),nodes(2,Nf2(1:12)),V(1,Nf2(1:12)),V(2,Nf2(1:12)))

sum(V(1,Nf2(1:12)))



nxu = reshape(xu(:,1),3,size(Nodes,2));
ns1 = reshape(s1(:,1),3,size(Nodes,2));
ns2 = reshape(s2(:,1),3,size(Nodes,2));
XB = reshape(x_bar(:,1),3,size(Nodes,2));
bb = reshape(b(:,1),3,size(Nodes,2));

figure
plot3(nxu(1,Nf2),nxu(2,Nf2),nxu(3,Nf2),'ok','MarkerFaceColor','b') 
hold on 
% plot3(ns1(1,Nf2),ns1(2,Nf2),ns1(3,Nf2),'ok','MarkerFaceColor','r') 
hold on 
% plot3(ns2(1,Nf2),ns2(2,Nf2),ns2(3,Nf2),'ok','MarkerFaceColor','g') 
hold on 
plot3(XB(1,Nf2),XB(2,Nf2),XB(3,Nf2),'ok','MarkerFaceColor','y') 
hold on 
%plot3(bb(1,Nf2),bb(2,Nf2),bb(3,Nf2),'ok','MarkerFaceColor','b') 

%% Test Local step Optimization 

ee = reshape(e,1,[]);


% Here is an simple example of how to use optLBFGS to solve optimization(min) problem. 
% we need to know the function value and gradient , as myfun return
% x0 = 0.5*ones(6,1);

% X(:,1)=[0,0,0]';
% X(:,2)=[0,0,1]';
% X(:,3)=[1,0,0]';
% X(:,4)=[0,1,0]';
% B1 = [X(:,2)-X(:,1),X(:,3)-X(:,1),X(:,4)-X(:,1)];
% 
% z0 = B+rand(3,3);
% z0 = reshape(z0,[],1); 


% x_curr = x_bar0;
% u_n = u0; 
% NodeNumber = Elements(:,iter);
% X_Node_globalIndex = get_Node_globalIndex(NodeNumber');
% X = reshape(x_n(X_Node_globalIndex),3,4);
% B = [X(:,2)-X(:,1),X(:,3)-X(:,1),X(:,4)-X(:,1)];
% U = [1,0.5,0.3,1;0,0,0,0.1;2,5,3,0]*0.001;
% z00 = B+[U(:,2)-U(:,1),U(:,3)-U(:,1),U(:,4)-U(:,1)];
% z0 = reshape(z00,[],1); 
% F0 = z00*inv(B); 
% F_bar = reshape(Di*x_n +u_n,3,3);
% [U,sigma_bar,V]=svd(F_bar);
% [U,sigma0,V]=svd(F0);
% sigma_bar = diag(sigma_bar)
% sigma0 = diag(sigma0)
% 
% taw = 1; 
% %z=z0;
% myfun2 = @(z) local_costfunction1_tetelemnt(X,z,miu,k,Di,Wi,x_curr,u_n);
% myfun3 = @(sigma) local_costfunction2_tetelemnt(sigma,sigma_bar,miu,k,taw);
% myfun= @(z)energyTetelement(X,z,miu,k)
% 
% [a,b]=myfun2(z0)
% [a,b]=myfun(z0)
% 
% 
% %
% maxIter = 100;
% % memsize <= length(x0)
% % % memsize : 5~20 will be fine
% memsize = 3; 
% 
% %[XM,FM,kM] =optLBFGS(myfun2,z0,maxIter,memsize)
% 
% options = optimoptions(@fminunc,'MaxIterations',10,'Algorithm','quasi-newton');
% options2 = optimoptions(@fminunc,'MaxIterations',10,'Algorithm','quasi-newton','SpecifyObjectiveGradien',true);
% 
% parpool('local',6)
% tic
% parfor iter=1:size(Elements,2)
% %[Xp,Fp]  = fminunc(myfun2,z0,options);
% [Xp2,Fp2]  = fminunc(myfun3,sigma0,options2);
% end
% toc
% 
% 
% 
% function x_bar = update_x_bar(x_n,v_n,Dt,M_inv,F_extr)
%  x_bar = x_n + v_n * Dt + M_inv * F_extr* Dt^2;
% end
% 
%%
%----- Check Generalized Force 
NodeNumber = Elements(:,1);
X_Node_globalIndex = get_Node_globalIndex(NodeNumber');
X = reshape(x_n(X_Node_globalIndex),3,4);
B = [X(:,2)-X(:,1),X(:,3)-X(:,1),X(:,4)-X(:,1)];
U = [1,0.5,0.3,1;0,0,0,0.1;2,5,3,0]*0.001;  % displacement of the nodes
z00 = B+[U(:,2)-U(:,1),U(:,3)-U(:,1),U(:,4)-U(:,1)];
F0 = z00*inv(B); 
[U,sigma0,V]=svd(F0);
sigma0 = diag(sigma0);

[Strain_Energy,Strain_Energy_gradiant0]=energyTetelement(X,reshape(z00,[],1),miu,k);

Di = get_reduction_localMatrix(X_Node_globalIndex,sizeX);
DDi = get_reduction_localMatrix2(X_Node_globalIndex,sizeX);

gi= reshape(Strain_Energy_gradiant0,[],1);
ggi= reshape(Strain_Energy_gradiant0',[],1);


a= rand(3,3);
gi= reshape(a,[],1);
ggi= reshape(a',[],1);
norm(Di'*gi-DDi'*ggi)

