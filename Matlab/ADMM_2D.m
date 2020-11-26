function [Current_configuration,S,S0] = ADMM_2D(x_desire,NumerLoop,Elements,Nodes,Ai,At,NBC,Nf1,Nf2,NP)
% This function get the information of Mesh and Boundary points and give
% the final position of nodes 
% x_desire : it is the desire location of boundary points 
% NumberLoop is the number of golobal step (Number of the outter loop if you see the code sudo code in the paper)
% Elements is the matrix containg the node numbers associated with each
% element 
% Nodes is the coordinate of nodes in cartesian coordinate 
% Ai is the vector of area of all elements 
% At is the total area of the mesh 
% NBC is nodes associated with fixed boundary condition 
% Nf1 is nodes associated with control points on top of the plate 
% Nf2 is nodes associated with control points on bottom of the plate 
% Np is nodes associated with inside target points 

% Mechanical properties of object 
E=3000;
rho=1000;
nu=0.49;
k = E/(3*(1-2*nu));
miu = E/(2*(1+nu));
%c10=0.5*miu;
%D11=2/k;
%Dt=0.005;
Dt=0.001;  % Dt is sampling time 

%% 
%Make Extrnal force vector 
% This section is necessaary when you eant to apply force on boundary
% condition not a prescribe position 

% % F_extr = zeros(size(Nodes,2)*2,1);
NodeNumber_whole = (1:1:size(Nodes,2));
% % 
Nf1_Node_globalIndex = get_Node_globalIndex_dircetion2D(Nf1,2);
Nf2_Node_globalIndex = get_Node_globalIndex_dircetion2D(Nf2,2);
% % 
% % 
% % F_extr(Nf1_Node_globalIndex,1) = F_extr (Nf1_Node_globalIndex,1)-10; % Apply -10 N on nodes on Face2 
% % F_extr(Nf2_Node_globalIndex,1) = F_extr (Nf2_Node_globalIndex,1)+10; % Apply -10 N on nodes on Face2 


%%
% Make C matrix and D vector  of  dirichlet B.C 
x_n = reshape(Nodes,[],1);
% each coordinate of a node is DOF. DOF= number of node *(2 or 3)
% Working with vector of DOF , nodes index should be transformed to DOF
% index
NBc_Node_globalIndex = get_Node_globalIndex2D(NBC); % Get the index associated with the nodes in DOF vector 
C_BC = sparse(size(NBc_Node_globalIndex,2),size(Nodes,2)*2);

for l=1:1:size(NBc_Node_globalIndex,2)
    C_BC(l,NBc_Node_globalIndex(l))=1;
end
d = x_n(NBc_Node_globalIndex);
% --- Pre computations and mass assignments ------------------------------
% D and W matrix are global materices used in equation 18 in paper 
sizeX = 2*size(Nodes,2);
D = sparse(size(Elements,2)*4,sizeX); % Multiplied with 4 beacuse deformation matrix 
                                       % associated with 2D problem is a
                                       % 2by2 matrix 
                                    
W = sparse(size(Elements,2)*4,size(Elements,2)*4);
M0 = sparse(sizeX,sizeX);
%[At,Ai] = area(mesh);
mm=At*rho/size(Nodes,2);
X_n = reshape(Nodes,[],1);

tic
for iter=1:1:size(Elements,2)
       NodeNumber = Elements(:,iter);
       Node_globalIndex = get_Node_globalIndex2D(NodeNumber');

       m = Ai(iter)*rho;
       for jj=1:size(Node_globalIndex,2)
          M0(Node_globalIndex(jj),Node_globalIndex(jj))= M0(Node_globalIndex(jj),Node_globalIndex(jj))+ m/4;
          %M0(Node_globalIndex(jj),Node_globalIndex(jj))=mm;
       end

    DDi = get_reduction_localMatrix2_2D(Node_globalIndex,sizeX);

    Binv_element_trans = inv(reshape(DDi*X_n,2,2));
    B2 = sparse(blkdiag(Binv_element_trans,Binv_element_trans));
    
    Wi = sparse(get_weight_2D(Ai(iter),k)); % To change the weighting matrix go to "get_weight_2D" function 
                                            % and change the weight
    
    D(1+(iter-1)*4:iter*4,:) = B2*DDi;
    W(1+(iter-1)*4:iter*4,1+(iter-1)*4:iter*4) = Wi;
    
end
toc

M_inv = inv(M0);

% In this section matixes in equation 18 are calculated
% They are fixed matrices though simulation
% To apply the fixed boundary condition augmented lagranigain form of
% equation 18 is solved which is "Augm_inv" matrix 
A = M0+Dt^2*D'*(W'*W)*D;
A_inv = inv(A);
C1 = Dt^2*D'*(W'*W);
Augm_Matrix = [A,C_BC';C_BC,sparse(size(NBc_Node_globalIndex,2),size(NBc_Node_globalIndex,2))];
Augm_inv = inv(Augm_Matrix);

wc = D'*(W'*W)*D; 

%%
% ------------- Solver Global + Local--------------------------
% Calculate x_bar 
X_n = reshape(Nodes,[],1);  % in the first values of all DOF 
x_n = reshape(Nodes,[],1);  % x_n is getting updated in the loop 
                            % it is the current position of DOF 
v_n = x_n*0;                % Initilize the velocity 

% If you are applying force on control points for x_bar you need to use
% this
%x_bar = update_x_bar(x_n,v_n,Dt,M_inv,F_extr);

% If you want do position control use this to update and initilize x_bar  
Index_control =[Nf1_Node_globalIndex , Nf2_Node_globalIndex ];
%%dx = rand(size(Nf1_Node_globalIndex,2)+size(Nf2_Node_globalIndex,2),1)*0.01;
x_bar = update_x_bar_position(x_n,v_n,Dt,x_desire,Index_control);

u = sparse(4,size(Elements,2));
z = sparse(4,size(Elements,2));

% taw =k;  % weightin in local step

% ADMM initialization 
% temporal parameters are used in optimization 
% curr_x and curr_u are vectors used in the loop 
curr_x = x_bar;    % Initilized eaqual to x_bar 
%curr_z = D*x_n;
curr_u = sparse(4,1);  %Initilized to Zero


options = optimoptions(@fminunc,'MaxIterations',20,'Algorithm','quasi-newton');
options2 = optimoptions(@fminunc,'MaxIterations',20,'Algorithm','quasi-newton','SpecifyObjectiveGradien',true);

% Note: Zi vector can be introduced to be the element of deformation
% gradient matrix F that has 4 element in case of 2D problem and 9 element
% in case of 3D.
% Beacuse we consider our material to be isotropic. The invarient of matrix
% F can be used as Zi vetor and they depends on singular values of matrix F
% So insted of a optimization with 4 , we can have an optimization with 2
% element in 2D case and insted of a optimization with 9 , we can have an 
% optimization with 3 element in 3D case



% Calculte initia sigma0 for the local step out side the loop 
% sigma0 is calculated based on the deriviation of element one
% Sigma0 is the initial value of vector sigma which are optimized in the
% loop 

NodeNumber = Elements(:,1);
X_Node_globalIndex = get_Node_globalIndex2D(NodeNumber');
X = reshape(x_n(X_Node_globalIndex),2,3);
B = [X(:,2)-X(:,1),X(:,3)-X(:,1)];
U = [1,0.5,0.3;0,0,0]*0.001;  % displacement of the nodes
z00 = B+[U(:,2)-U(:,1),U(:,3)-U(:,1)]; % Add an arbitrary displacement to 
                                       % current location of each nodes
                                       % element
F0 = z00*inv(B); 
[U,sigma0,V]=svd(F0);
sigma0 = diag(sigma0);                % Calculate the initial singular value 



% parpool('local',6)
for time_iter =1:1:NumerLoop % The global step 
 tic   
             if (time_iter >= 2)
             curr_z0 = curr_z;
             end
    
    % Local step in parallel computing 
              for iter=1:1:size(Elements,2)  % This is the local loop over the number of 
                                         % loops that needs to be accerelated
              
              NodeNumber = Elements(:,iter);
              Node_globalIndex = get_Node_globalIndex2D(NodeNumber');
              DDi = get_reduction_localMatrix2_2D(Node_globalIndex,sizeX);
              u_n_element = u(:,iter); % get ui (eq 27 and 28)
              
              Binv_element_trans = inv(reshape(DDi*X_n,2,2));
              B2 = blkdiag(Binv_element_trans,Binv_element_trans);
              
              
             %dummy1 = Di*x_n;
              dummy11 = B2*DDi*curr_x;  % this gives the vector elemnt of transpose(F) for 
                                        %  each element that is your zi

%               dummy2 = reshape(Di*x_n,3,3)* Binv_element + ...  
%               reshape(u_n_element,3,3); % if you want to work with F
%               Not necessary 
              
              dummy2 = reshape(dummy11,2,2)+ ...
                  reshape(u_n_element,2,2); % transpose(F)+ui == zi+ui
              
              % Compute the sigualr value of current F 
              F_bar = dummy2;
              % To avoid the element inversion, the sigular value are
              % modiifed if the determent of F < 0. 
              % This is done by "modifiedSVD_2D" function 
              [U,sigma_bar,V] = modifiedSVD_2D(F_bar'); 
              %[U,sigma_bar,V]=svd(F_bar);
              if sigma_bar(2,2)<0
                  sigma_bar(2,2)=-sigma_bar(2,2);
              end
               if sigma_bar(1,1)<0
                  sigma_bar(1,1)=-sigma_bar(1,1);
              end
                  sigma_bar = diag(sigma_bar);
              S0(:,iter,time_iter) = sigma_bar;
              ve = Ai(iter);
              a=(get_weight_2D(Ai(iter),k));
              taw = a(1,1)^2;
%             taw = k;
              
     
              myfun_2var = @(sigma) local_costfunction2_triangle(sigma,sigma_bar,miu,k,taw,ve);

              % Do optimization 
              % Optimize the cost funcrion based on eq 27 
              % Actually just optimizing the singular values and the whole
              % matrix F 
              [Sigma_op,Fp2] = fmincon(myfun_2var,sigma0,[],[],[],[],zeros(2,1),+Inf);
              %[Sigma_op,Fp2]  = fminunc(myfun_3var,sigma0,options);
               
              S(:,iter)= Sigma_op;
              
              % After finding the optimized value of singular values 
              % Update the value of Zi 
              z_n_element = reshape((U*diag(Sigma_op)*V')',[],1);
              % Update the global Z
              z(:,iter) = z_n_element;
              
              % update ui (eq 28)
              u_n_element_new = reshape(dummy2,[],1) - z_n_element;
              % update the global U 
              u(:,iter) = u_n_element_new;

              end
          % update curr_z and curr_u to be used in global step
          curr_z = reshape(z,[],1);
          curr_u = reshape(u,[],1);

% Global Step 
       
              % x_bar = update_x_bar(x_n,v_n,Dt,M_inv,F_extr);
             

              b = M0*x_bar + C1*(curr_z-curr_u);
              
              Augm_b = vertcat(b,d);
              v = Augm_inv*Augm_b;
              curr_x0= curr_x;
              curr_x = v(1:sizeX);

              v_n = Dt^-1*(curr_x-curr_x0);
              velocity_norm(time_iter)= norm(v_n);
              
              x_bar = update_x_bar_position(curr_x,v_n,Dt,x_desire,Index_control);
              x_bar(Nf1_Node_globalIndex);
              
             residual0(time_iter) = norm(W*(D*curr_x-curr_z));
             
             if (time_iter >= 2)
             residual1(time_iter) = norm(D'*W'*W*(curr_z-curr_z0));
             end
end

Current_configuration = curr_x;
end