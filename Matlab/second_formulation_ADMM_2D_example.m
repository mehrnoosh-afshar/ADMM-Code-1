%%%%%%%%%% Second Version of Code revise convergence 
clear all
%--------------- First make a 2D mesh ---------------------------------
model = createpde;
R1 = [3,4,0,0.2,0.2,0,0.05,0.05,-0.05,-0.05]';
% R1 = [3,4,0,0.2,0.2,0,0.1,0.1,0,0]';
gm = [R1];
sf = 'R1';
ns = char('R1');
ns = ns';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);
pdegplot(model,'EdgeLabels','on')
% axis equal
% xlim([-1.1,1.1])
% mesh = generateMesh(model,'GeometricOrder','Linear','Hmax',0.04);
mesh = generateMesh(model,'GeometricOrder','Linear');
pdeplot(model)


Elements = model.Mesh.Elements;
Nodes = model.Mesh.Nodes;

% Find Surface Nodes 
NBC = findNodes(mesh,'region','Edge',4);
Nf1 = findNodes(mesh,'box',[0.05 0.1],[0.05 0.05]);
Nf2 = findNodes(mesh,'box',[0.05 0.1],[-0.05 -0.05]);
NP = findNodes(mesh,'radius',[0.15 0],.004);
NM =[Nf1,Nf2];

figure
pdemesh(model)
hold on
plot(Nodes(1,NBC),Nodes(2,NBC),'ok','MarkerFaceColor','g') 
plot(Nodes(1,Nf1),Nodes(2,Nf1),'ok','MarkerFaceColor','r') 
plot(Nodes(1,Nf2),Nodes(2,Nf2),'ok','MarkerFaceColor','r') 

plot(Nodes(1,NP),Nodes(2,NP),'ok','MarkerFaceColor','b') 

%set(gca,'XTick',[], 'YTick', [])
%%
%%%%%%%%%% Second Version of Code revise convergence 

% 2_D example for positioning control by Jacobian calculation 
% based on Quasi-Static equations 
% clear all 
% XX=zeros(16,16);
% xx =zeros(16,16);
%--------------- First make a 2D mesh ---------------------------------
% model = createpde;
% R1 = [3,4,-1,1,1,-1,0.5,0.5,-0.75,-0.75]';
% gm = [R1];
% sf = 'R1';
% ns = char('R1');
% ns = ns';
% g = decsg(gm,sf,ns);
% geometryFromEdges(model,g);
% pdegplot(model,'EdgeLabels','on')
% % axis equal
% % xlim([-1.1,1.1])
% % mesh = generateMesh(model,'GeometricOrder','Linear','Hmax',0.04);
% mesh = generateMesh(model,'GeometricOrder','Linear');
% pdeplot(model)
% 
% 
% Elements = model.Mesh.Elements;
% Nodes = model.Mesh.Nodes;
% 
% % Find Surface Nodes 
% NBC = findNodes(mesh,'region','Edge',4);
% Nf1 = findNodes(mesh,'box',[-0.7 -0.1],[0.5 0.5]);
% Nf2 = findNodes(mesh,'box',[-0.7 -0.1],[-0.75 -0.75]);
% NP = findNodes(mesh,'radius',[0 -0.1],.05);
% NM =[Nf1,Nf2];
% 
% figure
% pdemesh(model)
% hold on
% plot(Nodes(1,NBC),Nodes(2,NBC),'ok','MarkerFaceColor','g') 
% plot(Nodes(1,Nf1),Nodes(2,Nf1),'ok','MarkerFaceColor','r') 
% plot(Nodes(1,Nf2),Nodes(2,Nf2),'ok','MarkerFaceColor','r') 
% 
% plot(Nodes(1,NP),Nodes(2,NP),'ok','MarkerFaceColor','b') 

%%
%--------------------- Mechanical Properties -----------------------------
E=3000;
rho=1000;
nu=0.49;
k = E/(3*(1-2*nu));
miu = E/(2*(1+nu));
c10=0.5*miu;
D11=2/k;
Dt=0.0045;

% Make Extrnal force vector 
F_extr = zeros(size(Nodes,2)*2,1);
NodeNumber_whole = (1:1:size(Nodes,2));

Nf1_Node_globalIndex = get_Node_globalIndex_dircetion2D(Nf1,2);
Nf2_Node_globalIndex = get_Node_globalIndex_dircetion2D(Nf2,2);


F_extr(Nf1_Node_globalIndex,1) = F_extr (Nf1_Node_globalIndex,1)-100; % Apply -10 N on nodes on Face2 
F_extr(Nf2_Node_globalIndex,1) = F_extr (Nf2_Node_globalIndex,1)+100; % Apply -10 N on nodes on Face2 



% Make C matrix and d vector  of  dirichlet B.C 
x_n = reshape(Nodes,[],1);
NBc_Node_globalIndex = get_Node_globalIndex2D(NBC);
C_BC = sparse(size(NBc_Node_globalIndex,2),size(Nodes,2)*2);

for l=1:1:size(NBc_Node_globalIndex,2)
    C_BC(l,NBc_Node_globalIndex(l))=1;
end
d = x_n(NBc_Node_globalIndex);

%% --- Pre computations and mass assignments ------------------------------
sizeX = 2*size(Nodes,2);
D = sparse(size(Elements,2)*4,sizeX);
W = sparse(size(Elements,2)*4,size(Elements,2)*4);
M0 = sparse(sizeX,sizeX);
[At,Ai] = area(mesh);
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
    
    Wi = sparse(get_weight_2D(Ai(iter),k));
    
    D(1+(iter-1)*4:iter*4,:) = B2*DDi;
    W(1+(iter-1)*4:iter*4,1+(iter-1)*4:iter*4) = Wi;
    
end
toc

M_inv = inv(M0);

A = M0+Dt^2*D'*(W'*W)*D;
A_inv = inv(A);
C1 = Dt^2*D'*(W'*W);
Augm_Matrix = [A,C_BC';C_BC,sparse(size(NBc_Node_globalIndex,2),size(NBc_Node_globalIndex,2))];
Augm_inv = inv(Augm_Matrix);

wc = D'*(W'*W)*D;


%% ------------- Solver Global + Local--------------------------

% Calculate x_bar 
X_n = reshape(Nodes,[],1);
x_n = reshape(Nodes,[],1);
v_n = x_n*0;
%x_bar = update_x_bar(x_n,v_n,Dt,M_inv,F_extr);

% If you want do position control 
Index_control =[Nf1_Node_globalIndex , Nf2_Node_globalIndex ];
dx = rand(size(Nf1_Node_globalIndex,2)+size(Nf2_Node_globalIndex,2),1)*0.01;
x_desire = [(0.031-0.000)*ones(size(Nf1_Node_globalIndex,2),1);(-0.031+0.000)*ones(size(Nf2_Node_globalIndex,2),1)];
x_bar = update_x_bar_position(x_n,v_n,Dt,x_desire,Index_control);

u = sparse(4,size(Elements,2));
z = sparse(4,size(Elements,2));

% taw =k;  % weightin in local step

% ADMM initialization 
% temporal parameters are used in optimization 
curr_x = x_bar;
%curr_z = D*x_n;
curr_u = sparse(4,1);


options = optimoptions(@fminunc,'MaxIterations',20,'Algorithm','quasi-newton');
options2 = optimoptions(@fminunc,'MaxIterations',20,'Algorithm','quasi-newton','SpecifyObjectiveGradien',true);

% Calculte initia sigma0 for the local step out side the loop 
% sigma0 is calculated based on the deriviation of element one

NodeNumber = Elements(:,1);
X_Node_globalIndex = get_Node_globalIndex2D(NodeNumber');
X = reshape(x_n(X_Node_globalIndex),2,3);
B = [X(:,2)-X(:,1),X(:,3)-X(:,1)];
U = [1,0.5,0.3;0,0,0]*0.001;  % displacement of the nodes
z00 = B+[U(:,2)-U(:,1),U(:,3)-U(:,1)];
F0 = z00*inv(B); 
[U,sigma0,V]=svd(F0);
sigma0 = diag(sigma0);

% parpool('local',6)
tic
for time_iter =1:1:500
 tic   
             if (time_iter >= 2)
             curr_z0 = curr_z;
             end
    
    % Local step in parallel computing 
          for iter=1:1:size(Elements,2)
              
              NodeNumber = Elements(:,iter);
              Node_globalIndex = get_Node_globalIndex2D(NodeNumber');
              DDi = get_reduction_localMatrix2_2D(Node_globalIndex,sizeX);
              u_n_element = u(:,iter);
              
              
              
              Binv_element_trans = inv(reshape(DDi*X_n,2,2));
              
              B2 = blkdiag(Binv_element_trans,Binv_element_trans);
              
              
                %dummy1 = Di*x_n;
                dummy11 = B2*DDi*curr_x;  % this gives the vector elemnt of transpose(F)

%               dummy2 = reshape(Di*x_n,3,3)* Binv_element + ...  
%                  reshape(u_n_element,3,3); % if you want to work with F
              
              dummy2 = reshape(dummy11,2,2)+ ...
                  reshape(u_n_element,2,2); % transpose(F)+un
              
              F_bar = dummy2;
              [U,sigma_bar,V] = modifiedSVD_2D(F_bar');
              %[U,sigma_bar,V]=svd(F_bar);
              sigma_bar = diag(sigma_bar);
              ve = Ai(iter);
               a=(get_weight_2D(Ai(iter),k));
               taw = a(1,1)^2;
     %         taw = k;
              myfun_2var = @(sigma) local_costfunction2_triangle(sigma,sigma_bar,miu,k,taw,ve);

              % Do optimization 
             
              [Sigma_op,Fp2] = fmincon(myfun_2var,sigma0,[],[],[],[],zeros(2,1),+Inf);
              %[Sigma_op,Fp2]  = fminunc(myfun_3var,sigma0,options);
               
              S(:,iter)= Sigma_op;
              
              z_n_element = reshape((U*diag(Sigma_op)*V')',[],1);
              z(:,iter) = z_n_element;
              
              u_n_element_new = reshape(dummy2,[],1) - z_n_element;
              u(:,iter) = u_n_element_new;

          end
          curr_z = reshape(z,[],1);
          curr_u = reshape(u,[],1);

% Global Step 
%        
              % x_bar = update_x_bar(x_n,v_n,Dt,M_inv,F_extr);
             

              b = M0*x_bar + C1*(curr_z-curr_u);
              
              Augm_b = vertcat(b,d);
              v = Augm_inv*Augm_b;
              curr_x0= curr_x;
              curr_x = v(1:sizeX);

              v_n = Dt^-1*(curr_x-curr_x0);
              velocity_norm(time_iter)= norm(v_n);
              
              x_bar = update_x_bar_position(curr_x,v_n,Dt,x_desire,Index_control);


             residual0(time_iter) = norm(W*(D*curr_x-curr_z));
             
             if (time_iter >= 2)
             residual1(time_iter) = norm(D'*W'*W*(curr_z-curr_z0));
             end
             xhistory (:,time_iter)= curr_x;
toc
end


toc
NP_Node_globalIndex = get_Node_globalIndex2D(NP);
%DX = curr_x(NP_Node_globalIndex,:)-X_True;
X_True2= curr_x(NP_Node_globalIndex,:);
% XXX=curr_x;
NP_Node_globalIndex = get_Node_globalIndex2D(NP);

fa=D'*(W'*W)*curr_u;
fz=fa(NP_Node_globalIndex);
all_index=(1:1:size(Nodes,2)*2);
norm(fa(NP_Node_globalIndex))
plot(residual0)
figure
plot(residual1)



xlabel('Number of iterations','FontSize',12,'FontWeight','bold')
ylabel('Sum of square errors (m)','FontSize',12,'FontWeight','bold')
p.LineWidth = 1.5;
set(gca,'fontweight','bold','fontsize',10)
%%  Show result on MESH 
model2 = createpde();
nodes = reshape(curr_x(:,1),2,size(Nodes,2));
figure;
plot(nodes(1,1:5),nodes(2,1:5),'ok','MarkerFaceColor','r') 
hold on 
plot(Nodes(1,:),Nodes(2,:),'ok','MarkerFaceColor','g') 
hold on 
plot(nodes(1,Nf1),nodes(2,Nf1),'ok','MarkerFaceColor','b') 


mesh2=geometryFromMesh(model2,nodes,Elements);
figure;
subplot(1,2,1)
pdemesh(model2);


%%
%LOOP_list = [50 , 100, 150, 300, 500];
LOOP_list = 1:6:500;

for kk=1:1:length(LOOP_list)
    LOOP = LOOP_list(kk);
UT = xhistory(:,LOOP) - x_n;
NT = 1:1:size(Nodes,2);
NT1_Node_globalIndex = get_Node_globalIndex_dircetion2D(NT,1);
NT2_Node_globalIndex = get_Node_globalIndex_dircetion2D(NT,2);
U1 = UT(NT1_Node_globalIndex);
U2 = UT(NT2_Node_globalIndex);
nodes0 = reshape(xhistory(:,LOOP),2,size(Nodes,2));

for ii=1:1:length(Elements)
    Node_index = Elements(:,ii);
    xn(:,ii) = nodes0(1,Node_index);
    yn(:,ii) = nodes0(2,Node_index)+0.05;
    un1(:,ii) = U1(Node_index);
    un2(:,ii) = U2(Node_index);
end

un_iter1(:,:,kk)=un1;
un_iter2(:,:,kk)=un2;
U1_iter(:,kk) = U1;
U2_iter(:,kk) = U2;

end
figure 
fill(xn,yn,abs(un1))
colormap jet
colorbar
xlim([0 0.22])
ylim([-0.005 0.105])
xlabel('X direction(m)','FontSize',12,'FontWeight','bold')
ylabel('Y direction(m)','FontSize',12,'FontWeight','bold')
colormap jet
c = colorbar;
c.Label.String = 'Node displacement in x directin U1 (m)';
c.FontSize =12;
c.FontWeight = 'bold'
c.Ruler.Exponent = -3;
%patch(xn,yn,'red')

% hold on 
% plot(Nodes(1,NP),Nodes(2,NP),'ok','MarkerFaceColor','b') 
% hold on 
% plot(nodes(1,NP),nodes(2,NP),'ok','MarkerFaceColor','r') 
% plot(nodes(1,Nf1),nodes(2,Nf1),'ok','MarkerFaceColor','r') 
% plot(nodes(1,Nf2),nodes(2,Nf2),'ok','MarkerFaceColor','r') 
% 
% 
% subplot(1,2,2)
% pdemesh(model);

%%
% Clorfull diagram 
% UT = xhistory(:,150) - x_n;
% NT = 1:1:size(Nodes,2);
% NT1_Node_globalIndex = get_Node_globalIndex_dircetion2D(NT,1);
% NT2_Node_globalIndex = get_Node_globalIndex_dircetion2D(NT,2);
% U1 = UT(NT1_Node_globalIndex);
% 
% U2 = UT(NT2_Node_globalIndex);
% x = curr_x(NT1_Node_globalIndex);
% y = curr_x(NT2_Node_globalIndex)+0.05;
% xv = min(x):0.001:max(x);                         % ‘x’ Vector For Interpolation
% yv = min(y):0.001:max(y);                         % ‘y’ Vector For Interpolation
% [X,Y] = meshgrid(xv,yv);                                      % Create Interpolation Grids
% Z = griddata(x, y, U1', X, Y);  
% figure
% contourf(X,Y,Z,15)
% colormap jet
% colorbar


%%
% filename1 = 'C:\Users\Mehrnoosh\Desktop\total-plate3.xlsx';
filename1 = 'C:\Users\Mehrnoosh\Desktop\total-plate3_u1u2.xlsx';

Data = readtable(filename1);
U01 =table2array(Data(:,5))*10^-3;
U02 =table2array(Data(:,6))*10^-3;
x0 = table2array(Data(:,2))*10^-3;
y0 = table2array(Data(:,3))*10^-3;
% index = 1:5:length(x0);
% xv0 = x0(index);                     
% yv0 = y0(index);                        
% [X0,Y0] = meshgrid(xv0,yv0);                                      
% Z0 = griddata(x0, y0, U01', X0, Y0);


% figure
% contourf(X0,Y0,Z0,20,'LineStyle','none')
% colormap jet
% colorbar


U01_2 = griddata(x0, y0, U01',nodes0(1,:),nodes0(2,:)+0.05,'nearest');
U02_2 = griddata(x0, y0, U02',nodes0(1,:),nodes0(2,:)+0.05,'nearest');

% shp = alphaShape(x0,y0,0.003);
% plot(shp)
% tri = alphaTriangulation(shp);
% triplot(tri,shp.Points(:,1),shp.Points(:,2))
% 
% for ii=1:1:length(tri)
%     Node_index = tri(ii,:);
%     xn(:,ii) = x0(Node_index);
%     yn(:,ii) = y0(Node_index);
%     un(:,ii) = U01(Node_index);
% end
% fill(xn,yn,abs(un))
% colormap jet
% colorbar
% xlim([0 0.22])
% ylim([-0.005 0.105])
% xlabel('X direction(m)','FontSize',12,'FontWeight','bold')
% ylabel('Y direction(m)','FontSize',12,'FontWeight','bold')
% colormap jet
% c = colorbar;
% c.Label.String = 'Node displacement in x directin U1 (m)';
% c.FontSize =12;
% c.FontWeight = 'bold'
for ii=1:1:length(Elements)
    Node_index = Elements(:,ii);
    xnn(:,ii) = nodes0(1,Node_index);
    ynn(:,ii) = nodes0(2,Node_index)+0.05;
    unn1(:,ii) = U01_2(Node_index);
    unn2(:,ii) = U02_2(Node_index);

end

figure 
fill(xnn,ynn,abs(unn1))
xlim([0 0.22])
ylim([-0.005 0.105])
xlabel('X direction(m)')
ylabel('Y direction(m)')

colormap jet
c = colorbar;
c.Label.String = 'Error between Abaques and proposed method for U1 (m)';
title ( 'Number of iteration = 50')

%% Multi figure 
figure 
fill(xnn,ynn,sqrt((unn1-un_iter1(:,:,1)).^2+(unn2-un_iter2(:,:,1)).^2))
xlim([0 0.22])
ylim([-0.005 0.105])
xlabel('X direction(m)','FontSize',12,'FontWeight','bold')
ylabel('Y direction(m)','FontSize',12,'FontWeight','bold')
colormap jet
c = colorbar;
c.Label.String = 'Displacement error field (m)';
title ( 'Number of iteration = 50')
c.FontSize =12;
c.FontWeight = 'bold'






figure 
fill(xnn,ynn,sqrt((unn1-un_iter1(:,:,2)).^2+(unn2-un_iter2(:,:,2)).^2))
xlim([0 0.22])
ylim([-0.005 0.105])
xlabel('X direction(m)','FontSize',12,'FontWeight','bold')
ylabel('Y direction(m)','FontSize',12,'FontWeight','bold')
colormap jet
c = colorbar;
c.Label.String = 'Displacement error field (m)';
title ( 'Number of iteration = 100')
c.FontSize =12;
c.FontWeight = 'bold'
figure 
fill(xnn,ynn,sqrt((unn1-un_iter1(:,:,3)).^2+(unn2-un_iter2(:,:,3)).^2))
xlim([0 0.22])
ylim([-0.005 0.105])
xlabel('X direction(m)','FontSize',12,'FontWeight','bold')
ylabel('Y direction(m)','FontSize',12,'FontWeight','bold')
colormap jet
c = colorbar;
c.Label.String = 'Displacement error field (m)';
title ( 'Number of iteration = 150')
c.FontSize =12;
c.FontWeight = 'bold'
figure 
fill(xnn,ynn,sqrt((unn1-un_iter1(:,:,4)).^2+(unn2-un_iter2(:,:,4)).^2))
xlim([0 0.22])
ylim([-0.005 0.105])
xlabel('X direction(m)','FontSize',12,'FontWeight','bold')
ylabel('Y direction(m)','FontSize',12,'FontWeight','bold')
colormap jet
c = colorbar;
c.Label.String = 'Displacement error field (m)';
title ( 'Number of iteration = 300')
c.FontSize =12;
c.FontWeight = 'bold'
figure 
fill(xnn,ynn,sqrt((unn1-un_iter1(:,:,5)).^2+(unn2-un_iter2(:,:,5)).^2))
xlim([0 0.22])
ylim([-0.005 0.105])
xlabel('X direction(m)','FontSize',12,'FontWeight','bold')
ylabel('Y direction(m)','FontSize',12,'FontWeight','bold')
colormap jet
c = colorbar;
c.Label.String = 'Displacement error field (m)';
title ( 'Number of iteration = 500')
c.FontSize =12;
c.FontWeight = 'bold'

%% Accuracy curve 

for kk=1:length(LOOP_list)
    ee(kk)= (mean(sqrt(U01_2'-U1_iter(:,kk)).^2+(U02_2'-U2_iter(:,kk)).^2));
end
p=plot(LOOP_list,ee*1000,'k')
xlabel('Number of iterations','FontSize',12,'FontWeight','bold')
ylabel('Mean of square errors (mm)','FontSize',12,'FontWeight','bold')
p.LineWidth = 1.5;
set(gca,'fontweight','bold','fontsize',10)
set(gcf,'paperunits','centimeters')
set(gcf,'position',[50 50 400 300])
%% Comparison result 

Z01 = griddata(x0, y0, U01', X, Y);
figure 
contourf(X,Y,abs(Z-Z01),15)
colormap jet
colorbar
%% Jacobian Estimation based On Quasi-Static equation 
Stiffnes = sparse (size(Nodes,2)*2,size(Nodes,2)*2);

% Build Global Stiffness Matrix form Local stifnesses 
N_triangle = [1,0;0,1;-1,-1]; % This is for triangle mesh
kk = E/(3*(1-2*nu));
miu = E/(2*(1+nu));
tic
for iter=1:1:size(Elements,2)
       NodeNumber = Elements(:,iter);
       Node_globalIndex = get_Node_globalIndex2D(NodeNumber');

       Initial_Node_Value = Nodes(:,NodeNumber);
       Current_Node_Value = reshape(curr_x(Node_globalIndex),2,3);
              Dx = Initial_Node_Value*N_triangle ; 
       Dx_inv = inv(Dx);
       F = Current_Node_Value*N_triangle* Dx_inv; 
       
         K_local = Ai(iter)*get_local_stiffness_triangle (Current_Node_Value, Initial_Node_Value,E,nu);
    %    K_local =get_local_stiffness_triangle (Current_Node_Value, Initial_Node_Value,E,nu);

       Je = (Dx)';
      %  K_local = Local_2D_stiffness(F,Je,kk,miu); % This is using FEM approach   
       for kk=1:1:length(Node_globalIndex)
           for hh=1:1:length(Node_globalIndex)
             Stiffnes(Node_globalIndex(kk),Node_globalIndex(hh))=Stiffnes(Node_globalIndex(kk),Node_globalIndex(hh))+...
                 K_local(kk,hh);
           end
       end
end
toc
% Modification Matrix Constraint 

Modified_Stiffness = ModificationStiffness_constraint (Stiffnes,Nf1_Node_globalIndex,Nf2_Node_globalIndex);



Stiffnes2 = Stiffnes;
Stiffnes2(:,NBc_Node_globalIndex)=[];
Stiffnes2(NBc_Node_globalIndex,:)=[];


rank(full(Stiffnes2))

NodeNumber_whole(NBC)=[];

% for i=1:length(NM)
%     NM_revised(i)= find(NodeNumber_whole==NM(i));
% end
% 
% for i=1:length(NP)
%     NP_revised(i)= find(NodeNumber_whole==NP(i));
% end

Nf1_Node_globalIndex2 = get_Node_globalIndex_dircetion2D(Nf1,2);
Nf2_Node_globalIndex2 = get_Node_globalIndex_dircetion2D(Nf2,2);
Nf1_Node_globalIndex1 = get_Node_globalIndex_dircetion2D(Nf1,1);
Nf2_Node_globalIndex1 = get_Node_globalIndex_dircetion2D(Nf2,1);
NF_Node_globalIndex2 = [Nf1_Node_globalIndex2,Nf2_Node_globalIndex2];
NF_Node_globalIndex_total = [Nf1_Node_globalIndex1,Nf1_Node_globalIndex2...
                             Nf2_Node_globalIndex1,Nf2_Node_globalIndex2];

NP_Node_globalIndex = get_Node_globalIndex2D(NP);



% NM_revised_globalIndex2 = get_Node_globalIndex_dircetion2D(NM_revised,1);
% NP_revised_globalIndex = get_Node_globalIndex2D(NP_revised);
%%
NodeNumber_whole = (1:1:size(Nodes,2));
Index = find(ismember(NodeNumber_whole,[Nf1,Nf2,NP,NBC]));
NodeNumber_whole(Index)=[];
NR_Node_globalIndex = get_Node_globalIndex2D(NodeNumber_whole);
Ordering_Vec = [NP_Node_globalIndex,NF_Node_globalIndex,NR_Node_globalIndex,NBc_Node_globalIndex,...
                Nf1_Node_globalIndex1,Nf2_Node_globalIndex1];
Stiffnes2 = Stiffnes2(Ordering_Vec,Ordering_Vec);

% figure 
% spy(Stifness)
% figure 
% spy(Stifness2)


Stiffnes2(:,end-length([NBc_Node_globalIndex,Nf1_Node_globalIndex1,Nf2_Node_globalIndex1]):end)=[];
Stiffnes2(end-length([NBc_Node_globalIndex,Nf1_Node_globalIndex1,Nf2_Node_globalIndex1]):end,:)=[];
Stiffnes2(length(NP_Node_globalIndex)+1:length(NP_Node_globalIndex)+...
    length(NF_Node_globalIndex),:)=[];
% for this specific 2D problem 
% Stifness(1:2,:)=[];
% Stifness(:,1:2)=[];

% _____________________
K1 = Stiffnes2(:,1:length(NP_Node_globalIndex));
K2 = Stiffnes2(:,length(NP_Node_globalIndex)+1:end);
Jacobian = -pinv(full(K2))*K1;
nm = length([Nf1_Node_globalIndex2,...
        Nf2_Node_globalIndex2]);
DX = X_True2-X_True1;

dcc=(Jacobian)*DX;
norm(dcc(1:nm)-dxx)/norm(dxx)
 
%%

% K1 = Stifness(:,1:length(NP_Node_globalIndex));
% K2 = Stifness(:,length(NP_Node_globalIndex)+1:end);

nm = length([Nf1_Node_globalIndex2,...
        Nf2_Node_globalIndex2]);
% A = K2(1:nm,1:nm);
% B = K2(1:nm,nm+1:end);
% C = K2(nm+1:end,1:nm);
% D = K2(nm+1:end,nm+1:end);

%Var_P = rand(length(NP_Node_globalIndex)-2,1);
Var_P=DX;

f=K1*Var_P;
a=f(1:nm);
b= f(nm+1:end);

size(D)
spy(A)

% Note this should be length(Np)=length(Nf1)+length(Nf2);
tic
Jacobian = -pinv(full(K2))*K1;
dc = Jacobian*DX;
toc
 
e=pinv(full(K2))*K2;
e(1:10,1:10)

tic 
e=B*inv(D);
Jacobian2 = pinv(full(A-e*C))*(a-e*b);
toc 

norm(dc(1:nm)-dxx)/norm(dxx)


pinv(full(A-e*C))*full(A-e*C)



%% Another Jacobian (This is the Jacobian I should use)
DX = X_True2-X_True1;
dxx= [ones(6,1)*-0.003;ones(6,1)*0.002];
dx2=[-0.003,0.002];
% Jf = sparse(size(Nodes,2)*2,length(NF_Node_globalIndex_total));
% for k=1:1:length(NF_Node_globalIndex_total)
%     Jf(NF_Node_globalIndex_total(k),k)=1;
% end
% Jf(NBc_Node_globalIndex,:)=[];



Jf = sparse(size(Nodes,2)*2,length(NF_Node_globalIndex2));
for k=1:1:length(NF_Node_globalIndex2)
    Jf(NF_Node_globalIndex2(k),k)=1;
end
Jf(NBc_Node_globalIndex,:)=[];


Jc = sparse(length(NF_Node_globalIndex2),size(Nodes,2)*2);
for k=1:1:length(NF_Node_globalIndex2)
    Jc(k,NF_Node_globalIndex2(k))=1;
end
Jc(:,NBc_Node_globalIndex)=[];

Jt = sparse(length(NP_Node_globalIndex),size(Nodes,2)*2);
for k=1:1:length(NP_Node_globalIndex)
    Jt(k,NP_Node_globalIndex(k))=1;
end
Jt(:,NBc_Node_globalIndex)=[];

%  Jc = sparse(length(NM_revised_globalIndex2),size(Stiffnes,2));
%  for k=1:1:length(NM_revised_globalIndex2)
%      Jc(k,NM_revised_globalIndex2(k))=1;
%  end
%  Jt = sparse(length(NP_revised_globalIndex),size(Stiffnes,2));
%  for k=1:1:length(NP_revised_globalIndex)
%      Jt(k,NP_revised_globalIndex(k))=1;
%  end
% 

Stiffnes00 = Stiffnes2;
K_inv = inv((Stiffnes00));
% Old version 
% S1 = inv(Jc*K_inv*Jc');
% S2 = Jt*K_inv*Jc';
% Jacobian3 = inv(full(S2*S1));


% New version 
S1 = Jc*K_inv*Jf;
AA = Jt*K_inv*Jf;
S2 = pinv(full(Jt*K_inv*Jf));
Jacobian3 = S1*S2;
Jprime= pinv(full(Jacobian3));
landa=0.001;
Jacobi_modi = Jprime'*inv(Jprime*Jprime'+landa^2*eye(size(Jprime*Jprime')));


Jacobian2 = AA*inv(S1);
kl = [1,0;1,0;1,0;1,0;1,0;1,0;...
      0,1;0,1;0,1;0,1;0,1;0,1;];
JJ=  Jacobian2*kl;
JJ_modi = JJ'*inv(JJ*JJ'+landa^2*eye(size(JJ*JJ')));

% Jacobi_modi2 = Jacobian3*inv(Jacobian3'*Jacobian3+landa^2*eye(38,38));

dcc=((Jacobian3))*DX;
dcc0=Jacobi_modi*DX;
dcc00=pinv(JJ)*DX(1:2);
% dcc00=JJ_modi*DX(1:2);


norm(dcc0-dxx)/norm(dxx)
norm(dcc00-dx2)/norm(dx2)


