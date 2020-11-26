% % Controll Loop 
% % Model 2D 
% clear all var 
% model = createpde;
% R1 = [3,4,0,2,2,0,0.5,0.5,-0.5,-0.5]';
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
% [At,Ai] = area(mesh);
% 
% % Find Surface Nodes 
% NBC = findNodes(mesh,'region','Edge',4);
% % Nf1 = findNodes(mesh,'box',[0.5 1],[0.5 0.5]);
% Nf1_0 = findNodes(mesh,'box',[0.5 0.8],[0.5 0.5]);
% Nf1_1 = findNodes(mesh,'box',[1.3 1.6],[0.5 0.5]);
% 
% % Nf2 = findNodes(mesh,'box',[0.5 1.05],[-0.5 -0.5]);
% Nf2_0 = findNodes(mesh,'box',[0.5 0.8],[-0.5 -0.5]);
% Nf2_1 = findNodes(mesh,'box',[1.3 1.6],[-0.5 -0.5]);
% Nf1 = [Nf1_0,Nf1_1];
% Nf2 = [Nf2_0,Nf2_1];
% %  NP = findNodes(mesh,'radius',[1.5 0],.04);
% NP = findNodes(mesh,'radius',[1 0],.04);
% 
% % NM =[Nf1,Nf2];
% 
% figure
% pdemesh(model)
% hold on
% plot(Nodes(1,NBC),Nodes(2,NBC),'ok','MarkerFaceColor','g') 
%  plot(Nodes(1,Nf1),Nodes(2,Nf1),'ok','MarkerFaceColor','r') 
% % plot(Nodes(1,Nf1_0),Nodes(2,Nf1_0),'ok','MarkerFaceColor','r') 
% % plot(Nodes(1,Nf1_1),Nodes(2,Nf1_1),'ok','MarkerFaceColor','r') 
% 
%  plot(Nodes(1,Nf2),Nodes(2,Nf2),'ok','MarkerFaceColor','r') 
% % plot(Nodes(1,Nf2_0),Nodes(2,Nf2_0),'ok','MarkerFaceColor','r') 
% % plot(Nodes(1,Nf2_1),Nodes(2,Nf2_1),'ok','MarkerFaceColor','r') 
% 
% 
% plot(Nodes(1,NP),Nodes(2,NP),'ok','MarkerFaceColor','b') 
% %set(gca,'XTick',[], 'YTick', [])
% 
% 
% Target_initial_position = [Nodes(1,NP),Nodes(2,NP)];
% Target_desired_position = [Nodes(1,NP)+0.0,Nodes(2,NP)+0.1];
% Number_intervals =1;
% Target_Movement_Intervals = (Target_desired_position - Target_initial_position)/Number_intervals;
%%
% clear all var 
% model = createpde;
% R1 = [3,4,0,2,2,0,0.5,0.5,-0.5,-0.5]';
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
% [At,Ai] = area(mesh);
% 
% % Find Surface Nodes 
% NBC = findNodes(mesh,'region','Edge',4);
% Nf1 = findNodes(mesh,'box',[0.5 1],[0.5 0.5]);
% %Nf1_0 = findNodes(mesh,'box',[0.5 0.8],[0.5 0.5]);
% %Nf1_1 = findNodes(mesh,'box',[1.3 1.6],[0.5 0.5]);
% 
% Nf2 = findNodes(mesh,'box',[0.5 1.05],[-0.5 -0.5]);
% %Nf2_0 = findNodes(mesh,'box',[0.5 0.8],[-0.5 -0.5]);
% %Nf2_1 = findNodes(mesh,'box',[1.3 1.6],[-0.5 -0.5]);
% %Nf1 = [Nf1_0,Nf1_1];
% %Nf2 = [Nf2_0,Nf2_1];
%   NP = findNodes(mesh,'radius',[1.5 0],.04);
% Nf3 = findNodes(mesh,'box',[2 2],[-0.3 0.3]);
% 
% % NP = findNodes(mesh,'radius',[1 0],.04);
% 
% % NM =[Nf1,Nf2];
% 
% figure
% pdemesh(model)
% hold on
% plot(Nodes(1,NBC),Nodes(2,NBC),'ok','MarkerFaceColor','g') 
%  plot(Nodes(1,Nf1),Nodes(2,Nf1),'ok','MarkerFaceColor','r') 
% % plot(Nodes(1,Nf1_0),Nodes(2,Nf1_0),'ok','MarkerFaceColor','r') 
% % plot(Nodes(1,Nf1_1),Nodes(2,Nf1_1),'ok','MarkerFaceColor','r') 
% 
%  plot(Nodes(1,Nf2),Nodes(2,Nf2),'ok','MarkerFaceColor','r') 
% % plot(Nodes(1,Nf2_0),Nodes(2,Nf2_0),'ok','MarkerFaceColor','r') 
% % plot(Nodes(1,Nf2_1),Nodes(2,Nf2_1),'ok','MarkerFaceColor','r') 
%  plot(Nodes(1,Nf3),Nodes(2,Nf3),'ok','MarkerFaceColor','r') 
% 
% 
% plot(Nodes(1,NP),Nodes(2,NP),'ok','MarkerFaceColor','b') 
% %set(gca,'XTick',[], 'YTick', [])
% 
% 
% Target_initial_position = [Nodes(1,NP),Nodes(2,NP)];
% Target_desired_position = [Nodes(1,NP)+0.1,Nodes(2,NP)+0.05];
% Number_intervals =1;
% Target_Movement_Intervals = (Target_desired_position - Target_initial_position)/Number_intervals;
%%
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
[At,Ai] = area(mesh);

% Find Surface Nodes 
NBC = findNodes(mesh,'region','Edge',4);
% Nf1 = findNodes(mesh,'box',[0.5 1],[0.5 0.5]);
Nf1_0 = findNodes(mesh,'box',[0.05 0.08],[0.05 0.05]);
Nf1_1 = findNodes(mesh,'box',[0.13 0.16],[0.05 0.05]);

% Nf2 = findNodes(mesh,'box',[0.5 1.05],[-0.5 -0.5]);
Nf2_0 = findNodes(mesh,'box',[0.05 0.08],[-0.05 -0.05]);
Nf2_1 = findNodes(mesh,'box',[0.13 0.16],[-0.05 -0.05]);
Nf1 = [Nf1_0,Nf1_1];
Nf2 = [Nf2_0,Nf2_1];
%  NP = findNodes(mesh,'radius',[1.5 0],.04);
NP = findNodes(mesh,'radius',[0.1 0],.004);

% NM =[Nf1,Nf2];

figure
pdemesh(model)
hold on
plot(Nodes(1,NBC),Nodes(2,NBC),'ok','MarkerFaceColor','g') 
 plot(Nodes(1,Nf1),Nodes(2,Nf1),'ok','MarkerFaceColor','r') 
% plot(Nodes(1,Nf1_0),Nodes(2,Nf1_0),'ok','MarkerFaceColor','r') 
% plot(Nodes(1,Nf1_1),Nodes(2,Nf1_1),'ok','MarkerFaceColor','r') 

 plot(Nodes(1,Nf2),Nodes(2,Nf2),'ok','MarkerFaceColor','r') 
% plot(Nodes(1,Nf2_0),Nodes(2,Nf2_0),'ok','MarkerFaceColor','r') 
% plot(Nodes(1,Nf2_1),Nodes(2,Nf2_1),'ok','MarkerFaceColor','r') 


plot(Nodes(1,NP),Nodes(2,NP),'ok','MarkerFaceColor','b') 
%set(gca,'XTick',[], 'YTick', [])


Target_initial_position = [Nodes(1,NP),Nodes(2,NP)];
Target_desired_position = [Nodes(1,NP)+0.0,Nodes(2,NP)+0.1];
Number_intervals =1;
Target_Movement_Intervals = (Target_desired_position - Target_initial_position)/Number_intervals;
%% ------------ Some pre-loop  calculations --------------------------------
E=3000;
nu=0.49;


Nf1_Node_globalIndex2 = get_Node_globalIndex_dircetion2D(Nf1,2);
Nf2_Node_globalIndex2 = get_Node_globalIndex_dircetion2D(Nf2,2);
% Nf3_Node_globalIndex1 = get_Node_globalIndex_dircetion2D(Nf3,1);
% Nf3_Node_globalIndex2 = get_Node_globalIndex_dircetion2D(Nf3,2);

Nf1_Node_globalIndex2 = get_Node_globalIndex_dircetion2D([Nf1_0,Nf1_1],2);
Nf2_Node_globalIndex2 = get_Node_globalIndex_dircetion2D([Nf2_0,Nf2_1],2);

NF_Node_globalIndex2 = [Nf1_Node_globalIndex2,Nf2_Node_globalIndex2];
% NF_Node_globalIndex2 = [Nf1_Node_globalIndex2,Nf2_Node_globalIndex2,Nf3_Node_globalIndex1];


NP_Node_globalIndex = get_Node_globalIndex2D(NP);
NBc_Node_globalIndex = get_Node_globalIndex2D(NBC);

%%
% This is Jacaobian Calculation 
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
% 
% this k1 is for first configuration 
% kl = [1,0;1,0;1,0;1,0;1,0;1,0;...
%        0,1;0,1;0,1;0,1;0,1;0,1;];

% for the second configuartion  

kl = [1,0,0,0;1,0,0,0;1,0,0,0;0,1,0,0;0,1,0,0;0,1,0,0;...
       0,0,1,0;0,0,1,0;0,0,1,0;0,0,0,1;0,0,0,1;0,0,0,1;];

%for third configuartion 
%  kl = [1,0,0,;1,0,0,;1,0,0,;1,0,0,;1,0,0,;1,0,0,;...
%         0,1,0;0,1,0;0,1,0;0,1,0;0,1,0;0,1,0;
%          0,0,1;0,0,1;0,0,1;0,0,1;0,0,1;0,0,1];
% Control Loop
NumerLoop = 50;
curr_x = reshape(Nodes,[],1);
Pgain = [1,0;0,1];
Igain = [1,0;0,1];
%%
% ll=1;
% for nn=1:Number_intervals
%     for kk=1:1:30
%      Jacobian = JacobiEstimation_2D(curr_x,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex);
%      DP =(Nodes(NP_Node_globalIndex)+nn*Target_Movement_Intervals)-curr_x(NP_Node_globalIndex)';
%      DC = Pgain*Jacobian*DP';
%      x_desire = [(curr_x(Nf1_Node_globalIndex2(1))+DC(1))*ones(size(Nf1_Node_globalIndex2,2),1);(curr_x(Nf2_Node_globalIndex2(1))+DC(2))*ones(size(Nf2_Node_globalIndex2,2),1)];  
%      
%      nodes = reshape(curr_x,2,size(Nodes,2));
%      curr_x = ADMM_2D(x_desire,NumerLoop,Elements,nodes,Ai,At,NBC,Nf1,Nf2,NP);
% 
%      history_x(:,kk)= curr_x;
% 
%      target_history1(:,kk) = curr_x(NP_Node_globalIndex);
%     end
% end
 
%%  MRAC (Model Refrence Adaptive Controller)
% P =eye(2)*0.05;
% Aref = -10*eye(2);
% Bref = eye(2);
% Cref= eye(2);
% Dref=eye(2)*0;
% sys = ss(Aref,Bref,Cref,Dref); 
% dt=0.1;
% Gx = 100*eye(2);
% Gr = 100*eye(2);
% Kx =-0.1*eye(2);
% Kr =-0.1*eye(2);
% for nn=1:Number_intervals
%    for   bb=1:1:5
%      Jacobian = JacobiEstimation_2D(curr_x,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex);
%      DP =(Nodes(NP_Node_globalIndex)+nn*Target_Movement_Intervals)-curr_x(NP_Node_globalIndex)';
%      B = inv(Jacobian);
%      X=curr_x(NP_Node_globalIndex);
%      t = 0:dt:bb*dt;
%      u =[1.5495;0.009]*ones(1,numel(t));
%      r = [1.5495;0.009];
%      [y,t,Xref0]=lsim(sys,u,t,Nodes(NP_Node_globalIndex));
%      Xref = Xref0(end,:);
%      e = X' - Xref;
%      Kx = Kx - Gx*X*e*P*B;
%      Kr = Kr - Gr*r*e*P*B;
%      DC = Kx*X + Kr*r;
%      
%      x_desire = [(curr_x(Nf1_Node_globalIndex2(1))+DC(1))*ones(size(Nf1_Node_globalIndex2,2),1);(curr_x(Nf2_Node_globalIndex2(1))+DC(2))*ones(size(Nf2_Node_globalIndex2,2),1)];
%      curr_x = ADMM_2D(x_desire,NumerLoop,Elements,nodes,Ai,At,NBC,Nf1,Nf2,NP);
%      
%      target_history(:,bb) = curr_x(NP_Node_globalIndex);
%     % bb = bb+1;
%    end
% end

%%
% Passivity Based Controller 

% bb=1;
% 
% for nn=1:Number_intervals
%      P =((Nodes(NP_Node_globalIndex)+nn*Target_Movement_Intervals)-curr_x(NP_Node_globalIndex)')';
% 
%     for kk=1:1:5
%      Jacobian = JacobiEstimation_2D(curr_x,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex);
%      DP =(Nodes(NP_Node_globalIndex)+nn*Target_Movement_Intervals)-curr_x(NP_Node_globalIndex)';
%          
%      DC = Jacobian*P;
%      C = 2*eye(2,2);
%      P = - 0.25*eye(2,2)*DP' + (eye(2,2)- C)*P
%      x_desire = [(curr_x(Nf1_Node_globalIndex(1))+DC(1))*ones(size(Nf1_Node_globalIndex,2),1);(curr_x(Nf2_Node_globalIndex(1))+DC(2))*ones(size(Nf2_Node_globalIndex,2),1)];
%      curr_x = ADMM_2D(x_desire,NumerLoop);
%      target_history2(:,bb) = curr_x(NP_Node_globalIndex);
%      bb = bb+1;
%     end
% end
%  

%% Linear MPC 


% alfa =0.9;
% r = 1;
% Ts = 0.1;
% H = 5;
% h = H*Ts;
% beta = exp(-r);
% 
% b = (h*alfa^h*log(alfa)-alfa^h+1)/(log(alfa)^2);
% a = (h^2*alfa^h-2*b)/log(alfa);
% c = (h*(alfa*beta)^h*log(alfa*beta)-(alfa*beta)^h+1)/(log(alfa*beta)^2);
% 
% Q = 0.5*eye(2,2);
% 
% nn=1;
% x_desire0 = [[(curr_x(Nf1_Node_globalIndex2(1))),(curr_x(Nf2_Node_globalIndex2(2)))]*[1,1,1,1,1,1,0,0,0,0,0,0;0,0,0,0,0,0,1,1,1,1,1,1]]';
%         
% Current_configuration = curr_x;
% 
%     for ii=1:1:35
%         target_history2(:,ii) = Current_configuration(NP_Node_globalIndex);
%         nodes = reshape(Current_configuration(:,1),2,size(Nodes,2));
%         history_x(:,ii)= Current_configuration;
%  
%      Jacobian = JacobiEstimation_2D(Current_configuration,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex);
%      DP =(Nodes(NP_Node_globalIndex)+nn*Target_Movement_Intervals)-Current_configuration(NP_Node_globalIndex)';
%      
%      J = (Jacobian);    
%      U = pinv(a*J + Jacobian'*Q)*(b-c)*DP';
%      DC(:,ii)=inv(J)*DP';
%      control(:,ii)=U(1:2,1);
% %     x_desire = [(curr_x(Nf1_Node_globalIndex2(1))+DC(1))*ones(size(Nf1_Node_globalIndex2,2),1);(curr_x(Nf2_Node_globalIndex2(1))+DC(2))*ones(size(Nf2_Node_globalIndex2,2),1)];
%      
%      
%      x_desire = x_desire0 + [[U(1,1),U(2,1)]*[1,1,1,1,1,1,0,0,0,0,0,0;0,0,0,0,0,0,1,1,1,1,1,1]]';
%      x_desire0 = x_desire;
%      Current_configuration = ADMM_2D(x_desire,NumerLoop,Elements,nodes,Ai,At,NBC,Nf1,Nf2,NP);
%      
%     end
% 
% 
% plot(target_history2(1,:))


%% Linear MPC without integrator  for first secanrio 
% Current_configuration = curr_x;
% 
%  x_desire0 = [[(Current_configuration(Nf1_Node_globalIndex2(1))),(Current_configuration(Nf2_Node_globalIndex2(2)))]*[1,1,1,1,1,1,0,0,0,0,0,0;0,0,0,0,0,0,1,1,1,1,1,1]]';
% %  
%  cc= history_x(:,ii);
%  x_desire0 = [[(cc(Nf1_Node_globalIndex2(1))),(cc(Nf2_Node_globalIndex2(2)))]*[1,1,1,1,1,1,0,0,0,0,0,0;0,0,0,0,0,0,1,1,1,1,1,1]]';
%  Current_configuration =cc;
% 
% 
% for ii=19:1:21
% target_history2(:,ii) = Current_configuration(NP_Node_globalIndex);
%  nodes = reshape(Current_configuration,2,size(Nodes,2));
%  history_x(:,ii)= Current_configuration;
% 
% Jacobian = JacobiEstimation_2D(Current_configuration,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex);
% 
% contrl_ability(ii)=rank(Jacobian);
% condition_number(ii)=cond(Jacobian);
% 
% % if condition_number(ii)>=2
% %     dd = condition_number(ii)/2;
% %     [u,s,v]=svd(Jacobian);
% %     s(2,2)=dd*s(2,2);
% %     Jacobian =  u*s*v;
% % end
% 
% % Jacobian = JacobiEstimation_2D(curr_x,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex);
% A = eye(2,2);
% B = (Jacobian);
% C = eye(2,2); 
% Np =5;
% Nc =5;
% 
% H = C;
% h = C; 
% F = C*A;
% f = C*A;
% for kk=2:Np
% h =h*A;
% H= vertcat(H,h);
% f = f*A;
% F= vertcat(F,f);
% end
% v=H*B;
% Phi=zeros(2*Np,2*Nc); %declare the dimension of Phi
% Phi(:,1:2)=v; % first column of Phi
% for i=2:Nc
% Phi(:,2*i-1:2*i)=[zeros(2*(i-1),2);v(1:2*(Np-i+1),:)]; %Toeplitz matrix
% end
% BarRs=[eye(2,2);eye(2,2);eye(2,2);eye(2,2);eye(2,2)];
% %BarRs=[eye(2,2);eye(2,2)];
% %BarRs=eye(2,2);
% 
% Phi_Phi= Phi'*Phi;
% Phi_F= Phi'*F;
% Phi_R=Phi'*BarRs;
% R = 10*eye(size(Phi_Phi));
% nn=1;
% setpoint =[(Nodes(NP_Node_globalIndex)+nn*Target_Movement_Intervals)]';
% %setpoint(2,1)=Current_configuration(NP_Node_globalIndex(2));
% %setpoint(1,1)=1.5777;
% 
% U =inv(Phi_Phi+R)*Phi'*(BarRs*setpoint-F*Current_configuration(NP_Node_globalIndex));
% 
% DC(:,ii)=inv(B)*(setpoint-Current_configuration(NP_Node_globalIndex));
% %dot(B*U(1:2,1),setpoint-Current_configuration(NP_Node_globalIndex))/(norm(B*U(1:2,1))*norm(setpoint-Current_configuration(NP_Node_globalIndex)));
% control(:,ii)=U(1:2,1);
% x_desire = x_desire0 + [[U(1,1),U(2,1)]*[1,1,1,1,1,1,0,0,0,0,0,0;0,0,0,0,0,0,1,1,1,1,1,1]]';
% Current_configuration = ADMM_2D(x_desire,NumerLoop,Elements,nodes,Ai,At,NBC,Nf1,Nf2,NP);
% x_desire0 = [[(Current_configuration(Nf1_Node_globalIndex2(1))),(Current_configuration(Nf2_Node_globalIndex2(2)))]*[1,1,1,1,1,1,0,0,0,0,0,0;0,0,0,0,0,0,1,1,1,1,1,1]]';
% 
% end 
% 
% % Current_configuration([Nf1_Node_globalIndex2,Nf2_Node_globalIndex2])
% Current_configuration(NP_Node_globalIndex)





%% %% Linear MPC without integrator  for second secanrio 
tic
x_desire0 = [[(curr_x(Nf1_Node_globalIndex2(1))),(curr_x(Nf1_Node_globalIndex2(4)))]*[1,1,1,0,0,0;0,0,0,1,1,1],...
              [(curr_x(Nf2_Node_globalIndex2(1))),(curr_x(Nf2_Node_globalIndex2(4)))]*[1,1,1,0,0,0;0,0,0,1,1,1]]';
Current_configuration = curr_x;
        
for ii=1:1:20
 target_history2(:,ii) = Current_configuration(NP_Node_globalIndex);
 nodes = reshape(Current_configuration(:,1),2,size(Nodes,2));
 history_x(:,ii)= Current_configuration;

Jacobian = JacobiEstimation_2D(Current_configuration,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex);
contrl_ability(ii)=rank(Jacobian);
condition_number(ii)=cond(Jacobian);
A = eye(2,2);
B = Jacobian;
C = eye(2,2); 
Np =5;
Nc =5;

H = C;
h = C; 
F = C*A;
f = C*A;
for kk=2:Np
h =h*A;
H= vertcat(H,h);
f = f*A;
F= vertcat(F,f);
end
v=H*B;
Phi=zeros(2*Np,4*Nc); %declare the dimension of Phi
Phi(:,1:4)=v; % first column of Phi
for i=2:Nc
Phi(:,4*(i-1)+1:4*i)=[zeros(2*(i-1),4);v(1:2*(Np-i+1),:)]; %Toeplitz matrix
end
M = 1*eye(2*Np,2*Np);
M(1:2,1:2)=[10,0;0,1];
M(3:4,3:4)=[10,0;0,1];
M(5:6,5:6)=[10,0;0,1];
M(7:8,7:8)=[10,0;0,1];
M(9:10,9:10)=[10,0;0,1];

BarRs=[eye(2,2);eye(2,2);eye(2,2);eye(2,2);eye(2,2)];
Phi_Phi= Phi'*M*Phi;
Phi_F= Phi'*M*F;
Phi_R=Phi'*M*BarRs;
R =1*eye(size(Phi_Phi));
nn=1;
setpoint = Target_desired_position';

U =inv(Phi_Phi+R)*(Phi_R*setpoint-Phi_F*Current_configuration(NP_Node_globalIndex));
control(:,ii)=U(1:4,1);
% x_desire = [[(curr_x(Nf1_Node_globalIndex2(1))+U(1,1)),(curr_x(Nf1_Node_globalIndex2(4))+U(2,1))]*[1,1,1,0,0,0;0,0,0,1,1,1],...
%              [(curr_x(Nf2_Node_globalIndex2(1))+U(3,1)),(curr_x(Nf2_Node_globalIndex2(4))+U(4,1))]*[1,1,1,0,0,0;0,0,0,1,1,1]]' ;

x_desire = x_desire0 + [[U(1,1),U(2,1)]*[1,1,1,0,0,0;0,0,0,1,1,1],[U(3,1),U(4,1)]*[1,1,1,0,0,0;0,0,0,1,1,1]]';
% x_desire0 = x_desire;
[Current_configuration,S,S0] = ADMM_2D(x_desire,NumerLoop,Elements,Nodes,Ai,At,NBC,Nf1,Nf2,NP);
x_desire0 = x_desire;    
end 
toc

U(1,1)=-0.01;
U(2,1)=0.01;
U(3,1)=0.02;
U(4,1)=0;

% Current_configuration([Nf1_Node_globalIndex2,Nf2_Node_globalIndex2])

%% Third configurstion 
% tic
% Current_configuration = curr_x;
%  x_desire0 = [[(Current_configuration(Nf1_Node_globalIndex2(1))),(Current_configuration(Nf2_Node_globalIndex2(1))),...
%               (Current_configuration(Nf3_Node_globalIndex1(1)))]...
%               *[1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0;...
%                 0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1]]';
% 
%         
% for ii=16:1:20
%  target_history3(:,ii) = Current_configuration(NP_Node_globalIndex);
%  nodes = reshape(Current_configuration(:,1),2,size(Nodes,2));
%  history_x(:,ii)= Current_configuration;
% 
% Jacobian = JacobiEstimation_2D(Current_configuration,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex);
% contrl_ability(ii)=rank(Jacobian);
% condition_number(ii)=cond(Jacobian);
% A = eye(2,2);
% B = Jacobian;
% C = eye(2,2); 
% Np =5;
% Nc =5;
% 
% H = C;
% h = C; 
% F = C*A;
% f = C*A;
% for kk=2:Np
% h =h*A;
% H= vertcat(H,h);
% f = f*A;
% F= vertcat(F,f);
% end
% v=H*B;
% Phi=zeros(2*Np,4*Nc); %declare the dimension of Phi
% Phi(:,1:3)=v; % first column of Phi
% for i=2:Nc
% Phi(:,3*(i-1)+1:3*i)=[zeros(2*(i-1),3);v(1:2*(Np-i+1),:)]; %Toeplitz matrix
% end
% M = 1*eye(2*Np,2*Np);
% M(1:2,1:2)=[30,0;0,1];
% M(3:4,3:4)=[30,0;0,1];
% M(5:6,5:6)=[30,0;0,1];
% M(7:8,7:8)=[30,0;0,1];
% M(9:10,9:10)=[30,0;0,1];
% 
% BarRs=[eye(2,2);eye(2,2);eye(2,2);eye(2,2);eye(2,2)];
% Phi_Phi= Phi'*M*Phi;
% Phi_F= Phi'*M*F;
% Phi_R=Phi'*M*BarRs;
% R =1*eye(size(Phi_Phi));
% nn=1;
% setpoint = Target_desired_position';
% 
% U =inv(Phi_Phi+R)*(Phi_R*setpoint-Phi_F*Current_configuration(NP_Node_globalIndex));
% control(:,ii)=U(1:3,1);
% % x_desire = [[(curr_x(Nf1_Node_globalIndex2(1))+U(1,1)),(curr_x(Nf1_Node_globalIndex2(4))+U(2,1))]*[1,1,1,0,0,0;0,0,0,1,1,1],...
% %              [(curr_x(Nf2_Node_globalIndex2(1))+U(3,1)),(curr_x(Nf2_Node_globalIndex2(4))+U(4,1))]*[1,1,1,0,0,0;0,0,0,1,1,1]]' ;
% 
% x_desire = x_desire0 + kl*U(1:3,1);
% % x_desire0 = x_desire;
% Current_configuration = ADMM_2D_third_version(x_desire,NumerLoop,Elements,Nodes,Ai,At,NBC,Nf1,Nf2,Nf3,NP);
% x_desire0 = x_desire;    
% end 
% toc
% 
% Current_configuration([Nf1_Node_globalIndex2,Nf2_Node_globalIndex2])


%%
nodes = reshape(Current_configuration(:,1),2,size(Nodes,2));
%nodes(2,:) = nodes(2,:)+0.5;
%Nodes(2,:) = Nodes(2,:)+0.5;

figure;
plot(nodes(1,:),nodes(2,:),'ok','MarkerFaceColor','r') 
hold on 
plot(Nodes(1,:),Nodes(2,:),'ok','MarkerFaceColor','g') 
hold on 
plot(nodes(1,Nf1_1),nodes(2,Nf2_1),'ok','MarkerFaceColor','b') 
hold on 
plot(nodes(1,NP),nodes(2,NP),'ok','MarkerFaceColor','y') 
hold on 
plot(Nodes(1,NP),Nodes(2,NP),'ok','MarkerFaceColor','b') 

%%
nodes2 = nodes;
Nodes2 = Nodes;
nodes2(2,:) = nodes2(2,:)+0.5;
Nodes2(2,:) = Nodes2(2,:)+0.5;
model2 = createpde();
mesh2=geometryFromMesh(model2,nodes2/10,Elements);
model3 = createpde();
mesh3=geometryFromMesh(model3,Nodes2/10,Elements);
figure;
pdemesh(model2);
hold on 
%pdemesh(model3,'EdgeColor', 'black')
b=pdegplot(model3)
b.Color = [0 0 0];
b.LineWidth = 1;
hold on 
p1=plot(Nodes2(1,NP)/10,Nodes2(2,NP)/10,'ok','MarkerFaceColor','r') 
hold on 
p2=plot(nodes2(1,NP)/10,nodes2(2,NP)/10,'ok','MarkerFaceColor','g') 
hold on 
p3=plot(nodes2(1,[Nf1])/10,nodes2(2,[Nf1])/10,'ok','MarkerFaceColor','b') 
hold on 
plot(nodes2(1,[Nf2])/10,nodes2(2,[Nf2])/10,'ok','MarkerFaceColor','b') 
xlabel('x(m)','FontSize',12,'FontWeight','bold')
ylabel('y(m)','FontSize',12,'FontWeight','bold')
title ('Second Scenario')
legend([p1 p2 p3],{'First location of target point','Final location of target point','Actuated points'})
ax = gca;
ax.FontWeight = 'bold';
%%
target_history4 =target_history2;
% for oo=1:5
% target_history3(:,30+oo)= target_history3(:,30);
% end
target_history4 =target_history4/10;
target_history4(2,:)=target_history4(2,:)+0.05;
%%
figure 
SamplingInstant = 0:1:20;
subplot(2,1,1)
xx=plot(SamplingInstant,target_history4(1,:),'black','LineWidth',1)

xlabel('Sampling Instant','FontSize',12,'FontWeight','bold')
ylabel('X position(m)','FontSize',12,'FontWeight','bold')
title ('Third Scenario','FontSize',12,'FontWeight','bold')
% ylim([0.1 0.1105])
ax=gca;
ax.FontWeight = 'bold';
subplot(2,1,2)
yy=plot(SamplingInstant,target_history4(2,:),'black','LineWidth',1)
xlabel('Sampling Instant','FontSize',12,'FontWeight','bold')
ylabel('Y position(m)','FontSize',12,'FontWeight','bold')
ax=gca;
ax.FontWeight = 'bold';

%% Linear MPC with Integrator 
curr_x0 = reshape(Nodes,[],1);
for ii=1:1:10

Jacobian = JacobiEstimation_2D(curr_x,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex);
A = [eye(2,2),zeros(2,2);eye(2,2),eye(2,2)];
B = [inv(Jacobian);inv(Jacobian)];
C = [zeros(2,2),eye(2,2)]; 
Np =5;
Nc =5;

H = C;
h = C; 
F = C*A;
f = C*A;
for kk=2:Np
h =h*A;
H= vertcat(H,h);
f = f*A;
F= vertcat(F,f);
end
v=H*B;
Phi=zeros(2*Np,2*Nc); %declare the dimension of Phi
Phi(:,1:2)=v; % first column of Phi
for i=2:Nc
Phi(:,2*i-1:2*i)=[zeros(2*(i-1),2);v(1:2*(Np-i+1),:)]; %Toeplitz matrix
end
BarRs=[[0,0,1,0;0,0,0,1];[0,0,1,0;0,0,0,1];[0,0,1,0;0,0,0,1];[0,0,1,0;0,0,0,1];[0,0,1,0;0,0,0,1]];
Phi_Phi= Phi'*Phi;
Phi_F= Phi'*F;
Phi_R=Phi'*BarRs;
R = 0*eye(size(Phi_Phi));
nn=1;
setpoint =[0;0;[(Nodes(NP_Node_globalIndex)+nn*Target_Movement_Intervals)]'];

DX = curr_x-curr_x0;
U =inv(Phi_Phi+R)*(Phi_R*setpoint-Phi_F*[DX(NP_Node_globalIndex);curr_x(NP_Node_globalIndex)]);
control(:,ii)=U(1:2,1);
x_desire = [[(curr_x(Nf1_Node_globalIndex2(1))+U(1,1)),(curr_x(Nf1_Node_globalIndex2(4))+U(2,1))]*[1,1,1,0,0,0;0,0,0,1,1,1],...
             [(curr_x(Nf2_Node_globalIndex2(1))+U(3,1)),(curr_x(Nf2_Node_globalIndex2(4))+U(4,1))]*[1,1,1,0,0,0;0,0,0,1,1,1]] ;
     curr_x0 = curr_x;
     curr_x = ADMM_2D(x_desire,NumerLoop);
     target_history2(:,ii) = curr_x(NP_Node_globalIndex);
end 