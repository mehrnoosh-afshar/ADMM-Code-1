%%%%%%%%%% Second Version of Code revise convergence 
clear all
%--------------- First make a 2D mesh ---------------------------------
model = createpde;
R1 = [3,4,0,0.2,0.2,0,0.05,0.05,-0.05,-0.05]';
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
Nf1_0 = findNodes(mesh,'box',[0.05 0.08],[0.05 0.05]);
Nf1_1 = findNodes(mesh,'box',[0.13 0.16],[0.05 0.05]);
Nf2_0 = findNodes(mesh,'box',[0.05 0.08],[-0.05 -0.05]);
Nf2_1 = findNodes(mesh,'box',[0.13 0.16],[-0.05 -0.05]);
Nf1 = [Nf1_0,Nf1_1];
Nf2 = [Nf2_0,Nf2_1];
NP = findNodes(mesh,'radius',[0.1 0],.004);
figure
pdemesh(model)
hold on
plot(Nodes(1,NBC),Nodes(2,NBC),'ok','MarkerFaceColor','g') 
plot(Nodes(1,Nf1),Nodes(2,Nf1),'ok','MarkerFaceColor','r') 
plot(Nodes(1,Nf2),Nodes(2,Nf2),'ok','MarkerFaceColor','r') 

%% ------------ Some pre-loop  calculations --------------------------------
E=3000;
nu=0.49;
Nf1_Node_globalIndex2 = get_Node_globalIndex_dircetion2D(Nf1,2);
Nf2_Node_globalIndex2 = get_Node_globalIndex_dircetion2D(Nf2,2);
Nf1_Node_globalIndex2 = get_Node_globalIndex_dircetion2D([Nf1_0,Nf1_1],2);
Nf2_Node_globalIndex2 = get_Node_globalIndex_dircetion2D([Nf2_0,Nf2_1],2);
NF_Node_globalIndex2 = [Nf1_Node_globalIndex2,Nf2_Node_globalIndex2];
NP_Node_globalIndex = get_Node_globalIndex2D(NP);
NBc_Node_globalIndex = get_Node_globalIndex2D(NBC);
NumerLoop = 50;
curr_x = reshape(Nodes,[],1);
%% Define the desire position of voundary point 
% U is the displacement of each boundary 
x_desire0 = [[(curr_x(Nf1_Node_globalIndex2(1))),(curr_x(Nf1_Node_globalIndex2(4)))]*[1,1,1,0,0,0;0,0,0,1,1,1],...
              [(curr_x(Nf2_Node_globalIndex2(1))),(curr_x(Nf2_Node_globalIndex2(4)))]*[1,1,1,0,0,0;0,0,0,1,1,1]]';
U(1,1)=-0.01;
U(2,1)=0.01;
U(3,1)=0.01;
U(4,1)=0;
x_desire = x_desire0 + [[U(1,1),U(2,1)]*[1,1,1,0,0,0;0,0,0,1,1,1],[U(3,1),U(4,1)]*[1,1,1,0,0,0;0,0,0,1,1,1]]';
%% ADMM function 
[Current_configuration,S,S0] = ADMM_2D(x_desire,NumerLoop,Elements,Nodes,Ai,At,NBC,Nf1,Nf2,NP);
%% Plot Data 
nodes = reshape(Current_configuration(:,1),2,size(Nodes,2));
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

%% Plot Mesh 
% for this elemnts should not collide with each other 
nodes2 = nodes;
Nodes2 = Nodes;
nodes2(2,:) = nodes2(2,:);
Nodes2(2,:) = Nodes2(2,:);
model2 = createpde();
mesh2=geometryFromMesh(model2,nodes2,Elements);
model3 = createpde();
mesh3=geometryFromMesh(model3,Nodes2,Elements);
figure;
pdemesh(model2);
hold on 
%pdemesh(model3,'EdgeColor', 'black')
b=pdegplot(model3)
b.Color = [0 0 0];
b.LineWidth = 1;
hold on 
p1=plot(Nodes2(1,NP),Nodes2(2,NP),'ok','MarkerFaceColor','r') 
hold on 
p2=plot(nodes2(1,NP),nodes2(2,NP),'ok','MarkerFaceColor','g') 
hold on 
p3=plot(nodes2(1,[Nf1]),nodes2(2,[Nf1]),'ok','MarkerFaceColor','b') 
hold on 
plot(nodes2(1,[Nf2]),nodes2(2,[Nf2]),'ok','MarkerFaceColor','b') 
xlabel('x(m)','FontSize',12,'FontWeight','bold')
ylabel('y(m)','FontSize',12,'FontWeight','bold')
title ('Second Scenario')
legend([p1 p2 p3],{'First location of target point','Final location of target point','Actuated points'})
ax = gca;
ax.FontWeight = 'bold';


          