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
E=3000; % These are defined in c++ but it is good if we have them as input to Mexfile
nu=0.49;

% These are the files that I have read in c++ from mat file
Forceindex = [Nf1';Nf2'];
Forcevalue = zeros(size(Forceindex,1),1);
Forcevalue(1:size(Nf1,2),:)= -20;
Forcevalue(size(Nf1,2)+1:end,:)=  20;
Forcevalue = [zeros(size(Forcevalue)),Forcevalue];
Bcindex = NBC';
Node_corr = Nodes';  % This is data.mat in c++ code 
Element_index = Elements'; % This is data2.mat in c++ code 


%% We want to put a Mexfile Here 
Final_Node_Cordinate = main(Forceindex,Forcevalue,Bcindex,Node_corr,Element_index);
% Final_Node_Cordinate is the result of optimization 

%% Plot Data 
Current_configuration = Final_Node_Cordinate ; 
nodes = reshape(Current_configuration,2,size(Nodes,2));
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