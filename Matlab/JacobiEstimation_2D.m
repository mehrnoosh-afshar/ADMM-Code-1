function Jacobian = JacobiEstimation_2D(curr_x,Nodes,Elements,E,nu,Ai,Jf,Jc,Jt,kl,NBc_Node_globalIndex)
% Jacobian Estimation based On Quasi-Static equation 
Stiffnes = sparse (size(Nodes,2)*2,size(Nodes,2)*2);

% Build Global Stiffness Matrix form Local stifnesses 
N_triangle = [1,0;0,1;-1,-1]; % This is for triangle mesh
kk = E/(3*(1-2*nu));
miu = E/(2*(1+nu));

for iter=1:1:size(Elements,2)
       NodeNumber = Elements(:,iter);
       Node_globalIndex = get_Node_globalIndex2D(NodeNumber');

       Initial_Node_Value = Nodes(:,NodeNumber);
       Current_Node_Value = reshape(curr_x(Node_globalIndex),2,3);
              Dx = Initial_Node_Value*N_triangle ; 
       Dx_inv = inv(Dx);
       F = Current_Node_Value*N_triangle* Dx_inv; 
       
         K_local = Ai(iter)*get_local_stiffness_triangle (Current_Node_Value, Initial_Node_Value,E,nu);

       Je = (Dx)';
       for kk=1:1:length(Node_globalIndex)
           for hh=1:1:length(Node_globalIndex)
             Stiffnes(Node_globalIndex(kk),Node_globalIndex(hh))=Stiffnes(Node_globalIndex(kk),Node_globalIndex(hh))+...
                 K_local(kk,hh);
           end
       end
end

Stiffnes(:,NBc_Node_globalIndex)=[];
Stiffnes(NBc_Node_globalIndex,:)=[];

% Jacobian calculation 

K_inv = inv((Stiffnes));

S1 = Jc*K_inv*Jf;
AA = Jt*K_inv*Jf;

Jacobian2 = AA*inv(S1);

JJ=  Jacobian2*kl;

Jacobian = JJ;


end