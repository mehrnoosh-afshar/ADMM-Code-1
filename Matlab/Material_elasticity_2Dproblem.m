function [Tangent_stifness , Second_Piola_Stress] = Material_elasticity_2Dproblem (C,kk,miu)
J1=sqrt(det(C));
I1_0 = trace(C);
C33 = I1_0/3 ;
%C33 = 1/(J1^2) ;

CR = zeros(3,3);
CR(1:2,1:2)=C;
CR(3,3)=C33;
  CR = C;
J=sqrt(det(CR));
I1 = trace(CR);


C_inv = inv(CR);

CT = tensor(C_inv);
IT = tensor(eye(size(CR)));
m = size(CR,2);
% for i=1:1:m
%     for j=1:1:m
%         for k=1:1:m
%             for l=1:1:m
%                 II(i,j,k,l)= 0.5*(C_inv(i,k)*C_inv(j,l)+C_inv(j,l)*C_inv(j,k));
%             end
%         end
%     end
% end

II = ttt(IT,IT);

C_had = 2*miu*J^(-2/3)*((1/3)*I1*II-(1/3)*ttt(IT,CT)-(1/3)*ttt(CT,IT)+(1/9)*I1*ttt(CT,CT));
C_p = kk*(J-1)*J*(ttt(CT,CT)-2*II);
C_k = kk *J^2*(ttt(CT,CT));

C_total = C_had + C_p + C_k; 
% Voigt notation 
% This Tanget Stiffness Matrix is Onlu for 2D problem 
Tangent_stifness = [C_total(1,1,1,1),C_total(1,1,2,2),C_total(1,1,1,2);
    C_total(1,1,2,2),C_total(2,2,2,2),C_total(2,2,1,2);
    C_total(1,1,1,2),C_total(2,2,1,2),C_total(1,2,1,2)];

S =  miu* J^(-2/3)*(eye(m,m)-(1/3)*I1*C_inv);
Second_Piola_Stress = [S(1:2,1:2),zeros(2,2);zeros(2,2),S(1:2,1:2)];

end