function [Strain_Energy,Strain_Energy_gradiant0]=energyTetelement(X,z,miu,k)
   % X and x is 3*4 matrix
   % z is a 1*9 matrix so that 
   % z(1)= x(2,:)-x(1,:);
   % z(2)= x(3,:)-x(1,:);
   % z(3)= x(4,:)-x(1,:);
   B = [X(:,2)-X(:,1),X(:,3)-X(:,1),X(:,4)-X(:,1)];
   F = [z(1:3),z(4:6),z(7:9)]*inv(B);
   %F = [x(2,:)-x(1,:),x(3,:)-x(1,:),x(4,:)-x(1,:)]*inv(B);
   
   % Cauchy- Green strain tensor
   C = F'*F;
   % Determine constans 
  
   J = det(F);
   %if (J > 0)
   C_bar = J^(-2/3)*C;
   I1 = trace(C); 
   I1_bar = trace(C_bar);
   I2_bar = 0.5*(I1_bar^2-trace(C_bar^2));
   Strain_Energy = 0.5*miu*(I1_bar-3)+0.5*k*(J-1)^2;
   Strain_Energy_gradiant0 = miu* J^(-2/3)*(F-(1/3)*I1*(inv(F))')+ k*(J-1)*J*(inv(F))';
   Strain_Energy_gradiant = reshape(Strain_Energy_gradiant0,[],1);
  % else 
  %     disp("det F is negative :: volume problem")
  % end
end