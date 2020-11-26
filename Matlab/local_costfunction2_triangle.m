function [Cost_function2,Gradient2 ]=local_costfunction2_triangle(sigma,sigma_bar,miu,k,taw,vi)

  
   % Determine constans 
   J=sigma(1)*sigma(2);
   N= sigma(1)^2+sigma(2)^2;
   n=(J)^(-2/3);
   I1_bar = n*(N);
   Strain_Energy = 0.5*miu*(I1_bar-3)+0.5*k*(J-1)^2;
%    Strain_Energy_gradiant =[miu*n*((-1/3)*N)*sigma(1)+sigma(1)+k*(J-1)*sigma(2)*sigma(3),...
%        miu*n*((-1/3)*N)*sigma(2)+sigma(2)+k*(J-1)*sigma(1)*sigma(3),...
%        miu*n*((-1/3)*N)*sigma(3)+sigma(3)+k*(J-1)*sigma(1)*sigma(2)]';  %
%       check the gradiant for triangle 
   
   
   t1 = sigma(1)-sigma_bar(1);
   t2 = sigma(2)-sigma_bar(2);
   
   second_term =  0.5*taw*((t1)^2+(t2)^2);
   %second_term_gradiant = taw*[t1,t2]';
%    Cost_function2 = vi*(Strain_Energy + second_term); % the first series of code was using this 
     Cost_function2 = vi*(Strain_Energy) + second_term;
   % Gradient2 = Strain_Energy_gradiant + second_term_gradiant;
   Gradient2=1;
end