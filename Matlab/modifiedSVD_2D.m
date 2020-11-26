function [U,sigma_bar,V] = modifiedSVD_2D(F_bar)

     [U,sigma_bar,V]=svd(F_bar);
     J = eye(2,2);
     J(2,2) = -1;
     if det(U) < 0 
         U = U*J;
         sigma_bar(2,2) = -sigma_bar(2,2);
     end
     
     if det(V) < 0 
         V = (J*V')';
         sigma_bar(2,2) = -sigma_bar(2,2);
     end
end