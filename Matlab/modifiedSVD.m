function [U,sigma_bar,V] = modifiedSVD(F_bar)

     [U,sigma_bar,V]=svd(F_bar);
     J = eye(3,3);
     J(3,3) = -1;
     if det(U) < 0 
         U = U*J;
         sigma_bar(3,3) = -sigma_bar(3,3);
     end
     
     if det(V) < 0 
         V = (J*V')';
         sigma_bar(3,3) = -sigma_bar(3,3);
     end
end