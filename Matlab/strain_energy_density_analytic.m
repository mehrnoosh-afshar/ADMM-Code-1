function Say = strain_energy_density_analytic (x,N_tri,Dx_inv,k,miu)
F = (x*N_tri)* Dx_inv; 
j = det(F);
i = trace(F'*F);
Say = 0.5*miu*(j^(-2/3)*i-3)+0.5*k*(j-1)^2;

end