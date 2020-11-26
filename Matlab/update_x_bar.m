function x_bar = update_x_bar(x_n,v_n,Dt,M_inv,F_extr)
 x_bar = x_n + v_n * Dt + M_inv * F_extr* Dt^2;
end