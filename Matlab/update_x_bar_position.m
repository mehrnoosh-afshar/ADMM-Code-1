function x_bar = update_x_bar_position(x_n,v_n,Dt,x_desire,Index_control)
 x_bar = x_n + v_n * Dt ;
%  DX = (x_n(Index_control)- x_desire);
 x_bar(Index_control) = x_desire;
end