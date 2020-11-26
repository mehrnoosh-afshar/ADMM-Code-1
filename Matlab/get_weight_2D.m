function Wi = get_weight_2D(vi,k)
     wi = 1*sqrt(k*vi); %THIS WAS GOOD CHOCE BUT SLOW
  
  %   wi = 0.6*sqrt(k*vi); 

     Wi = wi*eye(4,4);
end