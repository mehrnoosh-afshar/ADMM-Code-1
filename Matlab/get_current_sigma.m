function curr_sigma = get_current_sigma(curr_z,NumberofElemnts)
    for i=1:1:NumberofElemnts
        z =  curr_z((i-1)*9+1:i*9);
        Ft = reshape(z,3,3);
        [U,curr_sigma0,V] = modifiedSVD(Ft);
        curr_sigma((i-1)*3+1:i*3) =  diag(curr_sigma0);
    end
end