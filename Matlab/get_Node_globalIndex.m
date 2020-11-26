function Node_globalIndex = get_Node_globalIndex(NodeNumber)


       for i=1:1:size(NodeNumber,2)
           Node_globalIndex(3*i-2:3*i) = [3*NodeNumber(i)-2,3*NodeNumber(i)-1,3*NodeNumber(i)];
       end  
       
end 