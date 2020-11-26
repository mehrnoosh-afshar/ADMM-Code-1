function Node_globalIndex = get_Node_globalIndex2D(NodeNumber)


       for i=1:1:size(NodeNumber,2)
           Node_globalIndex(2*i-1:2*i) = [2*NodeNumber(i)-1,2*NodeNumber(i)];
       end  
       
end

