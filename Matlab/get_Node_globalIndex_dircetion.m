function Node_globalIndex_direction = get_Node_globalIndex_dircetion(NodeNumber,direction_index)
% for x direction direction_index=1 
% for y direction direction_index=2;
% for z direction direction_index=3;
l=1;

       for i=1:1:size(NodeNumber,2)
           
            switch direction_index
                case 1
                   Node_globalIndex_direction(l) = 3*NodeNumber(i)-2;
                   l=l+1;
                case 2
                   Node_globalIndex_direction(l) = 3*NodeNumber(i)-1;
                   l=l+1;
                case 3
                   Node_globalIndex_direction(l) = 3*NodeNumber(i);
                   l=l+1;
            otherwise
                disp('wrong direction_index')
            end
       end  
       
end 