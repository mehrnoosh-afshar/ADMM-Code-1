function Di = get_reduction_localMatrix(Node_globalIndex,sizeX)
     
     % Node_globalIndex is a vector telling the global index of elemnrt's
     % node 
     Ng = Node_globalIndex;
     Row = [1,1,2,2,3,3,4,4,5,5,6,6,...
            7,7,8,8,9,9];
     Colum = [Ng(1),Ng(4),Ng(2),Ng(5),Ng(3),Ng(6),...
              Ng(1),Ng(7),Ng(2),Ng(8),Ng(3),Ng(9),...
              Ng(1),Ng(10),Ng(2),Ng(11),Ng(3),Ng(12)];
     Value = [1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1];
     
     m = 9; n = sizeX;

     Di = -sparse(Row,Colum,Value,m,n);
     
end


