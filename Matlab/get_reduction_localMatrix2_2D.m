function DDi = get_reduction_localMatrix2_2D(Node_globalIndex,sizeX)
     
     % Node_globalIndex is a vector telling the global index of elemnrt's
     % node 
     Ng = Node_globalIndex;
     Row = [1,1,2,2,3,3,4,4];
        
     Colum = [Ng(1),Ng(3),Ng(1),Ng(5),...
              Ng(2),Ng(4),Ng(2),Ng(6)];
     Value = [1,-1,1,-1,1,-1,1,-1];
     
     m = 4; n = sizeX;

     DDi = -sparse(Row,Colum,Value,m,n);
     
end