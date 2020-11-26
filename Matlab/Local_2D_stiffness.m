function Element_triangle_stiffness = Local_2D_stiffness(F,Je,kk,miu) 

    C =F'*F;
    [Tangent_stifness , Second_Piola_Stress] = Material_elasticity_2Dproblem (C,kk,miu);
     
    DU = F -eye(size(F));
    BL0 = 1/det(Je)*[Je(2,2),0,-Je(1,2),0,-Je(2,2)+Je(1,2),0;
                 0,-Je(2,1),0,Je(1,1),0,Je(2,1)-Je(1,1);
                 -Je(2,1),Je(2,2),Je(1,1),-Je(1,2),Je(2,1)-Je(1,1),-Je(2,2)+Je(1,2)];
             
           
    BL1 = 1/det(Je)*[DU(1,1)*Je(2,2),DU(2,1)*Je(2,2),-DU(1,1)*Je(1,2),-DU(2,1)*Je(1,2),DU(1,1)*(-Je(2,2)+Je(1,2)),DU(2,1)*(-Je(2,2)+Je(1,2));
                 -DU(1,2)*Je(2,1),-DU(2,2)*Je(2,1),DU(1,2)*Je(1,1),DU(2,2)*Je(1,1),DU(1,2)*(Je(2,1)-Je(1,1)),DU(2,2)*(Je(2,1)-Je(1,1));
                 -DU(1,1)*Je(2,1)+DU(1,2)*Je(2,2),-DU(2,1)*Je(2,1)+DU(2,2)*Je(2,2),...
                  DU(1,1)*Je(1,1)-DU(1,2)*Je(1,2), DU(2,1)*Je(1,1)-DU(2,2)*Je(1,2),...
                  DU(1,1)*(Je(2,1)-Je(1,1))+DU(1,2)*(-Je(2,2)+Je(1,2)), DU(2,1)*(Je(2,1)-Je(1,1))+DU(2,2)*(-Je(2,2)+Je(1,2))];
    BL = BL0 + BL1;    
              
    BNL = 1/det(Je)*[Je(2,2),0,-Je(1,2),0,-Je(2,2)+Je(1,2),0;
                  -Je(2,1),0,Je(1,1),0,Je(2,1)-Je(1,1),0;
                  0,Je(2,2),0,-Je(1,2),0,-Je(2,2)+Je(1,2);
                  0,-Je(2,1),0,Je(1,1),0,Je(2,1)-Je(1,1)];
              
    KL = BL'*Tangent_stifness*BL;
    KN = BNL'*Second_Piola_Stress*BNL;

    Element_triangle_stiffness = KL+KN;


end