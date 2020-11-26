function Modified_Stiffness = ModificationStiffness_constraint (Stiffnes,NF1_globalindex,NF2_globalindex)
    
% Cloumn manipulation
Stiffnes(:,NF1_globalindex(1))= Stiffnes(:,NF1_globalindex(1))+sum(Stiffnes(:,NF1_globalindex(2):end),2);
Stiffnes(:,NF2_globalindex(1))= Stiffnes(:,NF2_globalindex(1))+sum(Stiffnes(:,NF2_globalindex(2):end),2);


% Row manipulation 
Stiffnes(NF1_globalindex(1),:)= Stiffnes(NF1_globalindex(1),:)+sum(Stiffnes(NF1_globalindex(2):end,:),1);
Stiffnes(NF2_globalindex(1),:)= Stiffnes(NF2_globalindex(1),:)+sum(Stiffnes(NF2_globalindex(2):end,:),1);


Modified_Stiffness = Stiffnes;

end