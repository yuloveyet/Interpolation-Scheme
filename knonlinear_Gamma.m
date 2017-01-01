function [ke_Gamma,Int_Gamma]=knonlinear_Gamma(ele,xxyy,D,h,ue,Gamma)
ke_Gamma=zeros(8,8);
% ke_l=zeros(8,8);
M_Gamma=zeros(4,4);
PK2_Gamma=zeros(3,1);
Int_Gamma=zeros(8,1);
gauss_point=1/sqrt(3);
 Id2=eye(2);
 GS=[gauss_point;-gauss_point];
 for i=1:2
     s=GS(i);
     for j=1:2
         t=GS(j);
         [~,B_Gamma,detJ,~,E_Gamma,~,G,~,~]=shape_Gamma(s,t,ele,xxyy,ue,Gamma);
         PK2_Gamma=D*E_Gamma;
         M_Gamma=[PK2_Gamma(1)*Id2 PK2_Gamma(3)*Id2
             PK2_Gamma(3)*Id2 PK2_Gamma(2)*Id2];
         ke_Gamma=ke_Gamma+h*detJ*(B_Gamma'*D*B_Gamma+G'*M_Gamma*G);
         Int_Gamma=Int_Gamma+h*detJ*B_Gamma'*D*E_Gamma;
%          ke_l=ke_l+h*detJ*B0'*D*B0;
     end
 end

end