function [ke_Gamma,Int_Gamma,dGamma]=gamma_nonlinear(ele,xx,yy,D,aa,bb,h,ue,nely,Gamma)
ke_Gamma=zeros(8,8);keL_Gamma=zeros(8,8);
M_Gamma=zeros(4,4);ML_Gamma=zeros(4,4);
PK2_Gamma=zeros(3,1);Sigma_Gamma=zeros(3,1);
Int_Gamma=zeros(8,1);IntL_Gamma=zeros(8,1);
dGamma=zeros(8,1);dLGamma=zeros(8,1);
gauss_point=1/sqrt(3);
 Id2=eye(2);
%% point1
s=-gauss_point;
t=-gauss_point;
[~,B_Gamma,~,~,E_Gamma,~,G,Bl,~]=gamma_shape(s,t,ele,xx,yy,aa,bb,ue,nely,Gamma);
% Nonlinear  
PK2_Gamma=D*E_Gamma;
M_Gamma=[PK2_Gamma(1)*Id2 PK2_Gamma(3)*Id2 
         PK2_Gamma(3)*Id2 PK2_Gamma(2)*Id2];
ke_Gamma1=B_Gamma'*D*B_Gamma+G'*M_Gamma*G;
Int_Gamma1=B_Gamma'*D*E_Gamma;
dGamma1=Bl'*D*E_Gamma+B_Gamma'*D*0.5*Bl*Gamma(ele)*ue+B_Gamma'*D*E_Gamma/Gamma(ele);

%% point2
s=gauss_point;
t=-gauss_point;
[~,B_Gamma,~,~,E_Gamma,~,G,Bl,~]=gamma_shape(s,t,ele,xx,yy,aa,bb,ue,nely,Gamma);
% Nonlinear  
PK2_Gamma=D*E_Gamma;
M_Gamma=[PK2_Gamma(1)*Id2 PK2_Gamma(3)*Id2 
         PK2_Gamma(3)*Id2 PK2_Gamma(2)*Id2];
ke_Gamma2=B_Gamma'*D*B_Gamma+G'*M_Gamma*G;
Int_Gamma2=B_Gamma'*D*E_Gamma;
dGamma2=Bl'*D*E_Gamma+B_Gamma'*D*0.5*Bl*Gamma(ele)*ue;

%% point3
s=gauss_point;
t=gauss_point;
[~,B_Gamma,~,~,E_Gamma,~,G,Bl,~]=gamma_shape(s,t,ele,xx,yy,aa,bb,ue,nely,Gamma);
% Nonlinear  
PK2_Gamma=D*E_Gamma;
M_Gamma=[PK2_Gamma(1)*Id2 PK2_Gamma(3)*Id2 
         PK2_Gamma(3)*Id2 PK2_Gamma(2)*Id2];
ke_Gamma3=B_Gamma'*D*B_Gamma+G'*M_Gamma*G;
Int_Gamma3=B_Gamma'*D*E_Gamma;
dGamma3=Bl'*D*E_Gamma+B_Gamma'*D*0.5*Bl*Gamma(ele)*ue;

%% point4
s=-gauss_point;
t=gauss_point;
[~,B_Gamma,detJ,~,E_Gamma,~,G,Bl,~]=gamma_shape(s,t,ele,xx,yy,aa,bb,ue,nely,Gamma);
% Nonlinear  
PK2_Gamma=D*E_Gamma;
M_Gamma=[PK2_Gamma(1)*Id2 PK2_Gamma(3)*Id2 
         PK2_Gamma(3)*Id2 PK2_Gamma(2)*Id2];
ke_Gamma4=B_Gamma'*D*B_Gamma+G'*M_Gamma*G;
Int_Gamma4=B_Gamma'*D*E_Gamma;
dGamma4=Bl'*D*E_Gamma+B_Gamma'*D*0.5*Bl*Gamma(ele)*ue;

%% Integration
% Nonlinear
ke_Gamma=h*detJ*(ke_Gamma1+ke_Gamma2+ke_Gamma3+ke_Gamma4); 
Int_Gamma=h*detJ*(Int_Gamma1+Int_Gamma2+Int_Gamma3+Int_Gamma4);
dGamma=h*detJ*(dGamma1+dGamma2+dGamma3+dGamma4);

end