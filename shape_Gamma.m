%% shape function (Iso-parametric Q-4 element)
function [B,B_Gamma,detJ,ele_strain,E_Gamma,C_Gamma,G,Bl,B0]=shape_Gamma(s,t,ele,xxyy,ue,Gamma)
 B0=zeros(3,8);
 Bl=zeros(3,8);Bl_Gamma=zeros(3,8);
 B=zeros(3,8);B_Gamma=zeros(3,8);
 G=zeros(4,8);
 A=zeros(3,4);A_Gamma=zeros(3,4);
 ele_strain=zeros(3,1);E_Gamma=zeros(3,1);C_Gamma=zeros(3,1);
 Theta=zeros(4,1);Theta_Gamma=zeros(4,1);
%% Interpolation function 
% N=[(1-s)*(1-t),(1+s)*(1-t),(1+s)*(1+t),(1-s)*(1+t)]/4; 
dNds=[t-1,1-t,1+t,-1-t]/4;
dNdt=[s-1,-1-s,1+s,1-s]/4;
%% Jacobi matrix
% J=zeros(2,2);
% J=[dNds;dNdt]*xxyy;
% detJ=det(J);
% inversJ=inv(J)*[dNds;dNdt];
J11=dNds*xxyy(:,1);
J12=dNds*xxyy(:,2);
J21=dNdt*xxyy(:,1);
J22=dNdt*xxyy(:,2);
% J=[J11 J12
%    J21 J22];
detJ=J11*J22-J12*J21;
% inv(J)=[J22 -J12
%         -J21 J11]/detJ;
% [dNdx;dNdy]=inv(J)*[dNds;dNdt];

dNdx=(J22*dNds-J12*dNdt)/detJ;
dNdy=(J11*dNdt-J21*dNds)/detJ;

B0=[dNdx(1)   0      dNdx(2)   0      dNdx(3)   0      dNdx(4)   0
       0    dNdy(1)    0     dNdy(2)    0     dNdy(3)    0     dNdy(4) 
    dNdy(1) dNdx(1)  dNdy(2) dNdx(2)  dNdy(3) dNdx(3)  dNdy(4) dNdx(4)];

G=[dNdx(1)   0      dNdx(2)   0      dNdx(3)   0      dNdx(4)   0
     0    dNdx(1)     0    dNdx(2)    0     dNdx(3)    0     dNdx(4)
   dNdy(1)   0      dNdy(2)   0      dNdy(3)   0      dNdy(4)   0 
     0    dNdy(1)    0     dNdy(2)    0     dNdy(3)    0     dNdy(4)];

Theta=G*ue;
Thetax=[Theta(1);Theta(2)];
Thetay=[Theta(3);Theta(4)];
A=[Thetax' 0 0;0 0 Thetay';Thetay' Thetax'];
Bl=A*G;
B=B0+Bl;

ele_strain=(B0+0.5*Bl)*ue;

Theta_Gamma=G*Gamma(ele)*ue;
Thetax_Gamma=[Theta_Gamma(1);Theta_Gamma(2)];
Thetay_Gamma=[Theta_Gamma(3);Theta_Gamma(4)];
A_Gamma=[Thetax_Gamma' 0 0;0 0 Thetay_Gamma';Thetay_Gamma' Thetax_Gamma'];
Bl_Gamma=A_Gamma*G;
B_Gamma=B0+Bl_Gamma;

E_Gamma=(B0+0.5*Bl_Gamma)*ue*Gamma(ele);
% C_Gamma=2*(B0+0.5*Bl_Gamma)*ue*Gamma(ele)+[1;1;0];
end