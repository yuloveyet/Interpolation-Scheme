%% element stiffness matrix and internal nodal force
function [ke_Gamma,fint_Gamma]=Cal_Gamma(ele,edofMat,xx,yy,aa,bb,h,U,nely,nu,plane,Gamma)
edofs=edofMat(ele,:);% element dofs
%% Gamma ue
ue=U(edofs);
%% element nodes coordinates
xxyy=zeros(4,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%
node=ele+ceil(ele/nely)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%
xxyy(4,1)=xx(node);
xxyy(3,1)=xxyy(4,1)+aa;
xxyy(2,1)=xxyy(3,1);
xxyy(1,1)=xxyy(4,1);
xxyy(4,2)=yy(node);
xxyy(3,2)=xxyy(4,2);
xxyy(2,2)=xxyy(3,2)+bb;
xxyy(1,2)=xxyy(2,2);

[D]=C(nu,plane);% unit elastic matrix
[ke_Gamma,fint_Gamma]=knonlinear_Gamma(ele,xxyy,D,h,ue,Gamma);
end