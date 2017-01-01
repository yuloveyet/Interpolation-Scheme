%%%% Geometrical Non-linear TOPOLOGY OPTIMIZATION CODE Nov, 2016 %%%%
clear;clc;
%% Input parameters
nelx=80;nely=20;%element numbers
a=1;b=0.25;h=0.1; % long beam size (m)
volfrac=0.2;penal=3; P=3;    

detaP=0.05; % contuation method
Beta=4;Eta=0.5; % Projection method
Beta_Int=500;threshold=0.01;% Interploation theme
T_Int=tanh(Beta_Int*threshold)+tanh(Beta_Int*(1-threshold));

ft=2;
rmin=nely/8; % filter length: numbers of elements surround
plane=1; % 1 plane strain 2 plane stress condition

Total_elem=nely*nelx;
Total_node=(nely+1)*(nelx+1);
Total_dofs=2*Total_node;

load=-5000e3;       %%%%%%%%%%%%% load (N) %%%%%%%%%%%%%%%%
steps=10;          %%%%%%%%% load step numbers %%%%%%%%%%%
stop=1e-5;         %%%%%%% convergence criterion %%%%%%%%%
%% Young Modulus (Pa)
E0 = 3e9;
Emin = 1e-9*E0;
nu = 0.4;
%% PREPARE FINITE ELEMENT ANALYSIS
% linear
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
ke_l = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% node coordinates
xx=zeros((nely+1)*(nelx+1),1);
yy=zeros((nely+1)*(nelx+1),1);
detaxx=linspace(0,a,nelx+1);  % x>0
detayy=-linspace(0,b,nely+1); % y<0
for j=1:nely+1
    for i=1:nelx+1
        xx((nely+1)*(i-1)+j)=detaxx(i);
        yy((nely+1)*(i-1)+j)=detayy(j)+b;         %% OFFSET y>0 %%
    end
end
%% element size
aa=detaxx(2);
bb=detayy(2);
%% DEFINE LOADS AND SUPPORTS (LONG BEAM:loaded downwards in the midpoint at the right edge)
F = sparse(2*(nely+1)*nelx+nely+2,1,load,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs =[1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
[dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
H = max(0,rmin-sqrt(dx.^2+dy.^2));
Hs = conv2(ones(nely,nelx),H,'same');
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
change = 1;
%% Iinitiation of MMA
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
m = 1;
n = nelx*nely;

onen    = ones(n,1);
onem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
%  a     = Column vector with the constants a_i in the terms a_i*z.
%  c     = Column vector with the constants c_i in the terms c_i*y_i.
%  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
a_mma   = zerom;
c_mma   = 1000*onem;
d_mma   = zerom;
a0      = 1;

xval    = x(:);
xold1   = xval;
xold2   = xold1;
%% START ITERATION
while change > 0.01
    tic% time start
    loop = loop + 1;
    xe=reshape(xPhys,nely*nelx,1);
    
    f=F/steps;% step load
    %%%%%%%%%%%%% Preallocation%%%%%%%%%%%%%
    ue=zeros(8,1);
    %     xe_bar=zeros(nely*nelx,1);
    %     Gamma=zeros(nely*nelx,1);
    ke_Gamma=cell(Total_elem,1);% cell
    fint_Gamma=cell(Total_elem,1);
    dGamma=cell(Total_elem,1);
    detaU=zeros(2*(nely+1)*(nelx+1),1);
    sK=zeros(64*Total_elem,1);
    sfint=zeros(8*Total_elem,1);
    ifint=zeros(8*Total_elem,1);
    %%%%%%%%%% Continuation scheme %%%%%%%%%%
%     if P<2
%         P=1+detaP*ceil(loop/2);
%     elseif P<3
%         P=2+detaP*ceil((loop-2/detaP)/5);
%     else
%         P=penal;
%     end
    if P<3
        Beta=4;
    elseif P==3 && Beta<64
        Beta=2*2^(ceil((loop-7/detaP+5)/10));
    else
        Beta=64;
    end
    
    T=tanh(Beta*Eta)+tanh(Beta*(1-Eta));
         %% Heaviside function 
            xe_bar=xe.^P-threshold;
            Gamma=(tanh(Beta_Int*threshold)+tanh(xe_bar.*Beta_Int))/T_Int;
    %% Geometrica Nonlinearity FE-ANALYSIS 
    for step=1:steps   % Incremental step
        Fext=f*step;% external nodal force
        Rmax=1;
        iter=0;
        while (Rmax>stop)
            iter=iter+1;
            % Start Parpool:
            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if isempty(poolobj)
                parpool
            end
            % Using Par for
            parfor ele=1:Total_elem%element number
                [ke_Gamma{ele},fint_Gamma{ele},dGamma{ele}]=Cal_Gamma(ele,edofMat,xx,yy,aa,bb,h,U,nely,nu,plane,Gamma);
            end

            for ele=1:Total_elem
                edof=edofMat(ele,:);% element dofs
                ifint((ele-1)*8+1:(ele)*8)=edof;
                %% Interpolation & SIMP
                sK((ele-1)*64+1:(ele)*64)=(Emin+xe(ele)^P*(E0-Emin))*(ke_l+Gamma(ele)^2*(ke_Gamma{ele}-ke_l*Gamma(ele)));
                sfint((ele-1)*8+1:(ele)*8)=(Emin+xe(ele)^P*(E0-Emin))*(ke_l*U(edof)+Gamma(ele)*(fint_Gamma{ele}-ke_l*U(edof)*Gamma(ele)^2));
            end
            
            K =sparse(iK,jK,sK); % global stiffness matrix
            Fint=sparse(ifint,ones(8*Total_elem,1),sfint);% global nodal force vector
            R=Fext-Fint;
            %% Newton Raphson Algorithm
            detaU(freedofs) =K(freedofs,freedofs)\R(freedofs);
            U=U+detaU;
            Rmax=max(abs(R(freedofs)));  %%%%%%%%% Boundary conditions %%%%%%%%%%%%%
             string=['subiter  ',int2str(iter), ' Rmax  ',sprintf('%12.8f',full(Rmax)),'  Rabs  ',sprintf('%12.8f', max(abs(full(R))))];
        disp(string)
        end
    end
    keyboard
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    c=F'*U;% end-compliance
    % Ajoint method
    L=zeros(2*(nely+1)*(nelx+1),1);
    L(freedofs)=K(freedofs,freedofs)\F(freedofs);
    dc=zeros(nely*nelx,1);
    dGammadxe=Beta_Int*sech(Beta_Int*(xe.^P-threshold)).^2*P.*xe.^(P-1)/T_Int;
    dfintdxe=zeros(8,1);
    dfintdGamma=zeros(8,1);
    dr=zeros(8,1);
    for ele=1:Total_elem
        edof=edofMat(ele,:);
        le=L(edof);
        dfintdxe=P*(E0-Emin)*xe(ele)^(P-1)*(ke_l*U(edof)+Gamma(ele)*(fint_Gamma{ele}-ke_l*U(edof)*Gamma(ele)^2));
        dfintdGamma=(Emin+xe(ele)^P*(E0-Emin))*(dGamma{ele}-3*Gamma(ele)^2*ke_l*U(edof));
        dr=-(dfintdxe+dfintdGamma.*dGammadxe);% chain rule
        dc(ele)=le'*dr;
    end
    
    dProjection=Beta*(sech(Beta*(xPhys-Eta))).^2/T;
    dc=reshape(dc,nely,nelx).*dProjection;
    dv = ones(nely,nelx).*dProjection;
    
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        %         dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
        dc = conv2(dc.*xPhys,H,'same')./Hs./max(1e-3,xPhys);
    elseif ft == 2
        dc = conv2(dc./Hs,H,'same');
        dv = conv2(dv./Hs,H,'same');
    end
    
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    
    %bounds for design variables
    move = 0.1;
    % column vector
    xmin = max(xval-move, zeron);
    xmax = min(xval+move, onen);
    low     = xmin;
    upp     = xmax;
    
    f0val = log(c);
    df0dx = dc(:)/c;
    df0dx2 = 0*df0dx;
    fval = (sum(xPhys(:)))/ volfrac/(nelx*nely)- 1;  % column vector
    dfdx = dv(:)'/ volfrac/(nelx*nely) ;    % (m * n)
    dfdx2 = 0*dfdx;
    
    %%%% The MMA subproblem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
        mmasub(m,n,loop,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a_mma,c_mma,d_mma);
    %%%% Some vectors are updated:
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
    xnew = reshape(xval, nely, nelx);
    
    if ft == 1
        xPhys = xnew;
    elseif ft == 2
%         x_tilde= conv2(xnew,H,'same')./Hs;
        xPhys= conv2(xnew,H,'same')./Hs;
        %% Projection filter
        xPhys=(tanh(Beta*Eta)+tanh((xPhys-Eta)*Beta))/T;
    end
        
    change = max(abs(xnew(:)-x(:)));
    x = xnew;  
    %% PLOT DENSITIES
%     [map] =brewermap(9,'*Blues');
%     colormap(map); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
 %% plot disp
  clf;
  [map] =brewermap(9,'*Blues');
  colormap(map);
  axis equal;axis tight; axis off
  for ely=1:nely
      for elx=1:nelx
          edofs=edofMat(nely*(elx-1)+ely,:);
          Ue=U(edofs);
          ly=(nely-(ely-1))*aa;lx=(elx-1)*aa;
          xxx=[Ue(1,1)+lx Ue(3,1)+lx+aa Ue(5,1)+lx+aa Ue(7,1)+lx];
          yyy=[Ue(2,1)+ly Ue(4,1)+ly Ue(6,1)+ly+aa Ue(8,1)+ly+aa];
%           aveUe=sqrt(sum(xxx.^2+yyy.^2))/4;
%           patch(xxx,yyy,-aveUe)
          patch(xxx,yyy,-xPhys(ely,elx))
          %       patch(xxx,yyy,-xe(nely*(elx-1)+ely))
      end
  end
  drawnow;
  t=toc;% time finish
    %% PRINT RESULTS
    fprintf('It.:%4i Obj.:%12.4f Vol.:%1.3f ch.:%1.3f sec.:%4.3f Rmax.:%3.3f p.:%1.2f B.:%2i\n',loop,c, ...
        mean(xPhys(:)),change,t,full(Rmax)/stop,P,Beta);
 

end


