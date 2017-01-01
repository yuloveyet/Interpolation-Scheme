%%%% Geometrical Non-linear TOPOLOGY OPTIMIZATION CODE in DTU Dec, 2016 %%%%
clear;clc;
% Start Parpool:
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool
end
movelimit=0.05;
%% Input parameters
nelx=80;nely=20;%element numbers
a=1;b=0.25;h=0.1; % long beam size (m)
vol=0;volfrac=0.5;penal=3; P=1;

detaP=0.1; % contuation method
Beta=0;Eta=0.5; % Projection method
Beta_Int=500;threshold=0.01;% Interploation theme
T_Int=tanh(Beta_Int*threshold)+tanh(Beta_Int*(1-threshold));

ft=2;% 1-sensitivity filter 2-density filter
rmin=nely/8; % filter length: numbers of elements surround
plane=1; % 1-plane strain 2-plane stress condition

Total_elem=nely*nelx;
Total_node=(nely+1)*(nelx+1);
Total_dofs=2*Total_node;

force=-300e3;       %%%%%%%%%%%%% load (N) %%%%%%%%%%%%%%%%
steps=10;          %%%%%%%%% load step numbers %%%%%%%%%%%
stop=1e-4;         %%%%%%% convergence criterion %%%%%%%%%
%% Young Modulus (Pa)
E0 = 3e9;
Emin = 1e-9*E0;
nu = 0.4;
%% NON-DESIGN DOMAIN (Void&Solid)
as=nelx/5;bs=nely/2;
vols=as*bs;
nondesign = zeros(nely,nelx);
% SQUARE DOMAIN
for i = 4*nelx/10+1:4*nelx/10+as
    for j = nely/4+1:nely/4+bs
        nondesign(j,i)=1;
    end
end
E1=E0;
%% PREPARE FINITE ELEMENT ANALYSIS
% linear
if plane==2 % plane stress
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
    A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
    %% h
    ke_l =h* 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
    
elseif plane ==1 % plane strain
    ke_l =h*[    7/(18*(nu+1))-1/(9*(2*nu-1)),     -1/(8*(2*nu-1)*(nu+1)),    1/(9*(2*nu-1))-5/(36*(nu+1)),-5/(24*(nu+1))-1/(12*(2*nu-1)),   1/(18*(2*nu-1))-7/(36*(nu+1)),      1/(8*(2*nu-1)*(nu+1)),-1/(18*(nu+1))-1/(18*(2*nu-1)),   5/(24*(nu+1))+1/(12*(2*nu-1))
        -1/(8*(2*nu-1)*(nu+1)),    7/(18*(nu+1))-1/(9*(2*nu-1)),   5/(24*(nu+1))+1/(12*(2*nu-1)),-1/(18*(nu+1))-1/(18*(2*nu-1)),      1/(8*(2*nu-1)*(nu+1)),   1/(18*(2*nu-1))-7/(36*(nu+1)),-5/(24*(nu+1))-1/(12*(2*nu-1)),    1/(9*(2*nu-1))-5/(36*(nu+1))
        1/(9*(2*nu-1))-5/(36*(nu+1)),   5/(24*(nu+1))+1/(12*(2*nu-1)),    7/(18*(nu+1))-1/(9*(2*nu-1)),      1/(8*(2*nu-1)*(nu+1)),-1/(18*(nu+1))-1/(18*(2*nu-1)),-5/(24*(nu+1))-1/(12*(2*nu-1)),   1/(18*(2*nu-1))-7/(36*(nu+1)),     -1/(8*(2*nu-1)*(nu+1))
        -5/(24*(nu+1))-1/(12*(2*nu-1)),-1/(18*(nu+1))-1/(18*(2*nu-1)),      1/(8*(2*nu-1)*(nu+1)),    7/(18*(nu+1))-1/(9*(2*nu-1)),   5/(24*(nu+1))+1/(12*(2*nu-1)),    1/(9*(2*nu-1))-5/(36*(nu+1)),     -1/(8*(2*nu-1)*(nu+1)),   1/(18*(2*nu-1))-7/(36*(nu+1))
        1/(18*(2*nu-1))-7/(36*(nu+1)),      1/(8*(2*nu-1)*(nu+1)),-1/(18*(nu+1))-1/(18*(2*nu-1)),   5/(24*(nu+1))+1/(12*(2*nu-1)),    7/(18*(nu+1))-1/(9*(2*nu-1)),     -1/(8*(2*nu-1)*(nu+1)),    1/(9*(2*nu-1))-5/(36*(nu+1)),-5/(24*(nu+1))-1/(12*(2*nu-1))
        1/(8*(2*nu-1)*(nu+1)),   1/(18*(2*nu-1))-7/(36*(nu+1)),-5/(24*(nu+1))-1/(12*(2*nu-1)),    1/(9*(2*nu-1))-5/(36*(nu+1)),     -1/(8*(2*nu-1)*(nu+1)),    7/(18*(nu+1))-1/(9*(2*nu-1)),   5/(24*(nu+1))+1/(12*(2*nu-1)),-1/(18*(nu+1))-1/(18*(2*nu-1))
        -1/(18*(nu+1))-1/(18*(2*nu-1)),-5/(24*(nu+1))-1/(12*(2*nu-1)),   1/(18*(2*nu-1))-7/(36*(nu+1)),     -1/(8*(2*nu-1)*(nu+1)),    1/(9*(2*nu-1))-5/(36*(nu+1)),   5/(24*(nu+1))+1/(12*(2*nu-1)),    7/(18*(nu+1))-1/(9*(2*nu-1)),      1/(8*(2*nu-1)*(nu+1))
        5/(24*(nu+1))+1/(12*(2*nu-1)),    1/(9*(2*nu-1))-5/(36*(nu+1)),     -1/(8*(2*nu-1)*(nu+1)),   1/(18*(2*nu-1))-7/(36*(nu+1)),-5/(24*(nu+1))-1/(12*(2*nu-1)),-1/(18*(nu+1))-1/(18*(2*nu-1)),      1/(8*(2*nu-1)*(nu+1)),    7/(18*(nu+1))-1/(9*(2*nu-1))];
end
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
        %          yy((nely+1)*(i-1)+j)=detayy(j);
        yy((nely+1)*(i-1)+j)=detayy(j)+b;         %% OFFSET y>0 %%
    end
end
%% element size
aa=detaxx(2);
bb=detayy(2);
%% DEFINE LOADS AND SUPPORTS (LONG BEAM:loaded downwards in the midpoint at the right edge)
F = sparse(2*(nely+1)*nelx+nely+2,1,force,2*(nely+1)*(nelx+1),1);
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
Tloop=0;
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
while change > 1e-3 || Beta<64
    tic% time start
    %% Preallocation
    ue=zeros(8,1);
    ke_Gamma=cell(Total_elem,1);% cell form
    fint_Gamma=cell(Total_elem,1);
    detaU=zeros(2*(nely+1)*(nelx+1),1);
    sK=zeros(64*Total_elem,1);
    sfint=zeros(8*Total_elem,1);
    ifint=zeros(8*Total_elem,1);
    %% Continuation method
    if P<2
        P=1+detaP*floor(Tloop/2);% P0=1
    elseif P<3
        P=2+detaP*floor((Tloop-2/detaP)/5);
    else
        P=penal;
    end
    Tloop= Tloop+1;
    
    if  P <3
        Beta=4;
    elseif (mod(Tloop,10)==0 || change < 1e-3) && Beta < 64 && P>=3
        Beta=Beta*2;
        if Beta > 64
            Beta=64;
        end
        loop=0;
        change=0.1;
    end
    loop = loop + 1;
    
    T=tanh(Beta*Eta)+tanh(Beta*(1-Eta));
    xe=reshape(xPhys,nely*nelx,1);
    f=F/steps;% step load
    %% Heaviside function
    xe_bar=xe.^P-threshold;
    Gamma=(tanh(Beta_Int*threshold)+tanh(xe_bar.*Beta_Int))/T_Int;
    %% SIMP
    SIMP=Emin+xe.^P*(E0-Emin);
    %% Geometrica Nonlinearity FE-ANALYSIS
    for step=1:steps   % Incremental step
        Fext=f*step;% external nodal force
        Rmax=1;
        while (Rmax>stop)
            
            % Using Par for
            parfor ele=1:Total_elem%element number
                [ke_Gamma{ele},fint_Gamma{ele}]=Cal_Gamma(ele,edofMat,xx,yy,aa,bb,h,U,nely,nu,plane,Gamma);
            end
            
            for ele=1:Total_elem
                edof=edofMat(ele,:);% element dofs
                ifint((ele-1)*8+1:(ele)*8)=edof;
                %% Interpolation
                sK((ele-1)*64+1:(ele)*64)=SIMP(ele)*(ke_l+Gamma(ele)^2*ke_Gamma{ele}-Gamma(ele)^2*ke_l);
                sfint((ele-1)*8+1:(ele)*8)=SIMP(ele)*(ke_l*U(edof)+Gamma(ele)*fint_Gamma{ele}-Gamma(ele)^2*ke_l*U(edof));
            end
            K =sparse(iK,jK,sK); % global stiffness matrix
            Fint=sparse(ifint,ones(8*Total_elem,1),sfint);% global nodal force vector
            R=Fext-Fint;
            %% Newton Raphson Algorithm
            detaU(freedofs) =K(freedofs,freedofs)\R(freedofs);
            U=U+detaU;
            Rmax=max(abs(R(freedofs)));  %%%%%%%%% Boundary conditions %%%%%%%%%%%%%
        end
    end
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    c=F'*U;% end-compliance
    L=zeros(2*(nely+1)*(nelx+1),1);
    L(freedofs)=K(freedofs,freedofs)\F(freedofs);% Ajoint method
    
    dc=zeros(nely*nelx,1);
    dfintdxe=zeros(8,1);
    dfintdGamma=zeros(8,1);
    dr=zeros(8,1);
    
    for ele=1:Total_elem
        edof=edofMat(ele,:);
        le=L(edof);
        dGammadxe=Beta_Int*sech(Beta_Int*(xe(ele)^P-threshold))^2*P*xe(ele)^(P-1)/T_Int;
        dfintdxe=P*(E0-Emin)*xe(ele)^(P-1)*(ke_l*U(edof)+Gamma(ele)*fint_Gamma{ele}-ke_l*U(edof)*Gamma(ele)^2);
        dfintdGamma=SIMP(ele)*(fint_Gamma{ele}+Gamma(ele)*ke_Gamma{ele}*U(edof)-2*Gamma(ele)*ke_l*U(edof));
        dr=-(dfintdxe+dfintdGamma*dGammadxe);% chain rule
        dc(ele)=le'*dr;
    end
    
    dProjection=Beta*(sech(Beta*(xPhys-Eta))).^2/T;
    dc=reshape(dc,nely,nelx).*dProjection;
    dv = ones(nely,nelx).*dProjection;
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc = conv2(dc.*xPhys,H,'same')./Hs./max(1e-3,xPhys);
    elseif ft == 2
        dc = conv2(dc./Hs,H,'same');
        dv = conv2(dv./Hs,H,'same');
    end
    % SQUARE DOMAIN
    dc(nondesign==1)=0;
    dv(nondesign==1)=0;
%     dc_sp(nondesign==1)=0;
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    
    %bounds for design variables
    move = movelimit;
    % column vector
    xmin = max(xval-move, zeron);
    xmax = min(xval+move, onen);
    low     = xmin;
    upp     = xmax;
    
    f0val = log(c);
    df0dx = dc(:)/c;
    df0dx2 = 0*df0dx;
    fval = (sum(xPhys(:))-vols)/ volfrac/(nelx*nely-vols)- 1;  % column vector
    dfdx = dv(:)'/ volfrac/(nelx*nely-vols) ;    % (m * n)
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
        x_tilde= conv2(xnew,H,'same')./Hs;
        %% Projection filter
        xPhys=(tanh(Beta*Eta)+tanh((x_tilde-Eta)*Beta))/T;
    end
    
     %% NONDESIGN
  xPhys(nondesign==1)=1;
  
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    vol=sum(xPhys(:));
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
    fprintf('It.:%4i l.:%4i Obj.:%2.4f Vol.:%1.3f ch.:%1.3f sec.:%4.3f Rmax.:%3.3f p.:%1.2f B.:%2i\n',Tloop,loop,log(c), ...
        (vol-vols)/nelx/nely,change,t,full(Rmax)/stop,P,Beta);
%     %% Create a folder
%     Name=['TopOpt'];fname='';
%     %% write data in a file
%     aa=strcat(Name,'/Results_',int2str(fname),'.txt');
%     fdes = fopen(aa, 'wt');
%     aa=[' It.: ' sprintf('%4i',Tloop) ...
%         ' Obj.: ' sprintf('%9.5f',c)...
%         %         ' Gray: ' sprintf('%9.6f  ',GraL) ...
%         ' Vol: ' sprintf('%9.6f  ',vol) ...
%         '  beta : ' sprintf('%9.6f',Beta) ...
%         %         '  Fac: ' sprintf('%9.6f',Fac) ...
%         ' p: ' sprintf('%9.6f',P) ...
%         %         '  fval: ' sprintf('%9.6f',fval(end)) ...
%         ' change: ' sprintf('%9.6f',change) ];
%     disp(aa);
%     fprintf(fdes, aa);
%     fprintf(fdes,'\n');
%     fclose(fname);
%     %%  save figure
%     set(gcf,'paperpositionmode','auto');
%     aa=strcat(Name,'/Design_MS');
% %     print(gcf,'-r300',aa,'-depsc');
% %     print(gcf,'-r300',aa,'-dpdf');
%     print(gcf,'-r300',aa,'-dpng');

end


