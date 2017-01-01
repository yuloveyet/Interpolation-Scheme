%% MATERIAL PROPERTIES
function [D]=C(nu,plane) 

if plane==1 % plane strain
    d=1/((1+nu)*(1-2*nu));
    D=d*[1-nu  nu    0
        nu  1-nu   0
        0    0  (1-2*nu)/2];
elseif plane==2 % plane stress
    d=1/(1-nu*nu); % unit elastic matrix
    D=d*[ 1  nu    0
        nu  1    0
        0  0 (1-nu)/2];
end

end