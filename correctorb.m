function [X0NEW,XFNEW,XX0,XXF,SS2,XXB0,XXBF,errori] = correctorb(X0,XF,XB0,XBF,STM,varyX0,checkXF,desXF,ril)
%CORRECTORB performes the differential correction to compute the QPOs in
%the ephemerides model
%
%   INPUT:
%       X0  SC initial state in the synodic RF
%       XF  SC final state in the synodic RF
%       XB0 non-central primary initial state in the synodic RF
%       XBF non-central primary final state in the synodic RF
%       STM complete State-Transition-Matrix
%       varyx0  string vector containing information about initial
%               variables allowed to change during differential correction
%       checkXF string vector containing information about terminal
%               constraints to be fulfilled during differential correction
%       desXF   vector containing the numerical values of constraints
%               included in checkXF (usually zeros)
%       ril relaxation parameter (0<ril<1)
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator
dx0v = zeros(1,7);
xfv = zeros(1,6);
if isempty(desXF), desXF = xfv; end
errori = zeros(1,6);
variables = {'x','y','z','vx','vy','vz','t'};
for k = 1:length(variables)
    dx0v(k) = 1;
    for j = 1:length(varyX0)
        if strcmp(variables{k},varyX0{j})
            dx0v(k) = 0;
        end
    end
end
for k = 1:length(variables)
    xfv(k) = 0;
    for j = 1:length(checkXF)
        if strcmp(variables{k},checkXF{j})
            xfv(k) = 1;
        end
    end
end

SS = [];
XX0 = [];
XXF = [];
XXB0 = [];
XXBF = [];
c = 0;
X0NEW = X0;
XFNEW = XF;
sav = [];
for k = 1:length(dx0v)
    if dx0v(k) == 0
        c = c+1;
        SS(:,c) = STM(:,k);
        XX0(c) = X0(k);
        if k<=6
            XXB0(c) = XB0(k);
        end
        sav(c) = k;
    end
end
c = 0;
SS2 = [];
sav2 = [];
for k = 1:length(xfv)
    if xfv(k) == 1
        c = c+1;
        SS2(c,:) = SS(k,:);
        XXF(c) = XF(k)-desXF(c);
        XXBF(c) = XBF(k);
        sav2(c) = k;
    end
end
if size(XXF,1) == 1, XXF = XXF'; end
if size(XX0,1) == 1, XX0 = XX0'; end
if size(XXBF,1) == 1, XXBF = XXBF'; end
if size(XXB0,1) == 1, XXB0 = XXB0'; end

CORR = SS2.'*(inv(SS2*SS2.'))*XXF*ril;

X0NEW(sav) = XX0-CORR;
XFNEW(sav2) = XXF;
errori(sav2) = XXF;

end