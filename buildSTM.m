function [STM,x,y,z,vx,vy,vz] = buildSTM(yout,rkj,bd,cb,mu)
%BUILDSTM returns the State Transition Matrix for a specific state in the
%ephemerides model
%
%   INPUT:
%   	yout    1x42 column vector containing 36 final values for the
%               monodromy matrix and 6 final values for the spacecraft
%               state (3 positions, 3 velocities)
%       rkj     matrix containing the relative positions of all celestial
%               bodies with respect to the central one
%       bd      boolean vector for the celestial bodies (1 present, 0
%               neglected)
%       cb      id identifying the central body (among bd)
%       mu      CR3BP specific mass parameter (smaller primary over total)
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator
temp = yout(end,1:36);
STM = reshape(temp,[6,6]);
x = yout(:,37);
y = yout(:,38);
z = yout(:,39);
vx = yout(:,40);
vy = yout(:,41);
vz = yout(:,42);

xf = yout(end,37);
yf = yout(end,38);
zf = yout(end,39);
vxf = yout(end,40);
vyf = yout(end,41);
vzf = yout(end,42);

uii = @(mui,muk,x1,rn) (mui+muk).*((3.*(x1.^2))./(rn.^5)-1./(rn.^3));
uij = @(mui,muk,x1,x2,rn) (mui+muk).*((3.*x1.*x2)./(rn.^5));
ag = @(mm,xi,ri,xk,rk) mm.*(xi./(ri.^3)-xk./(rk.^3));

rki = [xf yf zf];
rnki = norm(rki);
rkjt = reshape(rkj(:,end),[3 numel(bd)])';
rij = rkjt-rki;

ax = -(0+mu(cb))./(rnki.^3).*xf;
ay = -(0+mu(cb))./(rnki.^3).*yf;
az = -(0+mu(cb))./(rnki.^3).*zf;
for j = 1:length(bd)
    xij = rij(j,1);     % x posizione
    yij = rij(j,2);
    zij = rij(j,3);
    xkj = rkjt(j,1);
    ykj = rkjt(j,2);
    zkj = rkjt(j,3);
    rnij = norm(rij(j,:));
    rnkj = norm(rkjt(j,:));
    if bd(j) == 1
        if j~= cb
            ax = ax+ag(mu(j),xij,rnij,xkj,rnkj);
            ay = ay+ag(mu(j),yij,rnij,ykj,rnkj);
            az = az+ag(mu(j),zij,rnij,zkj,rnkj);
        end
    end
end

STM = [STM(1,:) vxf;
    STM(2,:) vyf;
    STM(3,:) vzf;
    STM(4,:) ax(end);
    STM(5,:) ay(end);
    STM(6,:) az(end)];
end

