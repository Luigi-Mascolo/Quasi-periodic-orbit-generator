function ydot = eqmstab(t,yy,bd,cb,rkj,mu,hb,tspan)
%EQMSTAB integrates the 42 ODEs to solve the stability problem for quasi
%periodic orbits in the full ephemeris model (JPL's DE430 ephemeris)
%It uses a fixed step size integrator (e.g. RK5) to be able to retrieve
%ephemeris data from input and avoid calling SPICE/MICE in the loop.
%
%Requires SPICE/MICE
%
% USE:
%   
% INPUT:
%   yy      1x42 vector of initial conditions in the synodic reference
%           frame. These will be transposed. The first 36 elements
%           compose the STM while the last 6 are the position and the
%           velocity of the spacecraft (state) at t0
%   bd      1x10 boolean vector to include or neglect a celestial body
%           1   Mercury
%           2   Venus
%           3   Earth
%           4   Mars
%           5   Jupiter
%           6   Saturn
%           7   Uranus
%           8   Neptune
%           9   Sun
%           10  Moon
%   cb      central body (among bd)
%   mu      1xN vector including the gravitational parameters of the
%           included bodies. Please note mu values will be
%           nondimensionalized wrt the cb specific gravitational parameter
%   hb      Halo body (among bd). For EML2, hb = 10 (Moon); for SEL2, hb =
%           3 (Earth)
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator

x   = yy(37);
y   = yy(38);
z   = yy(39);
u   = yy(40);
v   = yy(41);
w   = yy(42);

% ZIMOVAN 2017

uii = @(mui,muk,x1,rn) (mui+muk).*((3.*(x1.^2))./(rn.^5)-1./(rn.^3));
uij = @(mui,muk,x1,x2,rn) (mui+muk).*((3.*x1.*x2)./(rn.^5));
ag = @(mm,xi,ri,xk,rk) mm.*(xi./(ri.^3)-xk./(rk.^3));

rki = [x y z];
rnki = norm(rki);

m = 0;
for k = 1:length(tspan)
    if abs(t-tspan(k))<1e-6, m = k; break; end
end
% fprintf('%.4f %.4f %d %d\n',t,tspan(k),k,m)
rkjt = reshape(rkj(:,m),[3 numel(bd)])';
rij = rkjt-rki;

A41 = uii(0,mu(cb),x,rnki);
A42 = uij(0,mu(cb),x,y,rnki);
A43 = uij(0,mu(cb),x,z,rnki);
A51 = A42;
A52 = uii(0,mu(cb),y,rnki);
A53 = uij(0,mu(cb),y,z,rnki);
A61 = A43;
A62 = A53;
A63 = uii(0,mu(cb),z,rnki);

for j = 1:length(bd)
    xij = rij(j,1);     % x posizione 
    yij = rij(j,2);
    zij = rij(j,3);
    rnij = norm(rij(j,:));
    rnij = sqrt(rij(j,1).^2+rij(j,2).^2+rij(j,3).^2);
    if bd(j) == 1 
        if j ~= cb
            A41 = A41+uii(0,mu(j),xij,rnij);
            A42 = A42+uij(0,mu(j),xij,yij,rnij);
            A43 = A43+uij(0,mu(j),xij,zij,rnij);
            A51 = A42;
            A52 = A52+uii(0,mu(j),yij,rnij);
            A53 = A53+uij(0,mu(j),yij,zij,rnij);
            A61 = A43;
            A62 = A53;
            A63 = A63+uii(0,mu(j),zij,rnij);
        end
    end
end

U = [A41 A42 A43;
    A51 A52 A53;
    A61 A62 A63];

O   = zeros(3);
I   = eye(3);

A  = [ O, I;
    U, O];
PHI_t = zeros(6,6);
% PHI_t = reshape(yy(1:36),[6,6]);

for i=1:6
    for j = 1:6
        PHI_t(i,j)=yy(6*(i-1)+j);
    end
end

DPHI_t = A*PHI_t;

temp = zeros(36,1);     % rendo vettore la matrice d(PHI(t))/dt con le ODEs
for i=1:6
    for j = 1:6
        temp(6*(i-1)+j)=DPHI_t(i,j);
    end
end    

ax = -(0+mu(cb))./(rnki.^3).*x;
ay = -(0+mu(cb))./(rnki.^3).*y;
az = -(0+mu(cb))./(rnki.^3).*z;
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

ydot = [u; v; w; ax; ay; az];

ydot = [temp; ydot];
end

