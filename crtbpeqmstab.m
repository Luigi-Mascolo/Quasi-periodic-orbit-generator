function ydot = crtbpeqmstab(t,yy,mu)
%CRTBPEQMSTAB integrates the 42 ODEs composing the correction problem for
%periodic orbits in the CR3BP
%   INPUT:
%       yy      42 initial values (36 STMP, 6 spacecraft state)
%       mu      CR3BP specific mass parameter (smaller primary over total)
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

r1      = sqrt((x+mu).^2+y.^2+z.^2);     % m1-m3 distance
r2      = sqrt((x+-1+mu).^2+y.^2+z.^2);   % m2-m3 distance

Uxx = 1-(1-mu)/(r1^5)*(r1^2-3*(x+mu)^2)-mu/(r2^5)*(r2^2-3*(x-1+mu)^2);
Uxy = 3*y*((1-mu)/(r1^5)*(x+mu)+mu/(r2^5)*(x-1+mu));
Uxz = 3*z*((1-mu)/(r1^5)*(x+mu)+mu/(r2^5)*(x-1+mu));
Uyx = Uxy;
Uyy = 1-(1-mu)/(r1^5)*(r1^2-3*y^2)-mu/(r2^5)*(r2^2-3*y^2);
Uyz = 3*y*z*((1-mu)/(r1^5)+mu/(r2^5));
Uzx = Uxz;
Uzy = Uyz;
Uzz = -(1-mu)/(r1^5)*(r1^2-3*z^2)-mu/(r2^5)*(r2^2-3*z^2);

U = [Uxx Uxy Uxz;
    Uyx Uyy Uyz;
    Uzx Uzy Uzz];

O   = zeros(3);
I   = eye(3);
K   = [0, 2, 0;
    -2, 0, 0;
    0, 0, 0];

A  = [ O, I;
    U, K];
PHI_t = zeros(6,6);

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

ax = 2*v+x-(1-mu)*(x+mu)/(r1^3)-mu*(x-1+mu)/(r2^3);
ay = -2*u+y-(1-mu)*y/(r1^3)-mu*(y/r2^3);
az = -(1-mu)*z/(r1^3)-mu*z/(r2^3);

ydot = [u; v; w; ax; ay; az];

ydot = [temp; ydot];

end

