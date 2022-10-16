function [jc,V2,U2] = jacobiValue3D(X0,mu)
%JACOBIVALUE3D computes the jacobi value of a state in the CR3BP
%
%   INPUT:
%       X  - 1x6 Vector - Your current states [x,y,z,xDot,yDot,zDot]
%       mu - double     - Your systems nondimenisonal mass ratio
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator

x   = X0(1);
y   = X0(2);
z   = X0(3);
vx  = X0(4);
vy  = X0(5);
vz  = X0(6);
r   = sqrt(((x-1+mu).^2)+(y.^2)+(z.^2));
d   = sqrt(((x+mu).^2)+(y.^2)+(z.^2));
V2  = (vx.^2)+(vy.^2)+(vz.^2);
U2  = (x.^2)+(y.^2) +(2*(1-mu)./d)+(2*mu./r);

jc = (x.^2)+(y.^2) +(2*(1-mu)./d)+(2*mu./r)-((vx.^2)+(vy.^2)+(vz.^2));
end

