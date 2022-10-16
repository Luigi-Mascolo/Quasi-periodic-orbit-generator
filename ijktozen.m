function [r,th,ph,u,v,w] = ijktozen(x,y,z,vx,vy,vz)
%IJKTOZEN converts the spacecraft's state from cartesian IJK to ZEN
%coordinates
%
%   INPUT:
%       deg     flag to identify if operating in degrees or radians
%       x, y, z SC cartesian position components
%       vx, vy, vz  SC cartesian velocity components
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator

if (size(x,1)>1 && size(x,2)>1)||(size(x,1)==1&&size(x,2)==6)||(size(x,1)==6&&size(x,2)==1)
    if size(x,1)==6, x = x'; end
    y = x(:,2);
    z = x(:,3);
    vx = x(:,4);
    vy = x(:,5);
    vz = x(:,6);
    x = x(:,1);
end 

r = sqrt(x.^2+y.^2+z.^2);
th = atan2d(y,x);
ph = asind(z./r);
u = vx.*cosd(th).*cosd(ph)+vy.*sind(th).*cosd(ph)+vz.*sind(ph);
v = -vx.*sind(th)+vy.*cosd(th);
w = -(vx.*cosd(th).*sind(ph)+vy.*sind(th).*sind(ph)-vz.*cosd(ph));
end

