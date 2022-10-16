function [X,x,y,z,vx,vy,vz] = zentoijk(deg,r,th,ph,u,v,w)
%ZENTOIJK converts the spacecraft's state from ZEN to cartesian IJK
%coordinates
%
%   INPUT:
%       deg     flag to identify if operating in degrees or radians
%       r       SC position vector
%       th      SC right ascension
%       ph      SC declination
%       u       SC velocity in the radial direction
%       v       SC velocity in the tangential direction
%       w       SC velocity in the normal direction
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator

if (size(r,1)>1 && size(r,2)>1)||(size(r,1)==1&&size(r,2)==6)||(size(r,1)==6&&size(r,2)==1)
    if size(r,1)==6, r = r'; end
    th = r(:,2);
    ph = r(:,3);
    u = r(:,4);
    v = r(:,5);
    w = r(:,6);
    r = r(:,1);
end 
if deg == 2, th = th.*180./pi; ph = ph.*180./pi; end
x = r.*cosd(th).*cosd(ph);
y = r.*sind(th).*cosd(ph);
z = r.*sind(ph);

vx = u.*cosd(th).*cosd(ph)+v.*(-sind(th))+w.*(-cosd(th).*sind(ph));
vy = u.*sind(th).*cosd(ph)+v.*cosd(th)+w.*(-sind(th).*sind(ph));
vz = u.*sind(ph)+w.*cosd(ph);
X = [x y z vx vy vz]';
end

