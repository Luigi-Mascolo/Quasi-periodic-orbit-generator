function Xout = rot2j2000(Xin,r12v,mu,rconvlag,vconvlag,sys)
%ROT2J2000 converts an initial state X0 in a rotating synodic frame to a
%J2000 inertial frame defined by the instantaneous ephemeris of a two-body
%primary system.
%
%   INPUT:
%       Xin     Initial state in rotating RF
%       r12v    vector containing relative positioning between primary 2
%               and primary 1 gravitational bodies
%       v12v    vector containing relative velocities between primary 2
%               and primary 1 gravitational bodies
%       mu      CR3BP specific mass ration (smaller primary over total)
%       rconvlag, vconvlag  position and velocity conversion factors
%       sys     flag identifying the binary system (Earth-Moon / Sun-Earth)
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator


%1. Translation from barycenter to primary (1)
if strcmp(sys,'EM')
    Xin(1,:) = Xin(1,:)+mu;
elseif strcmp(sys,'SE')
    Xin(1,:) = Xin(1,:)-(1-mu);
end

%2. Dimensionalize spacecraft state
Xin(1:3,:) = Xin(1:3,:).*rconvlag;
Xin(4:6,:) = Xin(4:6,:).*vconvlag;

%3. Versors computation and rotation
Xout = zeros(size(Xin));
for k = 1:size(Xout,2)
    r12 = r12v(1:3,k);
    v12 = r12v(4:6,k);
    xv = r12/(norm(r12));
    h = cross(r12,v12);
    zv = h/(norm(h));
    yv = cross(zv,xv);
    if strcmp(sys,'SE')
        xv = -xv;
        yv = -yv;
        zv = -zv;
    end
    om = norm(h)/((norm(r12))^2);
    C = [xv yv zv];
    O = zeros(3,3);
    CD = om.*[C(1,2) -C(1,1) 0;
        C(2,2) -C(2,1) 0;
        C(3,2) -C(3,1) 0];
    CRI = [C O;
        CD C];
    Xout(:,k) = CRI*Xin(:,k);
end

%4. Non-dimensionalization of spacecraft state in J2000

Xout(1:3,:) = Xout(1:3,:)./rconvlag;
Xout(4:6,:) = Xout(4:6,:)./vconvlag;

end
