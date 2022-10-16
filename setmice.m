function [datafolder] = setmice()
%SETMICE includes the main information for NASA's MICE/SPICE setup for the
%current script
% 
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator

mr  = matlabroot;
for j=length(mr):-1:1
    if strcmp(mr(j),'\'), mr = mr(1:j-1); break; end
end

addpath(strcat(mr,'\mice\src\mice\'));
addpath(strcat(mr,'\mice\lib\'));
datafolder = strcat(mr,'\kernels');

% Load leapsecond kernel
cspice_furnsh( [datafolder, '\lsk\naif0012.tls.pc'] );
% Load (natural) satellite ephemeris kernel
cspice_furnsh( [datafolder, '\spk\satellites\jup310.bsp'] );
% Load planet ephemeris kernel
cspice_furnsh( [datafolder, '\spk\planets\de430.bsp'] );
% Load gravity contant kernel
cspice_furnsh( [datafolder, '\pck\gm_de431.tpc'] );
% Load planet constants kernel
cspice_furnsh( [datafolder, '\pck\pck00010.tpc'] );

end

