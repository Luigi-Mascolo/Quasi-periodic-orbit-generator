function [rkj,vkj,tdim,tadim,tfadim] = geteph(mode,et0,tf,str,ndim,IDS,bd,cb,NE,t00,conv)
%GETEPH retrieves the ephemerides for specific celestial bodies (state,
%i.e. positions and velocities, over time)
%
%   INPUT:
%       mode    select if eph. are retrieved from JPL SPICE/MICE (1) or
%               JPLEPH (fortran)
%       et0     starting epoch (elapsed time since J2000)
%       tf      duration
%       str     reference frame ID
%       ndim    dimensionalization parameter
%       IDS     gravitational body identifiers
%       bd      vector containing boolean values for celestial bodies
%       cb      identifier for the central body
%       NE      number of elements
%       t00     (fortran) identifier for t0
%       conv    struct containing conversion parameters
%
%   INTERNAL CALLS
%       zentoijk    ZEN to IJK RF conversion
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator

tconvelio = conv.te;
rconv = conv.rt;
vconv = conv.vt;


if mode == 1
    tf = tf*ndim*86400;
    t = linspace(et0,et0+tf,NE);
    for k = 1:length(bd)
        if bd(k) == 1
            X = cspice_spkezr(IDS{k},t,str,'none',IDS{cb});
        else
            X = zeros(6,length(t));
        end
        rkj(3*(k-1)+1:3*k,:) = X(1:3,:);
        vkj(3*(k-1)+1:3*k,:) = X(4:6,:);
    end
    tdim = t;
    tf = tf/ndim/86400;
    tadim = linspace(0,tf,NE);
    tfadim = tadim(end);
else
    t00 = t00*tconvelio*86400;
    t0f = t00+tf*ndim*86400;
    tot = linspace(t00,t0f,NE);
    for k = 1:length(bd)
        fp = fopen('time.dat','wt');
        fprintf(fp,'%d\n',NE);
        if bd(k) == 1 && k~=cb
            fprintf(fp,'%d\n',k);
            for j = 1:numel(tot)
                fprintf(fp,'%.36e\n',tot(j));
            end
            fclose(fp);
            [xa,xb] = system('geteph.exe < time.dat > outcmd.txt');
            X = dlmread('fort.11');
            X = zentoijk(2,X);
        else
            X = zeros(6,NE);
            fclose(fp);
        end
        rkj(3*(k-1)+1:3*k,:) = X(1:3,:).*rconv;
        vkj(3*(k-1)+1:3*k,:) = X(4:6,:).*vconv;
    end
    
    tdim = tot;
    tadim = linspace(0,tf,NE);
    tfadim = tadim(end);
end
end

