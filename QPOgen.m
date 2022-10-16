%% QPOgen - Quasi-Periodic Orbit generator
% QPOgen attempts at creating QPOs in a higher-fidelity model
% It uses ephemerides from JPL DE430 and has the following overall
% dependencies (in order of appearance)
%
%   darkbackground.m (graphical setting)
%   orbcorr.m
%       crtbpeqmstab.m
%   setmice.m   (SPICE/MICE kernels/files in the machine's program folder)
%   geteph.m
%       geteph.exe
%           EGM2008_to8_TideFree.txt
%           JPLEPH
%       zentoijk.m
%   rot2j2000.m
%   j20002rot.m
%   qpoeph.m
%       eqmstab.m
%       buildSTM.m
%       correctorb.m
%       jacobiValue3D.m
%   ijktozen.m
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator

clear; close all; fclose all; clc;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

fprintf('[1] Run\n[2] Help\n>>> ');
hr = input('');
if hr == 2
    fl = true;
    while fl
        clc
        fprintf('[ 1] orbcorr.m\n');
        fprintf('[ 2] crtbpeqmstab.m\n');
        fprintf('[ 3] setmice.m\n');
        fprintf('[ 4] geteph.m\n');
        fprintf('[ 5] zentoijk.m\n');
        fprintf('[ 6] rot2j2000.m\n');
        fprintf('[ 7] j20002rot.m\n');
        fprintf('[ 8] qpoeph.m\n');
        fprintf('[ 9] eqmstab.m\n');
        fprintf('[10] buildSTM.m\n');
        fprintf('[11] correctorb.m\n');
        fprintf('[12] jacobiValue3D\n');
        fprintf('[13] ijktozen.m\n');
        fprintf('[ 0] EXIT\n');
        ht = input('');
        switch ht
            case 1, help orbcorr
            case 2, help crtbpeqmstab
            case 3, help setmice
            case 4, help geteph
            case 5, help zentoijk
            case 6, help rot2j2000
            case 7, help j20002rot
            case 8, help qpoeph
            case 9, help eqmstab
            case 10, help buildSTM
            case 11, help correctorb
            case 12, help jacobiValue3D
            case 13, help ijktozen
            case 0, return;
            otherwise, warning('Program not found.');
        end
        fprintf('Press any key to continue.\n');
        inc = input('');
    end
end


md = pwd;
rd = strcat(md,'\resephx0');
fd = strcat(md,'\fig');
if ~exist('resephx0','dir'), mkdir('resephx0'); fprintf('[INFO] Missing result folder (resephx0). Creating...\n');
else, fprintf('[INFO] Existing result folder (resephx0).\n');
end
if ~exist('fig','dir'), mkdir('fig'); fprintf('[INFO] Missing image folder (fig). Creating...\n');
else, fprintf('[INFO] Existing image folder (fig).\n');
end
fprintf('[INFO] Pre-run checks/settings\n')
fprintf('[ACT] Save results (output results in resephx0)? [1] Y >>> ');
savres = input('');
if savres == 1, savres = true;
else, savres = false;
end

fprintf('[ACT] Overwrite results (if already existing)? [1] Y >>> ');
sovres = input('');
if sovres == 1, sovres = true;
else, sovres = false;
end

fprintf('[ACT] Save images (with overwrite, if existing)? [1] Y >>> ');
savim = input('');
if savim == 1, savim = true;
else, savim = false;
end
fprintf('DONE. Continue in ');
for k = 3:-1:1, fprintf('%d... ',k); pause(1); end
clc

NE = 1000;

fprintf('[INFO] ---------------- QPOs EPH BUILD ----------------\n');
fprintf('[INFO] --- Setting nondimensionalization parameters --- ');
mue = 398600.4415;      % Earth specific gravitational parameter [km3/s2]
mul = 4902.801076;      % Moon specific gravitational parameter [km3/s2]
mus = 1.32712440018e11; %
muj = 1.26686534e8;
rem = 384400;
res = 1.49597870691e8;
rconv = 6378.1363;
tconv = sqrt(rconv^3/(mue));
vconv = sqrt(mue/rconv);
rconvlag = rem;
tconvlag = sqrt(rconvlag^3/(mue+mul));
vconvlag = rconvlag/tconvlag;
rconvelio = 1.49597870691e8;
vconvelio = sqrt(1.32712440018e11/1.49597870691e8);
aconvelio = 1.32712440018d11/(1.49597870691e8^2);
tconvelio = vconvelio/aconvelio/86400;
fprintf('DONE\n');

conv.rl = rconvlag;
conv.tl = tconvlag;
conv.vl = vconvlag;
conv.re = rconvelio;
conv.te = tconvelio;
conv.ve = vconvelio;
conv.rt = rconv;
conv.tt = tconv;
conv.vt = vconv;

olderr = 1e6;
tolnrho = 1e-7;
fl = true;
fpr = true;
fix = 4;    %
t0 = 0;

fprintf('[ACT] Select bynary system:\n');
fprintf('\t[1] Earth-Moon L2 centered (EML2)\n');
fprintf('\t[2] Sun-Earth L2 centered (SEL2)\n>>> ');
scl2 = input('');
if scl2 == 1
    sys = 'EM';
    mut = mue+mul;
    mu  = mul/mut;
    a   = rem;
    [L2,~,exfl2]    = fzero(@(x) x-(1-mu)/((x+mu)^2)-mu/((x+(-1+mu))^2),1.1*(1-mu));
    Isp = 2000;
    tg = 15;
    Tlunasyn = 29.5305888531;
    Tnrhodim = Tlunasyn*2/9;
    AY = 4e4;
    strt1 = 'Earth-Moon barycentric RS';
    strt2 = 'Earth-Moon selenocentric RS';
else
    sys = 'SE';
    mut = mue+mus;
    mu  = mue/mut;
    a   = res;
    [L2,~,exfl2]    = fzero(@(x) x-(1-mu)/((x+mu)^2)-mu/((x+(-1+mu))^2),1.1*(1-mu));
    tg = 180;
    Isp = 3300;
    AY = 1e5;
    strt1 = 'Sun-Earth barycentric RS';
    strt2 = 'Sun-Earth geocentric RS';
end
clc
fprintf('----------- %s -----------\n',strt1)
fprintf('[ACT] Orbit type\n\t[1]Lyapunov\n\t[2]Halo\n>>> ');
typeqpo = input('');
if typeqpo~=1 && typeqpo ~=2, return; end
rconvlag = a;
tconvlag = sqrt(rconvlag^3/(mut));
vconvlag = rconvlag/tconvlag;
ndim = sqrt(rconvlag^3/(mut))/86400;    % conversione tempo adimensionale - giorni
Tlunasid = ndim*2*pi;
fprintf('[ACT] Insert guess orbital period (default: %.2f days) [days] >>> ',tg);
tf = input('');
if isempty(tf), tf = tg; end
tf = tf/ndim;
clc

paramr3bp.sys = sys;
paramr3bp.mu = mu;
paramr3bp.rc = rconvlag;
paramr3bp.tc = tconvlag;
paramr3bp.vc = vconvlag;
paramr3bp.LP = L2;
paramr3bp.nd = ndim;
paramr3bp.tg = tf;

reltol  = 1e-13;
abstol  = 1e-22;
opt113  = odeset('RelTol', reltol, 'AbsTol', abstol);
optpar.opt113 = opt113;
optpar.rilh = 0.5;
optpar.rilm = 0.2;
optpar.rill = 0.15;
optpar.wb = 1;
optpar.modmat = 1;
optpar.shooting = 'ss';
optpar.eph = 2;
ephm = optpar.eph;
optpar.nms = 10;

flag.iter = false;
flag.fig = savim;
flag.fmin = true;
flag.res = savres;
flag.tcor = false;
flag.iterplot = 0;
flag.scl2 = scl2;
flag.type = typeqpo;
flag.bck = 'dark';
flag.alt = false;
flag.locvid = true;
flag.cht = false;
flag.expl = 'n';

fprintf('[INFO] Analytic period orbit generation --- ');
K = (1-mu)/((L2+mu)^3)+mu/((L2-1+mu)^3);
om = 1/(sqrt(2))*sqrt(2-K+sqrt(9*(K^2)-8*K));
lam = 1/(sqrt(2))*sqrt(K-2+sqrt(9*(K^2)-8*K));
omz = sqrt(K);
c1 = (lam^2-2*K-1)/(2*lam);
c2 = (om^2+2*K+1)/(2*om);

paramr3bp.AY = AY;
paramr3bp.AX = AY./c2;
AX = AY./c2;

Ay = AY/rconvlag;
Ax = Ay/c2;
Az = 0/rconvlag;

phi = 0*pi/180;
psi = 0*pi/180;
t0 = 0;

t = linspace(0,2*pi,NE);

x0 = Ax*cos(om*t+phi);
y0 = Ay*sin(om*t+phi);
z0 = Az*cos(omz*t+psi);
vx0 = -Ax*om*sin(om*t+phi);
vy0 = -Ay*om*cos(om*t+phi);
vz0 = -Az*omz*sin(omz*t+psi);
fprintf('DONE\n');

fprintf('[INFO] Graphics settings --- ')
h = figure('renderer','painters','position',[100 100 1400 1000]); hold on; grid minor; view(30,30);
if scl2 == 1, LLL = 5e4/rconvlag;
elseif scl2 == 2, LLL = 5e4/rconvlag;
end
if strcmp(flag.bck,'dark')
    darkBackground(h,[.2 .2 .2]);
    if scl2 == 1
        subplot(1,2,1); hold on; grid minor; darkBackground(h,[.2 .2 .2]);
        subplot(1,2,2); hold on; grid minor; darkBackground(h,[.2 .2 .2]);
    end
end

if scl2 == 1, str2 = 'Moon';
else, str2 = 'Earth';
end
if scl2 == 1
    subplot(1,2,1);
    hm0 = plot3(1-mu,0,0,'w.');
    hmt = text(1-mu,0.02,0,str2,'fontsize',14,'color','w');
    hl2 = plot3(L2,0,0,'w.');
    hl2t = text(L2,0.02,0,'L_2','fontsize',14,'color','w');
    hp1 = plot3(x0+L2,y0,z0,'w-');
    xlim([1-mu-LLL/5 L2+LLL])
    ylim([-AY/rconvlag-LLL AY/rconvlag+LLL])
    if typeqpo == 1, zlim([-LLL LLL])
    elseif typeqpo == 2, zlim([-8e4/rconvlag 3*LLL])
    end
end

LE = (L2+LLL)-(1-mu-LLL/5);


if scl2 == 1, subplot(1,2,2); end

hm0 = plot3(0,0,0,'w.');
hmt = text(0,0.02,0,str2,'fontsize',14,'color','w');
hl2 = plot3(L2-(1-mu),0,0,'w.');
hl2t = text(L2-(1-mu),0.02,0,'L_2','fontsize',14,'color','w');
hp2 = plot3(x0+(L2-(1-mu)),y0,z0,'w-');
if scl2 == 1
    xlim([LLL+(L2-(1-mu))-LE LLL+(L2-(1-mu))])
else
    xlim([-LLL+(L2-(1-mu)) LLL+(L2-(1-mu))])
end
ylim([-AY/rconvlag-LLL AY/rconvlag+LLL])
if typeqpo == 1, zlim([-LLL LLL])
elseif typeqpo == 2, zlim([-8e4/rconvlag 3*LLL])
end
fprintf('DONE\n');

if scl2 == 1, subplot(1,2,1); end
fixval = 1;
ril = 0.5;

fprintf('[INFO] 3rd Order Richardson orbit correction (function orbcorr)\n')
X0 = [L2+x0(1) y0(1)/rconvlag z0(1) vx0(1) vy0(1) vz0(1)]';
X00 = X0;
tspan = linspace(0,tf/2,NE/2);
wb = 1;
flag.fig = 98;
optpar.rilh = 0.5;
fprintf('\t[INFO] Half orbit computation --- ');
[X0,~,XH,XTot,T,~,~] = orbcorr(X0,tspan,paramr3bp,optpar,fixval,[],flag);
fprintf('DONE\n');
tfg = 2*T(end);
tspan = linspace(0,tfg,NE);
fprintf('\t[INFO] Complete orbit computation --- ');
[X0,XF,X,XTot,T,~,tf] = orbcorr(X0,tspan,paramr3bp,optpar,fixval,[],flag);
fprintf('DONE\n');

if scl2 == 1, subplot(1,2,1);
    hp3 = plot3(X(:,1),X(:,2),X(:,3),'w-');
end
if scl2 == 1, subplot(1,2,2); end
hp4 = plot3(X(:,1)-(1-mu),X(:,2),X(:,3),'w-');


title('');

if scl2 == 1, cll = 0.2;
else, cll = 0.1;
end
if scl2 == 1
    hp1.Color(4) = cll;
    hp3.Color(4) = cll;
end
hp2.Color(4) = cll;
hp4.Color(4) = cll;

fprintf('[INFO] Set parameters for ephemerides and JPL SPICE\n');
optpar.rilh = 0.5;
mustar = mu;
fprintf('\t[INFO] Translation into primary RF --- ');
X0EPH = [X0(1)+mu; X0(2:6)];
fprintf('DONE\n');
fprintf('\t[INFO] Running SPICE/MICE kernel (function setmice) --- ');
setmice();
fprintf('DONE\n');
IDS = {'199','299','399','499','599','699','799','899','999','301','10'};
nameids = {'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto','Moon','Sun'};
str = 'J2000';

bd = zeros(1,numel(IDS));
if scl2 == 1
    bd([3 10 11]) = 1;
    cb  = 3;
    sb = 10;
else
    bd([3 10 11]) = 1;
    cb = 3;
    sb = 11;
end
fprintf('\t[INFO] List of included gravitational bodies:\n');
for k = 1:length(bd)
    if bd(k) == 1
        if k == cb, fprintf('\t\t(SPICE ID %s) \t%s (central body, primary 1)\n',IDS{k},nameids{k});
        elseif k==sb, fprintf('\t\t(SPICE ID %s) \t%s (binary body, primary 2)\n',IDS{k},nameids{k});
        else, fprintf('\t\t(SPICE ID %s) \t%s\n',IDS{k},nameids{k});
        end
    end
end
c = 0;
mub = [];
for k = 1:numel(bd)
    if bd(k) == 1
        switch k
            case 1, mub = [mub mumerc];
            case 2, mub = [mub muven];
            case 3, mub = [mub mue];
            case 4, mub = [mub mum];
            case 5, mub = [mub muj];
            case 6, mub = [mub musaturn];
            case 7, mub = [mub muu];
            case 8, mub = [mub munep];
            case 9, mub = [mub mupl];
            case 10, mub = [mub mul];
            case 11, mub = [mub mus];
        end
    end
end
for k = 1:length(bd)
    if bd(k) == 1
        c = c+1;
        mu(k) = mub(c);
    else
        mu(k) = 0;
    end
end
fprintf('\t[INFO] Gravitational parameter nondimentionalization wrt primary --- ');
if scl2 == 1
    mu = mu./(mu(3)); fprintf('Earth --- ');
else
    mu = mu./(mu(11)); fprintf('Sun --- ');
end
fprintf('DONE\n');

parameph.cb = cb;
parameph.sb = sb;
parameph.bd = bd;
parameph.sgp = mu;
mu0 = mu;
parameph.stage = 10;

parameph.IDS = IDS;
parameph.ref = str;
parameph.t00 = [162.033 162.153 162.274 162.394 162.514];
parameph.tolnrho = 1e-6;
parameph.vx0 = 0;

mub0 = mub;

if savres
    fpo = fopen('qpo.txt','at');
end

optpar.wb = 0;

fprintf('[INFO] Available starting dates (t0):\n')
date0 = {'2025 Oct 15 08:58:15.5433 UTC',...
    '2025 Oct 22 08:23:32.6902 UTC',...
    '2025 Oct 29 09:12:62.4800 UTC',...
    '2025 Nov 05 08:37:49.6269 UTC',...
    '2025 Nov 12 08:03:06.7738 UTC'};
for k = 1:numel(date0)
    fprintf('\t[%d] %s\n',k,date0{k});
end
fprintf('(Please do change manually starting dates for other t0)\n');

if scl2 == 1
    subplot(1,2,1); view(0,90); axis equal;
    subplot(1,2,2); view(0,90); axis equal;
else
    view(0,90); axis equal;
end

if scl2 == 1, subplot(1,2,1);
    xlim([1-mustar-LLL/5 L2+LLL])
    ylim([-AY/rconvlag-LLL/1.2 AY/rconvlag+LLL/1.2])
    if typeqpo == 1, zlim([-LLL LLL])
    elseif typeqpo == 2, zlim([-7.5e4/rconvlag LLL/2])
    end
end
title(strt1,'fontsize',12,'color','w');

LE = (L2+LLL)-(1-mustar-LLL/5);

if scl2 == 1, subplot(1,2,2); end
if scl2 == 1
    xlim([LLL+(L2-(1-mustar))-LE LLL+(L2-(1-mustar))])
else
    xlim([-8*LLL+(L2-(1-mustar)) 8*LLL+(L2-(1-mustar))])
end
ylim([-AY/rconvlag-LLL/1.2 AY/rconvlag+LLL/1.2])
if typeqpo == 1, zlim([-LLL LLL])
elseif typeqpo == 2, zlim([-7.5e4/rconvlag LLL/2])
end
title(strt2,'fontsize',12,'color','w');

if scl2 == 1, subplot(1,2,1); end

stage0 = 0;
no0 = 1;
X00S = X0;

if typeqpo == 1, fixval = 1;
else, fixval = 2; X0 = [X00S(1:2); -100/rconvlag; X00S(4:6)];
end

flag.fig = false;
tspan = linspace(0,tf/2,NE);
[X0,XF,XH,XTot,T,~,~] = orbcorr(X0,tspan,paramr3bp,optpar,fixval,[],flag);

X0 = XF;
tfg = 2*T(end);
tspan = linspace(0,tfg,NE);
[X0,XF,X,XTot,T,~,tf] = orbcorr(X0,tspan,paramr3bp,optpar,fixval,[],flag);
X0CR3BP = X0;
%
ht = plot3(X(:,1),X(:,2),X(:,3),'w-');
ht0 = plot3(X(1,1),X(1,2),X(1,3),'r.','markersize',12);

if size(X,1)~=6, X = X'; end
fprintf('[ACT] Select start date ID >>> ');
caso = input('');
t00 = parameph.t00(caso);

paramr3bp.kk = 1;
parameph.vx0 = 0;

X(1,:) = X(1,:)+mustar;
X(2,1) = 0;

strf = strcat('sx0_',num2str(t00,'%.3f'),'_',num2str(Isp,'%d'),'_',num2str(AY/1e5,'%.2f'),'.dat');
cd(rd);
AAA = dir;
flcnt = true;
for kf = 3:numel(AAA)
    if strcmp(strf,AAA(kf).name), flcnt = false; 
        fprintf('[INFO] File %s already found.\n',strf);
        if sovres, fprintf('[INFO] Overwriting\n');
        else, fprintf('[INFO] Overwriting not allowed. Stopping...\n');
        end
        break; 
    end
end
cd(md);

if flcnt || sovres
    
    fprintf('[INFO] Converting UTC to ET --- ');
    et0 = cspice_str2et(date0{caso});
    fprintf('DONE\n');
    fprintf('[INFO] Retrieve ephemerides (function geteph) --- ');
    [rkj,vkj,t,tadim] = geteph(ephm,et0,tf,str,ndim,IDS,bd,cb,NE,t00,conv);
    fprintf('DONE\n');
    
    r12 = [rkj((sb-1)*3+1:sb*3,:); vkj((sb-1)*3+1:sb*3,:)];
    X0SPRE = X(:,1);
    
    tfpre = tf;
    fprintf('[INFO] Conversion from rotating to J2000 frame (function rot2j2000) --- ');
    XI = rot2j2000(X,r12,mustar,rconvlag,vconvlag,sys);
    fprintf('DONE\n');
    X0IPRE = XI(:,1);
    
    if scl2 == 1
        XB = [rkj((sb-1)*3+1:sb*3,:)./rconvlag; vkj((sb-1)*3+1:sb*3,:)./vconvlag];
        XBS = j20002rot(XB,r12,mustar,rconvlag,vconvlag,sys);
        hm = plot3(XBS(1,1),XBS(2,1),XBS(3,1),'w.');
    end
    
    pause(1e-2);
    delete(ht); delete(ht0);
    
    rkj = rkj./rconvlag;
    vkj = vkj./vconvlag;
    X0 = XI(:,1);
    X0I = X0;
    
    y0 = zeros(42,1);
    PHI0 = eye(6);
    y0(37:42) = X0;
    
    tspan = linspace(0,tf,NE);
    
    fl = true;
    iter = 0;
    tf0 = tf;
    
    flag.iterplot = 0;
    
    fprintf('Continuing in ');
    for k = 5:-1:1, fprintf('%d... ',k); pause(1); end
    clc
    
    tic;
    
    fprintf('----- Quasi-Periodic Orbit generator ------\n');
    fprintf('--------------- MAIN CODE -----------------\n');
    fprintf('--- Mathematical Methods and Algorithms ---\n');
    fprintf('---- for Space Trajectory Optimization ----\n');
    fprintf('-- Mascolo, Luigi. Politecnico di Torino --\n');
    fprintf('-------------------------------------------\n');
    fprintf('[INFO] Selected start date %s\n',date0{caso});
    for stage = stage0:1:parameph.stage
        if stage==stage0, fprintf('\t[INFO] QPO in EPH binary frame (0%% perturbing bodies, PB)\n');
        else, fprintf('\t[INFO] %5.1f%% PB\n',stage/parameph.stage*100);
        end
        if stage/parameph.stage >= 0.1, flag.locvid = false; end
        flag.tcor = false;
        flag.fig = true;
        optpar.stage = stage;
        if abs(stage/parameph.stage*100)<1e-6 || abs(stage/parameph.stage*100)==50 || abs(stage/parameph.stage*100)== 100
            flag.iterplot = flag.iterplot+1;
            flag.locvid = true;
        end
        parameph.tolnrho = 1e-6;
        no = no0;
        
        tspan = linspace(0,tf0*no,NE);
        mu = parameph.sgp;
        for kj = 1:numel(bd)
            if bd(kj) == 1
                if kj == cb || kj == sb
                else
                    mu(kj) = mu(kj)./(parameph.stage).*stage;   % here it gradually increases third-body influences
                end
            end
        end
        parameph.sgp = mu;
        [rkj,vkj,t,tadim] = geteph(ephm,et0,tf0*no,str,ndim,IDS,bd,cb,NE,t00,conv);
        rkj = rkj./rconvlag;
        vkj = vkj./vconvlag;
        
        % available flags: x, y, z, vx, vy, vz, t
        % varyX0    includes variables allowed to vary (can be all)
        % checkXF   includes quantities to be zeroth (or made constant) at tf
        %           (cannot accept 't')
        if scl2 == 2
            varyX0 = {'x','vx','vy'};
            checkXF = {'y'};
        else
            if typeqpo == 1
                varyX0 = {'vy'};
                checkXF = {'y'};
            else
                varyX0 = {'x','vy'};
                checkXF = {'y'};
            end
        end
        desXF = zeros(size(checkXF));
        [XS,XI,XBS,XBI,STM,tf,rkj,vkj,X0I] = qpoeph(PHI0,X0I,rkj,vkj,parameph,paramr3bp,optpar,conv,tspan,t,et0,mustar,...
            varyX0,checkXF,desXF,no,flag,caso);
        X0S = XS(:,1);
        parameph.vx0 = X0S(4);
        XFS = XS(:,end);
        
        [r,th,ph,u,v,w] = ijktozen(XI);
        X0BS = XBS(:,1);
        XFBS = XBS(:,end);
        [rb,thb,phb,ub,vb,wb] = ijktozen(XBI);
        parameph.sgp = mu0;
        flag.tcor = true;
        
        if stage == parameph.stage
            fprintf('[INFO] Correcting tf and closing QPO... \n');
            fprintf('\t[INFO] Original guessed tf (GTF): %.4f days\n',tf0);
            flag.tcor = true;
            flk = true;
            stepk = 0.005;
            sg = 1;
            if scl2 == 2
                errold = XFS(1)-X0S(1);
                if errold>0, sg = -sg; end
            else
                errold = XFS(2);
                if errold>0, sg = -sg; end
            end
            
            no = no0;
            while flk
                tspan = linspace(0,tf0*no,NE);
                [rkj,vkj,t,tadim] = geteph(ephm,et0,tf0*no,str,ndim,IDS,bd,cb,NE,t00,conv);
                rkj = rkj./rconvlag;
                vkj = vkj./vconvlag;
                
                flag.alt = false;
                
                if scl2 == 2
                    varyX0 = {'x','vx','vy'};
                    checkXF = {'y'};
                else
                    if typeqpo == 1
                        varyX0 = {'vy'};
                        checkXF = {'x'};
                    else
                        varyX0 = {'vy'};
                        checkXF = {'x'};
                    end
                end
                
                desXF = zeros(size(checkXF));
                [XS,XI,XBS,XBI,STM,tf,rkj,vkj,X0I,XBS2,r122,rkj2,vkj2] = qpoeph(PHI0,X0I,rkj,vkj,parameph,paramr3bp,optpar,conv,tspan,t,et0,mustar,...
                    varyX0,checkXF,desXF,no,flag,caso);
                %                 X0I = XI(:,1);
                X0S = XS(:,1);
                XFS = XS(:,end);
                
                if scl2 == 2, errxf = XFS(1)-X0S(1); ERTOL = 5e-6;
                else, errxf = XFS(2); ERTOL = 1e-4;
                end
                if ((XFS(1)<X0S(1)&& sg<0) || (XFS(1)>X0S(1)&& sg>0)) && scl2 == 2
                    sg = -sg; stepk = stepk/1.5;
                end
                if scl2 == 1 && ((XFS(2)>0 && sg>0 || XFS(2)<0 && sg<0))
                    sg = -sg; stepk = stepk/1.5;
                end
                errold = errxf;
                if sg>0, strbf = 'forward ';
                else, strbf = 'backward';
                end
                fprintf('\t[INFO] %6.2f%% GTF. Search %s (step %.2e) - error % .2e\n',no*100,strbf,stepk,errxf);
                xs = XS(1,:);
                ys = XS(2,:);
                zs = XS(3,:);
                
                no = no+sg*stepk;
                if no~=1, flag.locvid = false; flag.cht = true; end
                if abs(errxf)<ERTOL, flk = false; fprintf('[INFO] CONVERGED. Closing...\n'); end
                %                     if no<2 && abs(errxf)<1e-2, no = no*2; end;
            end
            
            flag.locvid = true;
            if scl2 == 2
                varyX0 = {'x','vx','vy'};
                checkXF = {'y'};
            else
                if typeqpo == 1
                    varyX0 = {'vy'};
                    checkXF = {'y'};
                else
                    varyX0 = {'x','vy'};
                    checkXF = {'y'};
                end
            end
            desXF = zeros(size(checkXF));
            [XS,XI,XBS,XBI,STM,tf,rkj,vkj,X0I] = qpoeph(PHI0,X0I,rkj,vkj,parameph,paramr3bp,optpar,conv,tspan,t,et0,mustar,...
                varyX0,checkXF,desXF,no,flag,caso);
            X0S = XS(:,1);
            parameph.vx0 = X0S(4);
            XFS = XS(:,end);
        end
        
        % check if exist more bodies other than cb/sb
        ccc = 0;
        for ppl = 1:numel(mu0)
            if mu0(ppl)~=0, ccc = ccc+1; end
        end
        if ccc == 2, break; end
    end
    
    [r0,th0,ph0,u0,v0,w0] = ijktozen(XI(:,1));
    u0 = u0*vconvlag/vconv;
    v0 = v0*vconvlag/vconv;
    w0 = w0*vconvlag/vconv;
    
    if savres
        fprintf('----- ADDITIONAL OUTPUT: SAVING RESULTS -----\n');
        fprintf('[INFO] Update jacobi.txt file --- ');
        j0 = jacobiValue3D(XS(:,1),mustar);
        if exist('jacobi.txt','file'), fp = fopen('jacobi.txt','at');
        else, fp = fopen('jacobi.txt','wt');
        end
        fprintf(fp,'%.3f %1d %.6f\n',t00,AY/1e5,j0);
        fclose(fp);
        fprintf('DONE\n');
        
        fprintf('[INFO] Save sx0[...].dat initial state [r,th,ph,u,v,w] ---  ');
        strx0 = strcat('sx0_',num2str(t00,'%.3f'),'_',num2str(Isp),'_',num2str(AY/(1e5),'%.2f'),'.dat');
        fp = fopen(strx0,'wt');
        fprintf(fp,'% .24f\n% .24f\n% .24f\n% .24f\n% .24f\n% .24f\n',r0*rconvlag/rconv,th0*pi/180,ph0*pi/180,u0,v0,w0);
        fclose(fp);
        str0 = strcat(md,'\',strx0);
        strf = strcat(rd,'\',strx0);
        flt = true;
        while flt
            try movefile(str0,strf); flt = false;
            catch, flt = true;
            end
        end
        fprintf('DONE\n');
        
        fprintf('[INFO] Save traj_I_sx0[...].dat trajectory, inertial RF --- ');
        strx0 = strcat('traj_I_sx0_',num2str(t00,'%.3f'),'_',num2str(Isp),'_',num2str(AY/(1e5),'%.2f'),'.dat');
        fp = fopen(strx0,'wt');
        for klo = 1:size(XI,2)
            fprintf(fp,'%.16e %.16e %.16e %.16e %.16e %.16e %.16e\n',tspan(klo)*ndim,XI(1:3,klo).*rconvlag,XI(4:6,klo).*vconvlag);
        end
        fclose(fp);
        str0 = strcat(md,'\',strx0);
        strf = strcat(rd,'\',strx0);
        flt = true;
        while flt
            try movefile(str0,strf); flt = false;
            catch, flt = true;
            end
        end
        fprintf('DONE\n');
        
        fprintf('[INFO] Save eph_I_sx0[...].dat QPO eph., inertial RF ---    ');
        strx0 = strcat('eph_I_sx0_',num2str(t00,'%.3f'),'_',num2str(Isp),'_',num2str(AY/(1e5),'%.2f'),'.dat');
        fp = fopen(strx0,'wt');
        for klo = 1:size(rkj,2)
            fprintf(fp,'%.16e ',tspan(klo)*ndim);
            for kll = 1:size(rkj,1)
                fprintf(fp,'% .16e ',rkj(kll,klo)*rconvlag);
            end
            for kll = 1:size(rkj,1)
                fprintf(fp,'% .16e ',vkj(kll,klo)*vconvlag);
            end
            fprintf(fp,'\n');
        end
        fclose(fp);
        str0 = strcat(md,'\',strx0);
        strf = strcat(rd,'\',strx0);
        flt = true;
        while flt
            try movefile(str0,strf); flt = false;
            catch, flt = true;
            end
        end
        fprintf('DONE\n');
        
        fprintf('[INFO] Save traj_S_sx0[...].dat trajectory, synodic RF ---  ');
        strx0 = strcat('traj_S_sx0_',num2str(t00,'%.3f'),'_',num2str(Isp),'_',num2str(AY/(1e5),'%.2f'),'.dat');
        fp = fopen(strx0,'wt');
        for klo = 1:size(XS,2)
            fprintf(fp,'%.16e %.16e %.16e %.16e %.16e %.16e %.16e\n',tspan(klo)*ndim,XS(1:3,klo).*rconvlag,XS(4:6,klo).*vconvlag);
        end
        fclose(fp);
        str0 = strcat(md,'\',strx0);
        strf = strcat(rd,'\',strx0);
        flt = true;
        while flt
            try movefile(str0,strf); flt = false;
            catch, flt = true;
            end
        end
        fprintf('DONE\n');
        
        fprintf('[INFO] Updating qpo.txt file --- ');
        fprintf(fpo,'%.4e 000.000 % .8e % .8e % .8e % .8e % .8e % .8e % .8e %.4f \n',AY,X0SPRE,tfpre,tfpre*ndim);
        fprintf(fpo,'%.4e %7.3f % .8e % .8e % .8e % .8e % .8e % .8e % .8e %.4f \n',AY,t00,X0S,tf,tf*ndim);
        fprintf('DONE\n');
        
    end
    
    if savim
        fprintf('----- ADDITIONAL OUTPUT: SAVING FIGURES -----\n');
        fprintf('[INFO] Saving graphical output --- ');
        switch t00
            case 162.033, stri = 'c1';
            case 162.153, stri = 'c2';
            case 162.274, stri = 'c3';
            case 162.394, stri = 'c4';
            case 162.514, stri = 'c5';
        end
        if scl2 == 1
            str0 = strcat('eml2-',stri,'-',num2str(AY/(1e5),'%.2f'),'.jpeg');
        else
            str0 = strcat('sel2-',stri,'-',num2str(AY/(1e5),'%.2f'),'.jpeg');
        end
        h = figure(1);
        print(str0,'-djpeg','-r300');
        fl = true;
        while fl
            try
                movefile(strcat(md,'\',str0),...
                    strcat(fd,'\',str0)); 
                fl = false;
            catch, fl = true; 
            end
        end
        fprintf('DONE\n');
    end
    
end

if savres, fclose(fpo); end

fprintf('[INFO] Performing after run cleaning...\n');
c= 0;
if exist('outcmd.txt','file'); c = c+1; delete('outcmd.txt'); end
if exist('time.dat','file'); c = c+1; delete('time.dat'); end
if exist('in.txt','file'); c = c+1; delete('in.txt'); end
if exist('U30','file'); c = c+1; delete('U30'); end
if c>0, fprintf('%d temp file(s) eliminated.\n',c); end
fprintf('DONE.\n');

et = toc;
[h,m] = deal(0);
while et>=60
    et = et-60;
    m = m+1;
    if m>=60, m = m-60; h=h+1; end
end
fprintf('Elapsed time: %1d h %02d min %05.2f sec\n',h,m,et);

