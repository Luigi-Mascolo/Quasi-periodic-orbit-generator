function [X0,XF,X,yout,tout,iter,tf,STM] = orbcorr(X0,tspan,paramr3bp,optpar,fixd,fixc,flag)
%ORBCORR performs the CR3BP corrections to compute periodic orbits
%   INPUT:
%       X0      initial guess vector (x0, y0, z0, vx0, vy0, vz0) synodic RF
%       tspan   nondimensional integration period [0 tf]
%       paramr3bp      CR3BP parameters (e.g. mass ratio)
%       optpar  optimization parameters, i.e.
%           ril     relaxation parameter (0<ril<1)
%           opt113  integrator parameters for ode113
%           wb      waitbar flag (0 off, 1 on)
%           mode    correction mode (1 matrices, 2 direct equations)
%       fixd    vector identifying fixed values, i.e. fix = [2,4] means
%               that y and dot(x) cannot be changed during iterations in
%               the delta(X) vector. The STM is computed by taking all the
%               rows not included in fixd
%       fixc    vector identifying fixed values in the constraint vector.
%               I.e. upd = [1,3,5] means that x(T), z(T) and dot(y)(T)
%               should not be updated to convergence in the b vector.
%               THe STM is computed by taking all the columns not included
%               fixc
%
%   INTERNAL CALLS
%       crtbpeqmstab    CR3BP equations + STM
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory 
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022, 
% Politecnico di Torino. 
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator

ril = optpar.rilh;
opt113 = optpar.opt113;
wb = optpar.wb;
mode = optpar.modmat;
s = optpar.shooting;
fig = flag.fig;
fpr = flag.iter;
bck = flag.bck;
mu = paramr3bp.mu;

fix = fixd;
if nargin<7, fpr = false; fig = false; end
if nargin<5, fix = 0; end
if nargin<4, ril = 0.2; end

t0 = tspan(1);
tf = tspan(end);
y0 = zeros(42,1);
PHI0 = eye(6);
y0(37:42) = X0;
fl = true;
iter = 0;
tolnrho = 1e-8;

if wb == 1, h = waitbar(0,'Initializing'); end

if fig == 99 || fig == 98
    if fig == 98, delf = true;
    else, delf = false;
    end
    fig = true; nf = false;
else, nf = true;
end
ct = 'black';
cp = 'k.';
if fig
    if nf, figure; hold on; axis equal; end
    if strcmp(bck,'dark')
        darkBackground(h,[.2 .2 .2]); ct = 'white'; cp = 'w.';
    end
    title('CR3BP Orbit Computation','color',ct);
    [L2,~,~]    = fzero(@(x) x-(1-mu)/((x+mu)^2)-mu/((x-1+mu)^2),1.1*(1-mu));
    if flag.scl2 == 1, h2 = plot3(L2,0,0,cp); end
    view(30,15);
end

olderr = 1e6;
while fl
    iter = iter+1;
    for i=1:6
        for j=1:6
            y0(6*(i-1)+j)=PHI0(i,j);
        end
    end
    if strcmp(s,'ss')
        [tout,yout] = ode113(@(t,yy)crtbpeqmstab(t,yy,mu),tspan,y0,opt113);
    elseif strcmp(s,'ms')
        [tout,yout] = ode113(@(t,yy)crtbpeqmstab(t,yy,mu),tspan,y0,opt113);
    end
    
    temp = yout(end,1:36);
    STM = reshape(temp,[6,6])';
    
    x = yout(:,37);
    y = yout(:,38);
    z = yout(:,39);
    u = yout(:,40);
    v = yout(:,41);
    w = yout(:,42);
    
    if fig, hp = plot3(x,y,z,'w-'); hp.Color(4) = 0.5;
        if delf, pause(1e-1); delete(hp); end
    end
    
    xf = yout(end,37);
    yf = yout(end,38);
    zf = yout(end,39);
    uf = yout(end,40);
    vf = yout(end,41);
    wf = yout(end,42);
    XF = [xf yf zf uf vf wf];
    
    r1 = sqrt((x+mu).^2+y.^2+z.^2);     % m1-m3 distance
    r2 = sqrt((x-1+mu).^2+y.^2+z.^2);   % m2-m3 distance
    
    ax = 2.*v+x-(1-mu).*(x+mu)./(r1.^3)-mu.*(x-1+mu)./(r2.^3);
    ay = -2.*u+y-(1-mu).*y./(r1.^3)-mu.*(y./r2.^3);
    az = -(1-mu).*z./(r1.^3)-mu.*z./(r2.^3);
    
    axf = ax(end);
    ayf = ay(end);
    azf = az(end);
    
    if mode == 1 % matricial approach
        dX = [xf yf zf uf vf wf];
        dyf = yf;
        duf = uf;
        dwf = wf;
        PHIR = [STM(2,1) STM(2,3) STM(2,5) vf;
            STM(4,1) STM(4,3) STM(4,5) ax(end);
            STM(6,1) STM(6,3) STM(6,5) az(end)];
        DXF = [dyf; duf; dwf];
        
        XJ = [X0([1 3 5]); tf];
        if fix ~= 0
            if fix == 1
                PHIR = PHIR(:,2:end);
                XJ = [X0([3 5]); tf];
            elseif fix == 2
                PHIR = PHIR(:,[1 3 4]);
                XJ = [X0([1 5]); tf];
            elseif fix == 3
                PHIR = PHIR(:,[1 2 4]);
                XJ = [X0([1 3]); tf];
            elseif fix == 4
                PHIR = PHIR(:,1:3);
                XJ = X0([1 3 5]);
            elseif fix == 5
                PHIR = PHIR(:,2:3);
                XJ = X0([3 5]);
            end
        end
        XJ = XJ-PHIR.'*(inv(PHIR*PHIR.'))*DXF*ril;
        
        if fix == 0
            X0([1 3 5]) = XJ(1:3);
            tf = XJ(4);
        else
            if fix == 1, X0([3 5]) = XJ(1:2); tf = XJ(3);
            elseif fix == 2, X0([1 5]) = XJ(1:2); tf = XJ(3);
            elseif fix == 3, X0([1 3]) = XJ(1:2); tf = XJ(3);
            elseif fix == 4, X0([1 3 5]) = XJ(1:3);
            elseif fix == 5, X0([3 5]) = XJ(1:2);
            end
        end
        
    elseif mode == 2
        dX = -[xf yf zf uf vf wf]';
        Ac = [STM [uf; vf; wf; axf; ayf; azf]];
        take1 = 1:7;
        for k = 1:length(fixd)
            for j = 1:length(take1)
                if fixd(k)==take1(j)
                    take1(take1==fixd(k)) = [];
                end
            end
        end
        Ac = Ac(take1,:);
        take2 = 1:6;
        for k = 1:length(fixc)
            for j = 1:length(take2)
                if fixc(k)==take2(j)
                    take2(take2==fixc(k)) = [];
                end
            end
        end
        Ac = Ac(:,take2);
        dX = dX(take2);
        
        A = ([STM(4,3) STM(4,5); STM(6,3) STM(6,5)]-1/(vf)*[ax(end); az(end)]*[STM(2,3) STM(2,5)]);
        x = A\dX;
        dz = x(1);
        dvy = x(2);
        dtf = 1/(vf)*([STM(2,3) STM(2,5)]*x);
    end
    
    y0(37:42) = X0;
    tspan = linspace(t0,tf,numel(tspan));
    
    err = norm(DXF);
    if err<olderr, ril = ril+1.1e-3; end
    if err>olderr, ril = ril-1e-3; end
    if ril<=1e-4, ril = 1e-4; end
    olderr = err;
    if err>15, fprintf('Diverge.\n'); return; end
    if err<tolnrho, fl = false; end
    
    lt = log(tolnrho);
    le = log(err);
    
    if wb == 1, waitbar(le/lt,h,sprintf('%6.2f%%',100*le/lt)); end
    if fpr
        fprintf('%.6f %.6f %.6f %.6f %.6f %.6f %.4f %.4e %.6f\n',X0,tf,norm(err),ril)
        fprintf('%.6f %.6f %.6f %.6f %.6f %.6f %.4f %.4e %.6f\n',XF,tf,norm(err),ril)
    end
    
end
if wb == 1, delete(h); end
X = yout(:,37:42);
XF = [xf; yf; zf; uf; vf; wf];
end

