function [XS,X,XBS,XB,STM,tf,rkj,vkj,X0I,XBS2,r122,rkj2,vkj2] = qpoeph(PHI0,X0I,rkj,vkj,parameph,paramr3bp,optpar,conv,tspan,~,et0,...
    mustar,varyX0,checkXF,desXF,~,flag,caso)
%QPOEPH tries to compute a quasi-periodic orbit in an higher fidelity model
%based on planetary ephemeris.
%
%   INPUT:
%       PHI0    Initial state-transition matrix (eye matrix)
%       X0I     Vector containing SC inertial initial state
%       rkj     Matrix containing celestial bodies positions (wrt central)
%       vkj     Matrix containing celestial bodies velocities (wrt central)
%       parameph    struct with celestial bodies/optimization parameters
%       paramr3bp   CR3BP parameters (e.g. mass ratio)
%       optpar  integration/case specific parameters
%       conv    struct containing conversion parameters
%       tspan   nondimensional integration time
%       et0     starting epoch
%       mustar  CR3BP specific mass parameter (smaller primary over total)
%       varyx0  string vector containing information about initial
%               variables allowed to change during differential correction
%       checkXF string vector containing information about terminal
%               constraints to be fulfilled during differential correction
%       desXF   vector containing the numerical values of constraints
%               included in checkXF (usually zeros)
%       flag    struct containing cases choices
%       caso    number identifying the initial epoch among the listed ones
%
%   INTERNAL CALLS
%       ode2    fixed-step explicit Euler second order integrator
%       eqmstab differential corrector integrator (with STM)
%       buildSTM    function to rebuild the complete STM
%       correctorb  differential correction procedure
%       jacobiValue3D   function to compute the 3-dimensional Jacobi value
%
% Reference:
%
% Mascolo, Luigi. Mathematical Methods and Algorithms for Space Trajectory
% Optimization, unpublished doctoral dissertation as of 15 Oct 2022,
% Politecnico di Torino.
%
% https://github.com/Luigi-Mascolo/Quasi-periodic-orbit-generator


ffmin = flag.fmin;
fpr = flag.iter;
scl2 = flag.scl2;
type = flag.type;

bd = parameph.bd;
sb = parameph.sb;
cb = parameph.cb;
mu = parameph.sgp;
IDS = parameph.IDS;
str = parameph.ref;
t00 = parameph.t00(caso);
tolnrho = parameph.tolnrho;

rconvlag = paramr3bp.rc;
vconvlag = paramr3bp.vc;
sys = paramr3bp.sys;
ndim = paramr3bp.nd;
AY = paramr3bp.AY;
AYY = AY(paramr3bp.kk);

ril = optpar.rill;
mode = optpar.eph;
stage = optpar.stage;


lt = {'k-*','k-v','k-s','k-o','k-d','k-p','k-|','k-x','k->','k-<','k-^'};
if strcmp(flag.bck,'dark')
    lt = {'w-*','w-v','w-s','w-o','w-d','w-p','w-|','w-x','w->','w-<','w-^'};
end
if ~exist('flag.vary','var'), lt = {'w-','w-','w-','w-','w-','w-','w-','w-','w-','w-'}; end

iter = 0;
fl = true;
if nargin<17, desXF = []; end
t0 = tspan(1);
tf = tspan(end);
NE = numel(tspan);
ril0 = ril;
desXF0 = desXF;
cnup = 0;

fltempo = true;
cnlt = 0;
fleph = true;
XBS2 = [];
r122 = [];
if ~exist('flag.subvid','var'), flv = true;
else, if flag.subvid, flv = true; else, flv = false; end
end

while fl
    iter = iter+1;
    
    y0(1:36) = reshape(PHI0,[36 1]);
    y0(37:42) = X0I;
    if size(y0,1) == 1, y0 = y0'; end
    
    yout = ode2(@(t,yy)eqmstab(t,yy,bd,cb,rkj,mu,sb,tspan),tspan,y0);
    x = yout(:,37);
    y = yout(:,38);
    z = yout(:,39);
    vx = yout(:,40);
    vy = yout(:,41);
    vz = yout(:,42);
    X = [x y z vx vy vz]';
    XB = [rkj((sb-1)*3+1:sb*3,:); vkj((sb-1)*3+1:sb*3,:)];
    X0I = X(:,1);
    XFI = X(:,end);
    X0BI = XB(:,1);
    XFBI = XB(:,end);
    
    STM = buildSTM(yout,rkj,bd,cb,mu);
    eigstm = eig(STM(:,1:6));
    
    % RICONVERTIRE NEL SINODICO LE CONDIZIONI FINALI PER AVERE COME
    % INIZIALI
    %     r12 = cspice_spkezr(IDS{sb},t,str,'none',IDS{cb});
    r12 = [rkj((sb-1)*3+1:sb*3,:).*rconvlag; vkj((sb-1)*3+1:sb*3,:).*vconvlag];
    
    XS = j20002rot(X,r12,mustar,rconvlag,vconvlag,sys);
    X0S = XS(:,1);
    XFS = XS(:,end);
    XBS = j20002rot(XB,r12,mustar,rconvlag,vconvlag,sys);
    X0BS = XBS(:,1);
    XFBS = XBS(:,end);
    
    % USER EXPERIENCE
    % Hand-corrected quantities when transitioned to EPH
    % SEL2:
    %   - u0 slightly negative at first iteration
    %   - u0 increased (<1%) if first guess xf<x0 (it should be xf>x0)
    %   - u0 decreased (<0.1%) if first guess xf>>x0 (it should be xf>x0)
    % EML2:
    %   - u0 corrected by x0_barycentric/x0_synodic quantity (small)
    
    if iter == 1 && parameph.vx0 == 0
        if scl2 == 2, X0S(4) = -X0BS(4)*(paramr3bp.LP-(1-paramr3bp.mu));
        else
            X0S(4) = X0BS(4)*(X0BS(1)/X0S(1));
        end
    elseif iter == 1 && parameph.vx0 ~= 0, X0S(4) = parameph.vx0;
    end
    dd = abs(XFS(1)-X0S(1));
    if scl2 == 2
        if XFS(1)<X0S(1) && iter == 1
            X0S(4) = X0S(4)+dd/100;
        elseif XFS(1)>X0S(1)*1.001 && iter == 1, X0S(4) = X0S(4)-dd/1000;
        end
    end
    
    xs = XS(1,:);
    ys = XS(2,:);
    zs = XS(3,:);
    
    % EXPERIMENTAL
    % Double alternate correction (odd and even iterations fix different
    % final values)
    if flag.alt
        if mod(iter,2) == 0
            varyX0 = {'vy'};
            checkXF = {'vx'};
            desXF = zeros(size(checkXF));
        else
            varyX0 = {'vy'};
            checkXF = {'x'};
            desXF = zeros(size(checkXF));
        end
    end
    
    if any(contains(varyX0,'t'))
        X0S = [X0S; tf];
    end
    
    [rkj2,vkj2,~] = geteph(mode,et0,tf,str,ndim,IDS,bd,sb,NE,t00,conv);
    rkj2 = rkj2./rconvlag;
    vkj2 = vkj2./vconvlag;
    r122 = [rkj2((cb-1)*3+1:cb*3,:).*rconvlag; vkj2((cb-1)*3+1:cb*3,:).*vconvlag];
    XB2 = [rkj((cb-1)*3+1:cb*3,:); vkj((cb-1)*3+1:cb*3,:)];
    XBS2 = j20002rot(XB2,r122,mustar,rconvlag,vconvlag,sys);
    X0BS2 = XBS(:,1);
    XFBS2 = XBS(:,end);
    
    for k = 1:numel(checkXF)
        if strcmp(checkXF{k},'x')
            if scl2 == 1, desXF(k) = X0S(1)+(XFBS(1)-X0BS(1));
            else, desXF(k) = X0S(1)+(XFBS2(1)-X0BS2(1));
            end
        elseif strcmp(checkXF{k},'y'), desXF(k) = X0S(2);
        elseif strcmp(checkXF{k},'z'), desXF(k) = X0S(3);
        elseif strcmp(checkXF{k},'vx'), desXF(k) = X0S(4);
        elseif strcmp(checkXF{k},'vy'), desXF(k) = X0S(5);
        elseif strcmp(checkXF{k},'vz'), desXF(k) = X0S(6);
        end
    end
    
    [X0NEW,~,~,XXF] = correctorb(X0S,XFS,X0BS,XFBS,STM,varyX0,checkXF,desXF,ril);
    X0S = X0NEW;
    
    DXF = XXF;
    
    if any(contains(varyX0,'t'))
        tf = X0S(end);
        [rkj,vkj,t] = geteph(mode,et0,tf,str,ndim,IDS,bd,cb,NE,t00,conv);
        rkj = rkj./rconvlag;
        vkj = vkj./vconvlag;
        X0S = X0S(1:6);
    end
    
    err = norm(DXF);
    if iter == 1, errold = err;
    else
        if errold>err, ril = ril+ril0/1000;cnup = cnup+1;
            if cnup>10, cnup = 0; ril = ril+ril0/3; end
        else, ril = ril-ril0/1000; cnup = 0;
        end
        if ril<ril0/1000||ril<1e-3
            if ril0/1000>1e-3, ril = ril0/1000;
            else, ril = 1e-3;
            end
        end
        errold = err;
    end
    if fpr
        fprintf('% .6f % .6f % .6f % .6f % .6f % .6f %.4f %.4e %.2e\n',X0S,tf,norm(err),ril)
        fprintf('% .6f % .6f % .6f % .6f % .6f % .6f %.4f %.4e %.2e\n',XFS,tf,norm(err),ril)
    end
    if err>100, fprintf('Diverge.\n'); return; end
    
    %     r12 = cspice_spkezr(IDS{sb},t(1),str,'none',IDS{cb});
    r12 = [rkj((sb-1)*3+1:sb*3,:).*rconvlag; vkj((sb-1)*3+1:sb*3,:).*vconvlag];
    
    
    
    X0I = rot2j2000(X0S,r12,mustar,rconvlag,vconvlag,sys);
    if err<tolnrho && ffmin, fl = false; end
    if iter>5000 && err<1-5, tolnrho = tolnrho*10; end
    tspan = linspace(0,tf,NE);
    if scl2 == 1, cll = 0.5;
    else, cll = 0.25;
    end
    
    if flv
        if scl2 == 1, subplot(1,2,1); end
        %     if scl2 == 1 && (mod(iter,200) == 0 || iter == 1)
        if (flag.locvid) || (~fl && ~flag.locvid)
            ix = flag.iterplot;
            mk = lt{ix};
            mk = mk(end-1:end);
            if scl2 == 1
                subplot(1,2,1);
                hp1 = plot3(xs,ys,zs,lt{ix},'markerindices',1:round(numel(xs)/10):numel(xs),'markersize',6);
                hpp1 = plot3(xs(1),ys(1),zs(1),mk,'color','g','markersize',7);
                hpp1e = plot3(xs(end),ys(end),zs(end),mk,'color','r','markersize',7);
                subplot(1,2,2)
                hp2 = plot3(xs-XBS(1,:),ys,zs,lt{ix},'markerindices',1:round(numel(xs)/10):numel(xs),'markersize',6);
                hpp2 = plot3(xs(1)-XBS(1,1),ys(1),zs(1),mk,'color','g','markersize',7);
                hpp2e = plot3(xs(end)-XBS(1,end),ys(end),zs(end),mk,'color','r','markersize',7);
            else
                hp2 = plot3(xs-XBS2(1,:),ys,zs,lt{ix},'markerindices',1:round(numel(xs)/10):numel(xs),'markersize',6);
                hpp2 = plot3(xs(1)-XBS2(1,1),ys(1),zs(1),mk,'color','g','markersize',7);
                hpp2e = plot3(xs(end)-XBS2(1,end),ys(end),zs(end),mk,'color','r','markersize',7);
            end
            pause(1e-2)
            
            if scl2 == 1, delete(hp1); delete(hpp1); delete(hpp1e); end
            delete(hp2); delete(hpp2); delete(hpp2e);
        end
        
        if ~fl && flag.fig && ~flag.tcor
            if abs(stage/parameph.stage*100)<1e-6 || abs(stage/parameph.stage*100)==50 || abs(stage/parameph.stage*100)== 100
                ix = flag.iterplot;
                if scl2 == 1
                    subplot(1,2,1)
                    hp = plot3(xs,ys,zs,lt{ix},'markerindices',1:round(numel(xs)/10):numel(xs),'markersize',6);
                    hp.Color(4) = cll;
                    subplot(1,2,2)
                    hp = plot3(xs-XBS(1,:),ys,zs,lt{ix},'markerindices',1:round(numel(xs)/10):numel(xs),'markersize',6);
                else
                    hp = plot3(xs-XBS2(1,:),ys,zs,lt{ix},'markerindices',1:round(numel(xs)/10):numel(xs),'markersize',6);
                end
                hp.Color(4) = cll;
            end
            pause(1e-2)
        end
    end
    if type == 1
        strtype = 'CR3BP Lyap: ';
    else
        strtype = 'CR3BP Halo: ';
    end
    switch t00
        case 162.033, strt0 = '15/10/2025 08:58:15.5433 UTC';
        case 162.153, strt0 = '22/10/2025 08:23:32.6902 UTC';
        case 162.274, strt0 = '29/10/2025 09:12:62.4800 UTC';
        case 162.394, strt0 = '05/11/2025 08:37:49.6269 UTC';
        case 162.514, strt0 = '12/11/2025 08:03:06.7738 UTC';
    end
    
    X0N = [X0S; tf];
    for k = 1:numel(X0N)
        X0SD{k} = string(num2str(X0N(k),'%.5e'));
        if abs(X0N(k))<1e-8, X0SD{k} = '0.00000e+00'; end
    end
    %     X0SD = string(num2str(X0N,'%.6e'));
    XFN = XFS;
    for k = 1:numel(XFN)
        if abs(XFN)<1e-8, XFN(k) = 0; end
    end
    XFSD = string(num2str(XFN,'%.4e'));
    dx0v = zeros(1,7);
    xfv = zeros(1,6);
    variables = {'x','y','z','vx','vy','vz','t'};
    for k = 1:length(variables)
        dx0v(k) = 1;
        for j = 1:length(varyX0)
            if strcmp(variables{k},varyX0{j})
                dx0v(k) = 0;
            end
        end
    end
    for k = 1:length(variables)-1
        xfv(k) = 0;
        for j = 1:length(checkXF)
            if strcmp(variables{k},checkXF{j})
                xfv(k) = 1;
            end
        end
    end
    
    if flag.cht, dx0v(end) = 0; end
    variables = {'x','y','z','v_x','v_y','v_z','{\Delta}t'};
    strx01 = '';
    strx02 = '';
    for k = 1:numel(dx0v)
        if dx0v(k) == 0
            if isempty(strx01), strx01 = strcat('(',variables{k},{') '},X0SD{k});
            else, strx01 = strcat(strx01,{', '},'(',variables{k},{') '},X0SD{k});
            end
        else
            if isempty(strx02), strx02 = strcat('(',variables{k},{') '},X0SD{k});
            else, strx02 = strcat(strx02,{', '},'(',variables{k},{') '},X0SD{k});
            end
        end
    end
    strx01 = string(strx01);
    strx02 = string(strx02);
    strxf1 = '';
    strxf2 = '';
    for k = 1:numel(xfv)
        if xfv(k) == 1
            if isempty(strxf1), strxf1 = strcat('(',variables{k},{') '},XFSD{k});
            else, strxf1 = strcat(strxf1,{', '},'(',variables{k},{') '},XFSD{k});
            end
        else
            if isempty(strxf2), strxf2 = strcat('(',variables{k},{') '},XFSD{k});
            else, strxf2 = strcat(strxf2,{', '},'(',variables{k},{') '},XFSD{k});
            end
        end
    end
    strxf1 = string(strxf1);
    strxf2 = string(strxf2);
    
    j0 = jacobiValue3D(XS(:,1),mustar);
    
    if flag.cht
        sgtitle({
            sprintf('%s A_Y = %.2e km, t_0 = %s || %3d%% pert. bodies, err %.4e',strtype,AYY,strt0,round(stage/parameph.stage*100),err)
            sprintf('{\\color{green}{\\Delta}T = %.4f (%.4f days)}, Jacobi_C = %.6f',tf,tf*ndim,j0)
            sprintf('X_{0,var} = {\\color{green}[%s]}, X_{0,fix} = [%s]',strx01,strx02)
            sprintf('X_{f,chk} = {\\color{red}[%s]}, X_{0,free} = [%s]',strxf1,strxf2)
            sprintf('')
            },'color','white');
    else
        if strcmp(flag.expl,'n')
            sgtitle({
                sprintf('%s A_Y = %.2e km, t_0 = %s || %3d%% pert. bodies, err %.4e',strtype,AYY,strt0,round(stage/parameph.stage*100),err)
                sprintf('{\\Delta}T = %.4f (%.4f days), Jacobi_C = %.6f',tf,tf*ndim,j0)
                sprintf('X_{0,var} = {\\color{green}[%s]}, X_{0,fix} = [%s]',strx01,strx02)
                sprintf('X_{f,chk} = {\\color{red}[%s]}, X_{0,free} = [%s]',strxf1,strxf2)
                sprintf('')
                },'color','white');
        elseif strcmp(flag.expl,'x')
            sgtitle({
                sprintf('%s A_Y = %.2e km, t_0 = %s || %3d%% pert. bodies, err %.4e',strtype,AYY,strt0,round(stage/parameph.stage*100),err)
                sprintf('{\\Delta}T = %.4f (%.4f days), Jacobi_C = %.6f, \\color{cyan}\\rightarrow x search',tf,tf*ndim,j0)
                sprintf('X_{0,var} = {\\color{green}[%s]}, X_{0,fix} = [%s]',strx01,strx02)
                sprintf('X_{f,chk} = {\\color{red}[%s]}, X_{0,free} = [%s]',strxf1,strxf2)
                sprintf('')
                },'color','white');
        elseif strcmp(flag.expl,'z')
            sgtitle({
                sprintf('%s A_Y = %.2e km, t_0 = %s || %3d%% pert. bodies, err %.4e',strtype,AYY,strt0,round(stage/parameph.stage*100),err)
                sprintf('{\\Delta}T = %.4f (%.4f days), Jacobi_C = %.6f, \\color{cyan}\\rightarrow z search',tf,tf*ndim,j0)
                sprintf('X_{0,var} = {\\color{green}[%s]}, X_{0,fix} = [%s]',strx01,strx02)
                sprintf('X_{f,chk} = {\\color{red}[%s]}, X_{0,free} = [%s]',strxf1,strxf2)
                sprintf('')
                },'color','white');
        end
    end
    pause(1e-2)
    %     end
end
end

