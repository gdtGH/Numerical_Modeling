%% Solution of u_tt - c^2 Delta u = f with SEM-NI
%  Newmark time integration scheme with gamma=0.5 and beta=0.25

close all; clear; clc;

addpath Level_0
addpath Level_1
addpath Level_2


%% Initialize problem
dati;

%% Mesh generation
xy=[]; un=[]; D=[];

npdx  = nx + 1;
npdy  = ny + 1;
ldnov = npdx*npdy;
mn = ldnov;
ne = nex*ney;

% GLL weights and nodes + derivatives at GLL nodes
[x,wx] = xwlgl(npdx); [dx] = derlgl(x,npdx);
[y,wy] = xwlgl(npdy); [dy] = derlgl(y,npdy);

%nov = connectivity matrix
[nov] = cosnov_2d(npdx,nex,npdy,ney);
noe = nov(ldnov,ne);  %total degress of freedom

% Mesh generation
gammax = []; gammay = [];
[xx,yy,jacx,jacy,xy,ww,ifro]=mesh_2d(xa,xb,ya,yb,cb,nex,ney,npdx,npdy,...
                                     nov,x,wx,y,wy,gammax,gammay);

% Plot elements
for j = 1 : size(nov,2)
    plot(xy(nov([1 npdx npdx*npdy npdx*npdy-nx 1],j),1),...
        xy(nov([1 npdx npdx*npdy npdx*npdy-nx 1],j),2),'k'); 
    hold on;
end


% Plot grid nodes except dirichlet ones
for i = 1 : size(xy,1)
    if(ifro(i)==0) 
        plot(xy(i,1),xy(i,2),'ro'); hold on;
    elseif(ifro(i)==31)
        plot(xy(i,1),xy(i,2),'b*'); hold on;
    end
end


%% Matrix assembling
%stiffness
A = stiff_2d_se(npdx,nex,npdy,ney,nov,wx,dx,jacx,wy,dy,jacy);
A = v_c^2*A; %scaled by acoustic velocity
%mass
M = spdiags(ww,0,noe,noe);

% Volume force
fvol_0 = ff(xy(:,1),xy(:,2),0).* ww;

% Neumann boundary condition
fneu_0 = 0.*fvol_0;
[fneu_0] = neumannbc(xa,xb,ya,yb,fneu_0,h,jacx,jacy,wx,wy,xy,ifro,nov,0);

% Total forces : volume + boundary
f_0 = fvol_0 + fneu_0;


% Lists of internal, boundary, interface nodes. All lists are referred to
% global ordering on Omega:
% lbor: list of Dirichlet boundary nodes
% lint: list of internal nodes
% lintint: list of internal nodes which are internal to spectral elements
% lgamma: list of internal nodes which are on the interface between spectral elements

%% Setting Dirichlet boundary conditions on both matrix and r.h.s

[ldir,lint,lintint,lgamma,ifro]=liste(ifro,nov);

ub=zeros(noe,1);   %lifting of the Dirichlet datum
u0=zeros(noe,1);   %initial condition for u
v0=zeros(noe,1);   %initial condition for u_t

for i=1:noe
    u0(i)=u0ex(xy(i,1),xy(i,2));
    v0(i)=u1ex(xy(i,1),xy(i,2));
end

z0 = [u0;v0];

%% system matrix and rhs for newmark
% K = 1/dt^2 * M + 0.25*A;
% b = (1/dt^2 * M - 0.25*A)*u0 + 1/dt*M*v0  + 1/4*f_1 + 1/4*f_0;

gamma = 0.5;
gammat = 1-gamma;
beta = 0.25;
betat = (0.5-beta);

H1 = [M + beta*dt^2*A  M*0;  gamma*dt*A  M];
H0 = [M - dt^2*betat*A  dt*M; -gammat*dt*A  M]; 
K = H1;

%% Time Loop of newmark with gamma=1/2 and beta=1/4
if (param(8)~=0)
    fig=figure('Name','SEM-NI solution of u_{tt} - c^2 \Delta u = f','Visible','on');
    subplot(2,1,1)
    scatter3(xy(:,1),xy(:,2),u0,20,u0,'filled'); view(2); colorbar;
    xlim([xa xb]); ylim([ya yb]);  caxis([-1,1]);
    xlabel('x'); ylabel('y'); title(['solution at time t  ', num2str(t),' s.']);
    subplot(2,1,2)
    xplot = [xa:0.005:xb]'; yplot = 0*xplot;
    idx = knnsearch(xy,[xplot yplot]);
    plot(xplot,u0(idx),'Linewidth',2,'Color','b'); colorbar;
    ylim([-2 2]); xlim([xa xb]); caxis([-1,1]);
    title('solution at y=0'); xlabel('x'); ylabel('p(x,t)');
    pause(0.15);
end




fprintf('============================================================\n')
fprintf('Starting time-loop ... \n');
fprintf('============================================================\n')



for t = dt : dt : T-dt
    
    fprintf('time = %5.3e \n',t);

    % Volume force
    fvol_1 = ff(xy(:,1),xy(:,2),t).* ww;
    
    % Neumann boundary condition
    fneu_1 = 0.*fvol_1;
    [fneu_1] = neumannbc(xa,xb,ya,yb,fneu_1,h,jacx,jacy,wx,wy,xy,ifro,nov,t);
    
    % Total forces : volume + boundary
    f_1 = fvol_1 + fneu_1;

    % General time step
    rhs = H0*z0 + [dt^2*beta*f_1 + dt^2*betat*f_0; ...
                   dt*gamma*f_1 + dt*gammat*f_0];

    % apply dirichlet conditions
    for i=1:noe
        if(ifro(i)==1)
            ub(i) = g(xy(i,1),xy(i,2),t);
        end
    end
    
    ub_tot = [ub; 0*ub];
    
    %apply bc for the rhs
    if ~isempty(ldir)
        b = rhs - H1*ub_tot;
        for i = 1 : length(ldir)
            K(ldir(i),:) = 0;
            K(:,ldir(i)) = 0;
            K(ldir(i),ldir(i)) = 1;
            b(ldir(i)) = 0;
        end
    end
    
    z1 = K\b;
    %update solution when Diri is not zero
    if ~isempty(ldir)
        z1 = z1 + ub_tot;
    end
    u2 = z1(1:length(u0));    


    % plot the solution
    if (param(8)~=0)
        % old plot
        % [ha] = plot_sem_2d(fig,'surf',nex,ney,x,xx,jacx,y,yy,jacy,xy,ww,nov,u2,param(9));
        % fast plot
        gcp = subplot(2,1,1);
        scatter3(xy(:,1),xy(:,2),u2,20,u2,'filled'); view(2); colorbar; %axis equal;
        xlim([xa xb]); ylim([ya yb]);  caxis([-1,1]);
        xlabel('x'); ylabel('y');  title(['solution at time t  ', num2str(t),' s.'])
        gcp2 = subplot(2,1,2);
        plot(xplot,u2(idx),'Linewidth',2,'Color','b'); colorbar;
        ylim([-2 2]); xlim([xa xb]); caxis([-1,1]);
        title('solution at y=0'); xlabel('x'); ylabel('p(x,t)');
        pause(0.15);
        
    end
    
    % update for the next iteration
    z0 = z1;
    f_0 = f_1;
    
end

%% Error analysis
% errors are evaluated at the final time

if (param(4)==1)
    [err_inf,err_h1,err_l2] = errors_2d(x,wx,dx,xx,jacx,y,wy,dy,yy,jacy,...
        xy,ww,nov,ub,uex,uex_x,uex_y,param,T);
    fprintf('nx=%d,nex=%d,err_inf=%11.4e, err_h1=%11.4e,err_l2=%11.4e \n',...
        nx,nex,err_inf,err_h1,err_l2)
end





