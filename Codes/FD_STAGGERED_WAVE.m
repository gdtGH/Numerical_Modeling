function FD_STAGGERED_WAVE(~)

clearvars; close all; clc;

% Solution of the problem:
%  rho v_t = sigma_x  in (0,L) x (0,T)
%  sigma_t = mu v_x   in (0,L) x (0,T)
%  v(0,t)  = phi(t)   in (0,T)
%  v(L,t)  = psi(t)   in (0,T)
%  sigma(0,t) = s_0   in (0,L)
%  v(0,t)     = v_0   in (0,L)
%
% with the staggered finite difference scheme
% v = u_t and sigma = mu u_x
% u is the solution of rho u_tt - (mu u_x)_x = 0
%
% To be done:
% - add external force f to the system of eqs
% - add boundary conditions for sigma
%
%
%  Model problem: Plane wave solution in (0,1)
%  u(x,t) = sin(2pix)sin(2pit)
%  v(x,t) = 2pi sin(2pix)cos(2pit)
%  sigma(x,t) = mu 2pi cos(2pix)sin(2pit)


% Parameters
L = 1;              % Length of the domain
T = 2;              % Total simulation time
mu = 1;             %
rho = 1;
dx = 0.01;          % Spatial step size
dt = 0.01;         % Temporal step size

% Discretization
xv = 0 : dx/2 : L;    % Spatial grid
tv = 0 : dt/2 : T; % Temporal grid for v

nx = length(xv);
nt = length(tv);

%grid points velocity
[Xv,Yv] = meshgrid(xv(1:2:end),tv(2:2:end));
% grid points stress
[Xs,Ys] = meshgrid(xv(2:2:end),tv(1:2:end));

if (nx < 50 && nt <50)
    figure(1)
    plot(Xv,Yv,'b*')
    hold on; plot(Xs,Ys,'ro')
    xlabel('x'); ylabel('t'); title('staggered grid')
%     legend('velocity', 'stress');
end

% u(x,t) = sin(2pix)sin(2pit);
v_ex = @(t,x) 2*pi*cos(2*pi*t)*sin(2*pi*x);      %u_t
sig_ex = @(t,x) mu*2*pi*cos(2*pi*x)*sin(2*pi*t); %mu u_x

v     = NaN(length(tv),length(xv));
sigma = NaN(length(tv)+1,length(xv));

v(1,:)     = v_ex(-dt/2,xv);               % v(x,-dt/2) blue points
sigma(2,2:2:nx) = sig_ex(0,xv(2:2:nx));    % sigma(x,0) red  points

for j = 3 : nt + 1
    if mod(j,2) ~= 0
        for i = 3 : 2 : nx-2
            v(j,i) = v(j-2,i) + dt/dx*1/rho*(sigma(j-1,i+1) - sigma(j-1,i-1));
        end
        %boundary conditions for v
        v(j,1)   = v_ex((j-2)*dt/2,0);
        v(j,end) = v_ex((j-2)*dt/2,L);

    else
        for i = 2 : 2 : nx-1
            sigma(j,i) = sigma(j-2,i) + dt/dx*mu*(v(j-1,i+1) - v(j-1,i-1));
        end
    end

end


% Reconstruct velocity at blue points from the total matrix
v_grid = [];
for i = 2 : nt
    vec = v(i,~isnan(v(i,:)));
    if ~isempty(vec)
        v_grid = [v_grid; vec];
    end
end

% Reconstruct stress at red points from the total matrix
sig_grid = [];
for i = 1 : nt+1
    vec = sigma(i,~isnan(sigma(i,:)));
    if ~isempty(vec)
        sig_grid = [sig_grid; vec];
    end
end


% surface plot of velocity and stress
figure(3)
subplot(1,2,1)
surf(Xv,Yv,v_grid,'EdgeColor','none');
xlabel('x'); ylabel('t'); title('v = \rho u_t')
subplot(1,2,2)
surf(Xs,Ys,sig_grid,'EdgeColor','none');
xlabel('x'); ylabel('t'); title('\sigma = \mu u_x')


% plot of velocity and stress at final time
% comparison between analytical and approximated solutions
figure(2)
subplot(1,2,1)
plot(xv(2:2:end-1),sigma(end,(2:2:end-1))); hold on;
plot(xv,sig_ex(T,xv));
legend('\sigma_h','\sigma_{ex}');
xlabel('x'); title('\sigma = \mu u_x at T');
ylim([-2*pi,2*pi]); grid on;
subplot(1,2,2)
plot(xv(1:2:end),v(end,(1:2:end))); hold on;
plot(xv,v_ex(T-dt/2,xv));
legend('v_h','v_{ex}');
xlabel('x'); title('v = \rho u_t at T');
ylim([-2*pi,2*pi]); grid on;




