% FINITE VOLUME - BURGER's EQUATION
% CENTRAL FLUX 
clearvars; clc; close all;

% number of grid cells
N = 800;
% grid nodes
% x = linspace(-1,1,N+1);
x = linspace(-2,2,N+1);
% cell lenght
dx = x(2:end) - x(1:end-1);
% primitive of initial codition
u_int = -1/(2*pi) * cos(2*pi*x);
% u_int = 0.*(x>0) - x.*(x<0); 
%initialization of average values
u0 = (u_int(2:end)-u_int(1:end-1))./dx;
figure(1);
plot([x(1:end-1); x(2:end)], [u0; u0], 'LineWidth',2);
% grid on; hold on;
% xp = [0:0.01:1]; yp = sin(2*pi*xp);
% plot(xp,yp,'k','LineWidth',1.5);
u = u0 + 0.1;


dt = 0.01;
T = 1;
t_end = 0;

% for final plot
uplot = NaN(N,T/dt);
uplot(:,1) = u';
j = 2;

for t = dt : dt : T
    % Runge-Kutta scheme for the ode du/dt = f(t,u)
    [t,u] = ode45(@ddtFiniteVolume,[0,dt],u);
    u = u(end,:);
    uplot(:,j) = u;
    plot([x(1:end-1); x(2:end)], [u;u], 'k');
    ylim([-1.5,1.5]);
    xlabel('x'); ylabel('u(x,t)'); 
    t_end = t_end + dt;
    title(['Burgers equation at t = ', num2str(t_end)]);
    drawnow;
    j = j + 1;    
end

figure(2)
xplot = x(1) - dx(1)/2 : dx(1) : x(end) - dx(1)/2;
tplot = 0:dt:T;
surf(tplot,xplot,uplot,'EdgeColor','none')
xlabel('t'); ylabel('x'); title('u(x,t)'); colorbar;


function [dudt] = ddtFiniteVolume(t,u)
    N = length(u);
    x = linspace(0,1,N+1);
    dx = x(2:end) - x(1:end-1);
    f = u.^2/2;   %-> bugers flux (modify here for other fluxes)
    dfdu = u;
    % f = u.^2./(u.^2 + (1-u).^2);
    % dfdu = (2.*u.*(u.^2 + (1-u).^2) - u.^2.*(2.*u - 2*(1-u)))/(u.^2 + (1-u).^2).^2;
    %% contant flux reconstruction
    % f_int = (f(1:end-1) + f(2:end))/2;  %constant flux reconstruction
    %% upwind
    speed =  (f(1:end-1) - f(2:end))./ (u(1:end-1) - u(2:end)); % speed of shock if u_left different form u_right
    special_case = abs((u(1:end-1) - u(2:end))) < 1E-8;
    speed(special_case) = dfdu(special_case);
    fL = f(1:end-1);
    fR = f(2:end);
    f_int = fL .* (speed > 0) + fR .* (speed < 0); % upwind
    %% Gudunov
    % ul_gt_ur = u(1:end-1) - u(2:end) > 0;
    % fL = f(1:end-1);
    % fR = f(2:end);
    % f_int = max(fL,fR).* ul_gt_ur + min(fL,fR) .* (~ul_gt_ur);  


    f_int = [0; f_int; 0];
    dudt = (f_int(1:end-1)-f_int(2:end))./dx';
end
