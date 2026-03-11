% Solution and error of the transport equation
%
%    u_t + a*u_x = 0
%
% x in [a,b]  and  t in [0,T] with the
% initial data  u(x,0) = u_0(x) with different methods
%   Forward Euler
%   Lax-Friedrichs
%   Lax-Wendroff
%   Upwind


clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial data
u0 = @(x) sin(8*pi*x);
% u0 = @(x) exp(-100*(x-0.5).^2);
% u0 = @(x) (x>=0.4).*(x<=0.6);



% velocity
a = @(x,t) 0.5.*(t<=pi/6) + 2.*(t>pi/6).*(t<=1) + 0.*x;
% a = @(x,t) -1 + 2.*t + 0.*x;
% a = @(x,t)  1 + 0.*t.*x;

% exact solution for computing the error
% u = @(x,t) sin(8*pi*(x-a(x,t).*t));
% u = @(x,t) 0.*x.*t;
u = @(x,t) u0(x-a(x,t)*t);


% Intervals
T = 1;
I = [0 1];

% Number of time steps
NT = 1000;  %dt = T/NT
% Number of space steps
NX = 500;   %h = b/NX 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% menu
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag = menu('Please select a method',...
%     'Forward-Euler/Central',...
%     'Lax-Friedrichs',...
%     'Lax-Wendroff',...
%     'Upwind');
% 
% disp('Please wait ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = T/NT;
dx = (I(2)-I(1))/NX;

x = I(1) : dx : I(2);
t = 0 : dt : T;

lambda = dt/dx;

% Initial conditions
SOL = zeros(NX+1,NT+1);
for j = 1:NX+1
    SOL(j,1) = u0(I(1) + (j-1)*dx);
end

%     figure;
for n = 1:NT
    for j = 2:NX
        
        if flag == 1 % FE scheme
            SOL(j,n+1) = SOL(j,n) - 0.5*lambda*a(x(j),t(n))*SOL(j+1,n)...
                                  + 0.5*lambda*a(x(j-1),t(n))*SOL(j-1,n);
        elseif flag == 2 % LF scheme
            SOL(j,n+1) = 0.5*(SOL(j+1,n)+SOL(j-1,n)) ...
                              - 0.5*lambda*a(x(j),t(n))*SOL(j+1,n)...
                              + 0.5*lambda*a(x(j-1),t(n))*SOL(j-1,n) ;
        elseif flag == 3 % LW scheme
            SOL(j,n+1) = SOL(j,n) - 0.5*lambda*a(x(j),t(n))*(SOL(j+1,n)-SOL(j-1,n)) ...
                              + 0.5*lambda^2*a(x(j),t(n)).^2*(SOL(j+1,n)-2*SOL(j,n)+SOL(j-1,n));
        elseif flag == 4 % UPWIND
            SOL(j,n+1) = SOL(j,n) - 0.5*lambda*a(x(j),t(n))*(SOL(j+1,n)-SOL(j-1,n)) ...
                              + 0.5*lambda*abs(a(x(j),t(n)))*(SOL(j+1,n)-2*SOL(j,n)+SOL(j-1,n));
                          
        end
    end
    SOL(1,n+1) = u(I(1),(n+1)*dt);
    % periodic condition
    % SOL(end,n+1) = SOL(1,n+1);
    % linear extrapolation
    SOL(end,n+1) = 2*SOL(end-1,n+1)-SOL(end-2,n+1);
end

% Computation of the L2-error for t=T
SOL_EX = zeros(NX+1,1);
for j = 1:NX+1
    SOL_EX(j,1) = u(I(1) + (j-1)*dx,NT*dt);
end
L2_ERR   = norm(SOL_EX-SOL(:,n+1),2)*dx^0.5
L1_ERR   = norm(SOL_EX-SOL(:,n+1),1)*dx
LINF_ERR = norm(SOL_EX-SOL(:,n+1),Inf)


% Visualization
time = ones(NX+1,1)*linspace(0,T,NT+1);
space = linspace(I(2),I(1),NX+1)'*ones(1,NT+1);
surf(time,space,SOL,'EdgeColor','none');
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view([-90 90]);






