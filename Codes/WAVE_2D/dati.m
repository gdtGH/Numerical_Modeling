%% Test case 1 with analytical solution
% v_c = 1;                       %Velocity c
% uex     = @(x,y,t)      sin(4*pi*x).*sin(4*pi*y).*sin(2*pi*t);   %Analytical solution u
% uex_x   = @(x,y,t) 4*pi*cos(4*pi*x).*sin(4*pi*y).*sin(2*pi*t);   %Analytical gradient u_x
% uex_y   = @(x,y,t) 4*pi*sin(4*pi*x).*cos(4*pi*y).*sin(2*pi*t);   %Analytical gradient u_y
% 
% u0ex  = @(x,y)      0*sin(4*pi*x).*sin(4*pi*y);             %Initial condition u
% u1ex  = @(x,y)   2*pi*sin(4*pi*x).*sin(4*pi*y);             %Initial condition u_t
% 
% ff   = @(x,y,t) sin(2*pi*t).*4*pi^2.*sin(4*pi*x)...
%                              .*sin(4*pi*y)*(8*v_c^2 - 1);   %Forcing term space
%            
% 
% g    = @(x,y,t)   sin(4*pi*x).*sin(4*pi*y).*sin(2*pi*t);                      %Dirichlet condition
% h    = @(x,y,t)  sin(2*pi*t).*[-4*pi*sin(4*pi*x).*cos(4*pi*y); ...  %bottom
%                                 4*pi*cos(4*pi*x).*sin(4*pi*y);  ...  %right
%                                 4*pi*sin(4*pi*x).*cos(4*pi*y); ...   %top
%                                -4*pi*cos(4*pi*x).*sin(4*pi*y)];     %left    %Neumann condition
% 
%                         
% % Domain's bounds and boundary conditions
% % Omega = (xa,xb) x (ya,yb)
% xa = 0; xb = 1;
% ya = 0; yb = 1;
% 
% % boundary conditions on the sides of \partial\Omega
% % numeration is 1-bottom, 2-right, 3-up, 4-left
% % n = neumann, d = dirichlet
% cb = ['nnnn'];
% 
% 
% % Discretization parameters
% nex = 20;           % number of elements in the x direction 
% ney = 20;           % number of elements in the y direction 
% nx = 4;             % polynomial degree in each element along x-direction
% ny = 4;             % polynomial degree in each element along y-direction
% 
% t = 0;               %Initial time
% dt = 0.01;          %Time step 
% T = 1;               %Final time
%                         
                         
%% Test case 2 : channel with Dirichlet condition on the left edge
%  and sound hard boundary conditions on the other edges
% 
v_c = 1;                      %Velocity c

uex     = @(x,y,t) 0.*x.*y.*t;  %Analytical solution u
uex_x   = @(x,y,t) 0.*x.*y*t;   %Analytical gradient u_x
uex_y   = @(x,y,t) 0.*x.*y*t;   %Analytical gradient u_y

u0ex  = @(x,y) 0.*x.*y;       %Initial condition u
u1ex  = @(x,y) 0.*x.*y;       %Initial condition u_t

ff   = @(x,y,t) 0.*x.*y.*t;   %Forcing term

%Gauss wavelet
fp    = 5;                    % peak frequency [Hz]
t0    = 0.2;                  % delay [s]
alpha = 2*(pi*fp)^2;          % deviation
fac   = 2*pi*fp*sqrt(exp(1));   % scaling
%

g    = @(x,y,t) fac*(t-t0).*exp(-alpha*(t-t0).^2);  % Dirichlet condition

h    = @(x,y,t) [0.*x.*y.*t; 0.*x.*y.*t;...
                 0.*x.*y.*t; 0.*x.*y.*t];           % Neumann condition

% Domain's bounds and boundary conditions
% Omega = (xa,xb) x (ya,yb)
xa = 0; xb = 0.25;
ya = -0.10; yb = 0.10;

% boundary conditions on the sides of \partial\Omega
% numeration is 1-bottom, 2-right, 3-up, 4-left
% n = neumann, d = dirichlet
cb = ['nnnd'];


% Discretization parameters
nex = 20;             % number of elements in the x direction 
ney = 10;             % number of elements in the y direction 
nx = 6;              % polynomial degree in each element along x-direction
ny = 6;             % polynomial degree in each element along y-direction

t = 0;              %Initial time
dt = 0.005;          %Time step 
T = 1;             %Final time


%% Don't modify the following parameters
param    = zeros(20,1);
param(1) = 1;       % 1=SEM-NI,  2= Patching
param(2) = 0;       % 0=no reordering, 1=CM ordering, 2=AMD ordering
param(3) = 1;       % 1 = solve linear system by Choleski fact.
param(4) = 1;       % computes errors
param(5) = 0;       % 0 exact norms, 1= discrete norms
param(6) = nx*2;    % nq for LG quadrature formulas
param(7) = 0;       % 0=absolute errors, 1=relative errors
param(8) = 1;       % 0 no plot, 1 surf
param(9) = (nx+1);  % nodes used to plot numerical solution

