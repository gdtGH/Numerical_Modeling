% Solution and error of the Burger's equation
%
%    u_t + u u_x = 0
%
% x \in I=[a,b]  and  t \in [0,T] with the
% initial data  u(x,0) = u_0(x) 

clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial data
% u0 = @(x) 1./(1 + x.^2);
u0 = @(x)  0.25*x.*(x<=0.5) + 0.25*(1-x).*(x>0.5);
% u0 = @(x)  exp(-(x - 3).^2);
% u0 = @(x)  0 + 1.*(x>=2.5); 
% u0 = @(x)  0 + 1.*(x>=0.4).*(x<=0.6);
% u0 = @(x)  (x+1).*exp(-x/2);

% exact solution for computing the error
u = @(x,t) 0.*x.*t;

% Flux
F = @(x) 0.5*x.^2;

% Intervals
T = 10;
% I = [-1,1];
I = [0 2];
% I = [0 1];

% Number of time steps
NT = 2000;
% Number of space steps
NX = 1500;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% menu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag = menu('Please select a method',...
%     'Lax-Friedrichs')
%
% disp('Please wait ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = T/NT;
dx = (I(2)-I(1))/NX;
lambda = dt/dx;

% Initial conditions
SOL = zeros(NX+1,NT+1);
for j = 1:NX+1
    SOL(j,1) = u0(I(1) + (j-1)*dx);
end

% figure;
% x_space = linspace(I(1),I(2),NX+1);

for n = 1:NT
    for j = 2:NX
        
        if flag == 1 % LF scheme
            SOL(j,n+1) = 0.5*(SOL(j+1,n)+SOL(j-1,n)) ...
                              - 0.5*lambda*(F(SOL(j+1,n))-F(SOL(j-1,n)));
        elseif flag == 2 % LW scheme
            Ujph = 0.5 * (SOL(j,n) + SOL(j+1,n)) ...
                   - 0.5*lambda*(F(SOL(j+1,n)) - F(SOL(j,n))); 
            Ujmh = 0.5 * (SOL(j-1,n) + SOL(j,n)) ...
                   - 0.5*lambda*(F(SOL(j,n)) - F(SOL(j-1,n))); 
            
            SOL(j,n+1) = SOL(j,n) - lambda*(F(Ujph)-F(Ujmh));
                          
        end
    end
%     SOL(1,n+1) =  u(I(1),(n+1)*dt);
    
    %constant extrapolation
    SOL(1,n+1) = SOL(2,n+1);
    SOL(end,n+1) = SOL(end-1,n+1);

%         plot(x_space,SOL(:,n+1));
%         pause(0.15);
end

% close;
% Computation of the L2-error for t=T
SOL_EX = zeros(NX+1,1);
for j = 1:NX+1
    SOL_EX(j,1) = u(I(1) + (j-1)*dx,NT*dt);
end
L2_ERR   = norm(SOL_EX-SOL(:,n+1),2)*dx^0.5;
L1_ERR   = norm(SOL_EX-SOL(:,n+1),1)*dx;
LINF_ERR = norm(SOL_EX-SOL(:,n+1),Inf);


% Visualization
time = linspace(0,T,NT+1);
space = linspace(I(1),I(2),NX+1);
[SPACE, TIME] = meshgrid(space,time);
surf(time,space,SOL,'EdgeColor','none');
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view(2);






