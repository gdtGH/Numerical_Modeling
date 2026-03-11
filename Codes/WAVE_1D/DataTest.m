%% Data for the following PDE
% ro*u_{tt} - (mu*u_{x})_x = f  in (a,b) x (0,T] 
% u(0) = u0;
% v(0) = v0;
% + b.c in {a,b} x (0,T]
function [Data] = DataTest(TestName)

if strcmp(TestName,'Test1')
    %% Input Data Test1
    
    Data.name = TestName;
    Data.domain = [0,2*pi];
    % D = Dirichlet; N = Neumann; R = Robin; P = periodic
    % A = Absorbing
    Data.boundary = 'AA';
    Data.T = 2;
    Data.dt = 0.0001; %time step
    
    % polynomial degree
    Data.p = 1;

    % Parameters and external forces 
    Data.k  = 1; 
    Data.mu = 1;
    Data.ro = 1;
    Data.c = sqrt(Data.mu/Data.ro);
    Data.a = 1;
    Data.b = 1;
    Data.alfa = 1;

    Data.w  = Data.k*Data.c;
    
    Data.force = @(x,t) (Data.k^2 - Data.w^2).* cos(Data.k*x - Data.w*t);
    Data.gD1  =  @(t) cos(Data.w*t);  % Dirichlet condition x = 0
    Data.gD2  =  @(t) cos(Data.k*Data.domain(2) - Data.w*t);  % Dirichlet condition x = 1
    Data.gN1  =  @(t) Data.mu*Data.k.*sin(Data.w*t);  % Neumann condition x = 0
    Data.gN2  =  @(t) -Data.mu*Data.k.*sin(Data.k*Data.domain(2) - Data.w*t);  % Dirichlet condition x = 1
    Data.gR1  =  @(t) Data.k.*sin(Data.w*t) - Data.a*cos(Data.w*t);  % Robin condition x = 0
    Data.gR2  =  @(t) -Data.k*sin(Data.k*Data.domain(2) - Data.w*t) ...
                      + Data.b*cos(Data.k*Data.domain(2) - Data.w*t);  % Robin condition x = 1
    
    Data.gI1  =  @(t) - (Data.alfa*Data.w + Data.c*Data.k).*sin(Data.w*t);
    Data.gI2  =  @(t)   (Data.alfa*Data.w - Data.c*Data.k)...
                          .*sin(Data.k*Data.domain(2) - Data.w*t);

    % Initial conditions
    Data.u0 = @(x) cos(Data.k*x);
    Data.v0 = @(x) Data.w*sin(Data.k*x);
    
    % Exact solution for error analysis
    Data.uex = @(x,t) cos(Data.k*x - Data.w*t);
    Data.graduex = @(x,t) - Data.k*sin(Data.k*x - Data.w*t);
    
elseif strcmp(TestName,'Test2')
    %% Input Data Test1
    
    Data.name = TestName;
    Data.domain = [0,2];
    % D = Dirichlet; N = Neumann; R = Robin; P = periodic
    % A = Absorbing
    Data.boundary = 'AA';
    Data.T = 2;
    Data.dt = 0.0001; %time step
    
    % polynomial degree
    Data.p = 1;

    % Parameters and external forces 
    Data.k  = 1; 
    Data.mu = 1;
    Data.ro = 1;
    Data.c = sqrt(Data.mu/Data.ro);
    Data.a = 1;
    Data.b = 1;
    Data.alfa = 0.9;

    Data.w  = Data.k*Data.c;
    
    Data.force = @(x,t) 0.*t.*x;
    Data.gD1  =  @(t) 0.*t;
    Data.gD2  =  @(t) 0.*t;
    Data.gN1  =  @(t) 0.*t;
    Data.gN2  =  @(t) 0.*t;
    Data.gR1  =  @(t) 0.*t;
    Data.gR2  =  @(t) 0.*t;
                      
    
    Data.gI1  =  @(t) 0.*t;
    Data.gI2  =  @(t) 0.*t;

    % Initial conditions
    Data.u0 = @(x) exp(-100*(x-1).^2);
    Data.v0 = @(x) 0.*x;
    
    % Exact solution for error analysis
    Data.uex = @(x,t) 0*x.*t;
    Data.graduex = @(x,t) 0.*x.*t;
    
end
