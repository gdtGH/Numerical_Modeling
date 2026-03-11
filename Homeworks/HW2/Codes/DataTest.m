%% DataTest — test-case factory
%  Returns the Data struct for the requested test.
%  Called by RunMain.m (legacy) and RunMainSWE.m (HW2).
%
%  HW2 test cases:
%    'HW2_P4'   Travelling-wave, periodic BC  (Point 4)
%    'HW2_P5a'  Reflecting walls,  wall  BC   (Point 5a)
%    'HW2_P5b'  Travelling-wave + friction,   (Point 5b)
%
function [Data] = DataTest(TestName)

% =========================================================================
%  LEGACY TEST CASES  (HW1 wave equation, used by RunMain.m)
% =========================================================================

if strcmp(TestName,'Test1')

    Data.name    = TestName;
    Data.domain  = [0, 2*pi];
    Data.boundary = 'AA';
    Data.T  = 2;
    Data.dt = 0.0001;
    Data.p  = 1;

    Data.k   = 1;
    Data.mu  = 1;
    Data.ro  = 1;
    Data.c   = sqrt(Data.mu/Data.ro);
    Data.a   = 1;
    Data.b   = 1;
    Data.alfa = 1;
    Data.w   = Data.k * Data.c;

    Data.force = @(x,t) (Data.k^2 - Data.w^2) .* cos(Data.k*x - Data.w*t);
    Data.gD1   = @(t) cos(Data.w*t);
    Data.gD2   = @(t) cos(Data.k*Data.domain(2) - Data.w*t);
    Data.gN1   = @(t)  Data.mu*Data.k.*sin(Data.w*t);
    Data.gN2   = @(t) -Data.mu*Data.k.*sin(Data.k*Data.domain(2) - Data.w*t);
    Data.gR1   = @(t)  Data.k.*sin(Data.w*t) - Data.a*cos(Data.w*t);
    Data.gR2   = @(t) -Data.k*sin(Data.k*Data.domain(2) - Data.w*t) ...
                      + Data.b*cos(Data.k*Data.domain(2) - Data.w*t);
    Data.gI1   = @(t) -(Data.alfa*Data.w + Data.c*Data.k).*sin(Data.w*t);
    Data.gI2   = @(t)  (Data.alfa*Data.w - Data.c*Data.k) ...
                          .*sin(Data.k*Data.domain(2) - Data.w*t);

    Data.u0 = @(x) cos(Data.k*x);
    Data.v0 = @(x) Data.w*sin(Data.k*x);

    Data.uex     = @(x,t) cos(Data.k*x - Data.w*t);
    Data.graduex = @(x,t) -Data.k*sin(Data.k*x - Data.w*t);

elseif strcmp(TestName,'Test2')

    Data.name    = TestName;
    Data.domain  = [0, 2];
    Data.boundary = 'AA';
    Data.T  = 2;
    Data.dt = 0.0001;
    Data.p  = 1;

    Data.k   = 1;
    Data.mu  = 1;
    Data.ro  = 1;
    Data.c   = sqrt(Data.mu/Data.ro);
    Data.a   = 1;
    Data.b   = 1;
    Data.alfa = 0.9;
    Data.w   = Data.k * Data.c;

    Data.force = @(x,t) 0.*t.*x;
    Data.gD1   = @(t) 0.*t;
    Data.gD2   = @(t) 0.*t;
    Data.gN1   = @(t) 0.*t;
    Data.gN2   = @(t) 0.*t;
    Data.gR1   = @(t) 0.*t;
    Data.gR2   = @(t) 0.*t;
    Data.gI1   = @(t) 0.*t;
    Data.gI2   = @(t) 0.*t;

    Data.u0 = @(x) exp(-100*(x-1).^2);
    Data.v0 = @(x) 0.*x;

    Data.uex     = @(x,t) 0*x.*t;
    Data.graduex = @(x,t) 0.*x.*t;

% =========================================================================
%  HW2 — SHALLOW WATER EQUATIONS (SEM + theta-method)
% =========================================================================

elseif strcmp(TestName,'HW2_P4')
    %----------------------------------------------------------------------
    %  Point 4 — Right-travelling Gaussian, periodic BC.
    %  Exact solution:  eta(x,t) = eta0( mod(x - c*t , L) )
    %                     q(x,t) = c * eta(x,t)
    %----------------------------------------------------------------------
    Data.name     = 'HW2_P4';
    Data.domain   = [0, 1];
    Data.boundary = 'PP';           % periodic–periodic
    Data.T        = 0.5;
    Data.dt       = 1e-3;           % may be overridden in RunMainSWE

    Data.p        = 4;              % SEM polynomial degree (default)

    % Physical parameters
    Data.H = 1;
    Data.g = 9.81;
    Data.c = sqrt(Data.g * Data.H); % wave speed
    Data.L = Data.domain(2) - Data.domain(1);

    % Initial conditions
    Data.eta0 = @(x) exp(-50*(x - 0.5).^2);
    Data.q0   = @(x) Data.c .* exp(-50*(x - 0.5).^2);

    % Exact solution at time t (periodic wrap via mod)
    Data.uex     = @(x,t) exp(-50*(mod(x - Data.c*t, Data.L) - 0.5).^2);
    Data.graduex = @(x,t) -100*(mod(x - Data.c*t, Data.L) - 0.5) ...
                          .* exp(-50*(mod(x - Data.c*t, Data.L) - 0.5).^2);

    % Zero forcing
    Data.force = @(x,t) 0*x;

    % Friction flag (off for P4)
    Data.use_friction = 0;
    Data.gamma        = 0;

    % These flags are set in RunMainSWE, NOT here
    % Data.calc_errors, Data.visual_graph, Data.snapshot, Data.surf

elseif strcmp(TestName,'HW2_P5a')
    %----------------------------------------------------------------------
    %  Point 5a — Reflecting walls.
    %  BC:  q(0,t) = q(L,t) = 0   (Dirichlet, no-flux walls)
    %       eta: natural BC (free at walls)
    %  Exact solution: not available in closed form.
    %----------------------------------------------------------------------
    Data.name     = 'HW2_P5a';
    Data.domain   = [0, 1];
    Data.boundary = 'WW';           % wall–wall (reflecting)
    Data.T        = 0.5;
    Data.dt       = 1e-3;

    Data.p = 4;

    Data.H = 1;
    Data.g = 9.81;
    Data.c = sqrt(Data.g * Data.H);
    Data.L = Data.domain(2) - Data.domain(1);

    % Initial conditions: centred Gaussian at rest
    Data.eta0 = @(x) exp(-50*(x - 0.5).^2);
    Data.q0   = @(x) 0*x;            % q=0 initially (compatible with wall BC)

    % No exact solution
    Data.uex     = @(x,t) 0*x;
    Data.graduex = @(x,t) 0*x;

    Data.force = @(x,t) 0*x;

    Data.use_friction = 0;
    Data.gamma        = 0;

elseif strcmp(TestName,'HW2_P5b')
    %----------------------------------------------------------------------
    %  Point 5b — Friction term  gamma * q  in the momentum equation.
    %  q_t + gH eta_x + gamma*q = 0
    %  Periodic BC, same initial condition as P4.
    %  No closed-form exact solution with friction.
    %----------------------------------------------------------------------
    Data.name     = 'HW2_P5b';
    Data.domain   = [0, 1];
    Data.boundary = 'PP';
    Data.T        = 0.5;
    Data.dt       = 1e-3;

    Data.p = 4;

    Data.H = 1;
    Data.g = 9.81;
    Data.c = sqrt(Data.g * Data.H);
    Data.L = Data.domain(2) - Data.domain(1);

    Data.eta0 = @(x) exp(-50*(x - 0.5).^2);
    Data.q0   = @(x) Data.c .* exp(-50*(x - 0.5).^2);

    Data.uex     = @(x,t) 0*x;   % no exact solution with friction
    Data.graduex = @(x,t) 0*x;

    Data.force = @(x,t) 0*x;

    % Friction
    Data.use_friction = 1;
    Data.gamma        = 1;        % friction coefficient

else
    error('DataTest: unknown test name ''%s''.', TestName);
end

end % function DataTest