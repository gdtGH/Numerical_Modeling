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

if strcmp(TestName,'HW2_P4_GAUSS')
    %----------------------------------------------------------------------
    % Point 4 — periodic travelling Gaussian pulse
    %----------------------------------------------------------------------
    Data.name     = 'HW2_P4_GAUSS';
    Data.domain   = [0, 1];
    Data.boundary = 'PP';
    Data.T        = 0.5;
    Data.dt       = 1e-3;
    Data.p        = 4;

    Data.H = 1;
    Data.g = 9.81;
    Data.c = sqrt(Data.g * Data.H);
    Data.L = Data.domain(2) - Data.domain(1);

    Data.test_label = 'Gaussian pulse';
    Data.test_tag   = 'gauss';

    % initial conditions
    Data.eta0 = @(x) exp(-50*(x - 0.5).^2);
    Data.q0   = @(x) Data.c .* exp(-50*(x - 0.5).^2);

    % exact solution
    Data.uex     = @(x,t) exp(-50*(mod(x - Data.c*t, Data.L) - 0.5).^2);
    Data.graduex = @(x,t) -100*(mod(x - Data.c*t, Data.L) - 0.5) ...
                          .* exp(-50*(mod(x - Data.c*t, Data.L) - 0.5).^2);

    % optional exact discharge
    Data.qex = @(x,t) Data.c .* Data.uex(x,t);

    Data.force = @(x,t) 0*x;

    Data.use_friction = 0;
    Data.gamma        = 0;

elseif strcmp(TestName,'HW2_P4_SMOOTH')
    %----------------------------------------------------------------------
    % Point 4 — periodic smooth travelling wave
    %----------------------------------------------------------------------
    Data.name     = 'HW2_P4_SMOOTH';
    Data.domain   = [0, 1];
    Data.boundary = 'PP';
    Data.T        = 0.5;
    Data.dt       = 1e-3;
    Data.p        = 4;

    Data.H = 1;
    Data.g = 9.81;
    Data.c = sqrt(Data.g * Data.H);
    Data.L = Data.domain(2) - Data.domain(1);

    Data.test_label = 'Smooth periodic wave';
    Data.test_tag   = 'smooth';

    % initial conditions
    Data.eta0 = @(x) cos(2*pi*x);
    Data.q0   = @(x) Data.c .* cos(2*pi*x);

    % exact solution
    Data.uex     = @(x,t) cos(2*pi*(x - Data.c*t));
    Data.graduex = @(x,t) -2*pi*sin(2*pi*(x - Data.c*t));

    % optional exact discharge
    Data.qex = @(x,t) Data.c .* Data.uex(x,t);

    Data.force = @(x,t) 0*x;

    Data.use_friction = 0;
    Data.gamma        = 0;

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