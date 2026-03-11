%% Data for the following PDE

function [Data] = DataTest(TestName)

if strcmp(TestName,'HW1_P2')
    
    Data.name   = TestName;
    Data.domain = [0, 1];
    % N = Neumann in x=0 ; I = Impedance in x=L
    Data.boundary = 'NI';

    % Physical parameters
    Data.ro    = 1;
    Data.mu    = @(x) ones(size(x));   % µ(x) = 1
    Data.omega = 5*pi;
    Data.vel   = 1;   % vel=1 => omega^2/vel^2 = omega^2, consistent with hw system

    % Source and boundary data
    Data.force = @(x) 0.*x;              % f = 0 (derived in report)
    Data.gN1   = @(omega,ro,vel) -omega;  % gN = 5π (Neumann at x=0)
    Data.gA    = 1i;                     % gA = i  (Impedance at x=L)
    
    % Dummy gN2/gD (not used for NI, but needed to avoid errors)
    Data.gN2   = @(omega,ro,vel) 0;
    Data.gD1   = @(omega,ro,vel) 0;
    Data.gD2   = @(omega,ro,vel) 0;

    % Exact solution
    Data.uex     = @(x,omega,ro,vel) sin(omega.*x);
    Data.graduex = @(x,omega,ro,vel) omega.*cos(omega.*x);

elseif strcmp(TestName,'HW1_P3')

    Data.name   = TestName;
    Data.domain = [0, 1];
    % N = Neumann at x=0 ; I = Impedance at x=L
    Data.boundary = 'NI';

    % Physical parameters — µ piecewise constant, discontinuous at x=0.5
    Data.ro    = 1;
    Data.mu    = @(x) 4.*(x <= 0.5) + 1.*(x > 0.5);
    Data.omega = 5*pi;
    Data.vel   = 1;

    % Source and boundary data
    Data.force = @(x) 0.*x;              % f = 0
    Data.gN1   = @(omega,ro,vel) -1;     % gN = 1 (stored with sign convention)
    Data.gA    = 0;                       % gA = 0 (perfectly absorbing at x=L)

    % Dummy unused fields
    Data.gN2   = @(omega,ro,vel) 0;
    Data.gD1   = @(omega,ro,vel) 0;
    Data.gD2   = @(omega,ro,vel) 0;

    % No exact solution available
    Data.uex     = @(x,omega,ro,vel) 0.*x;
    Data.graduex = @(x,omega,ro,vel) 0.*x;

elseif strcmp(TestName,'HW1_P4')

    Data.name     = TestName;
    Data.domain   = [0, 1];
    % N = Neumann at x=0 (natural) ; D = Dirichlet at x=L
    Data.boundary = 'ND';

    % Physical parameters — homogeneous medium
    Data.ro    = 1;
    Data.mu    = @(x) ones(size(x));
    Data.omega = 1;    % dummy, not used in eigenvalue analysis
    Data.vel   = 1;

    % All sources = 0
    Data.force = @(x) 0.*x;
    Data.gN1   = @(omega,ro,vel) 0;   % homogeneous Neumann at x=0
    Data.gN2   = @(omega,ro,vel) 0;
    Data.gD1   = @(omega,ro,vel) 0;
    Data.gD2   = @(omega,ro,vel) 0;   % homogeneous Dirichlet at x=L

    % No exact solution needed
    Data.uex     = @(x,omega,ro,vel) 0.*x;
    Data.graduex = @(x,omega,ro,vel) 0.*x;

elseif strcmp(TestName,'HW1_P5')

    L_phys    = 1;
    L_PML     = 0.2;
    sigma_max = 100;   % damping strength

    Data.name   = TestName;
    Data.domain = [0, L_phys + L_PML];
    % N = Neumann at x=0 ; D = Dirichlet at x=L+L_PML
    Data.boundary = 'ND';

    Data.ro    = 1;
    Data.omega = 5*pi;
    Data.vel   = 1;

    % Damping profile σ(x) — parabolic, zero in physical domain
    sigma  = @(x) sigma_max .* ((x - L_phys)./L_PML).^2 .* (x > L_phys);

    % Complex stretching s(x) = 1 - i*σ(x)/ω
    s_func = @(x) 1 - 1i .* sigma(x) ./ (5*pi);

    % Physical µ: same as P3 on [0,1], extended as µ=1 in PML
    mu_phys = @(x) 4.*(x <= 0.5) + 1.*(x > 0.5);

    % PML-modified coefficients:
    % µ̃(x) = µ(x)/s(x)   [complex stiffness]
    % ρ̃(x) = ρ(x)*s(x)   [complex density, ρ=1 so ρ̃=s(x)]
    Data.mu      = @(x) mu_phys(x) ./ s_func(x);
    Data.ro_func = @(x) s_func(x);

    % Source and BCs
    Data.force = @(x) 0.*x;
    Data.gN1   = @(omega,ro,vel) -1;   % same excitation as P3
    Data.gN2   = @(omega,ro,vel) 0;
    Data.gD1   = @(omega,ro,vel) 0;
    Data.gD2   = @(omega,ro,vel) 0;    % homogeneous Dirichlet at end of PML
    Data.gA    = 0;

    Data.uex     = @(x,omega,ro,vel) 0.*x;
    Data.graduex = @(x,omega,ro,vel) 0.*x;

end

           



