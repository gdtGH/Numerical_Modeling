function [A,M]=Matrix1D(Data,femregion)
%% [A,M] = Matrix1D(Data,Femregion)
%==========================================================================
% Assembly of the stiffness matrices A and rhs f
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          femregion   : (struct)  see CreateFemregion.m
%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffnes matrix
%          M           : (sparse(ndof,ndof) real) mass matrix


fprintf('Assembling the matrices M and A... \n');


% connectivity infos
ndof         = femregion.ndof;         % degrees of freedom
nln          = femregion.nln;          % local degrees of freedom
ne           = femregion.ne;           % number of elements
connectivity = femregion.connectivity; % connectivity matrix


% shape functions
[basis] = ShapeBasis;

% quadrature nodes and weights for integrals
[nodes_1D, w_1D] = Quadrature(2);

% evaluation of shape bases on quadrature nodes
[Phi,GradPhi] = EvalShapeBasis(basis,nodes_1D);


% Assembly begin ...
A = sparse(ndof,ndof);  % Global Stiffness matrix
M = sparse(ndof,ndof);  % Global mass matrix

for ie = 1 : ne
     
    iglo = connectivity(1:nln,ie);
    [BJ, nodes_1D_phys] = GetJacobian(femregion.coord(iglo,:), nodes_1D);

    % Evaluate mu at quadrature nodes
    if isfield(Data, 'mu')
        mu_loc = Data.mu(nodes_1D_phys(:));
    else
        mu_loc = ones(length(nodes_1D), 1);
    end

    % Evaluate rho at quadrature nodes (PML: complex, else = 1)
    if isfield(Data, 'ro_func')
        ro_loc = Data.ro_func(nodes_1D_phys(:));
    else
        ro_loc = ones(length(nodes_1D), 1);
    end

    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%
    [A_loc] = Stiffness(GradPhi, w_1D, nln, BJ, mu_loc);
    A(iglo,iglo) = A(iglo,iglo) + A_loc;
    
    %=============================================================%
    % MASS MATRIX
    %=============================================================%
    [M_loc] = Mass(Phi, w_1D, nln, BJ, ro_loc);
    M(iglo,iglo) = M(iglo,iglo) + M_loc;   

end

% Add impedance term to stiffness matrix if right BC is Impedance
if strcmp(Data.boundary(2), 'I')
    mu_L = Data.mu(Data.domain(2));       % µ(L)
    beta = sqrt(Data.ro * mu_L);          % β = sqrt(ρ·µ(L))
    A(end,end) = A(end,end) + 1i * Data.omega * beta;
end


