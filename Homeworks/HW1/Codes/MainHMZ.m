function [Errors,Solutions,Femregion,Data] = MainHMZ(Data, nEl)

fprintf('============================================================\n')
fprintf(['Solving test ', Data.name, ' with ',num2str(nEl),' elements \n']);

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = CreateMesh(Data,nEl);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[Femregion] = CreateFemregion(Data,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES and RIGHT-HAND SIDE
%==========================================================================

[A_nbc,M_nbc] = Matrix1D(Data,Femregion);


%==========================================================================
% BUILD FINITE ELEMENTS RHS
%==========================================================================

[b_nbc] = Rhs1D(Data,Femregion);


%==========================================================================
% SOLVE LINEAR SYSTEM - COMPUTE BOUNDARY CONDITIONS
%==========================================================================

% Solve the reduced system
K_nbc =  - A_nbc + Data.omega^2/Data.vel^2 * M_nbc;

if (strcmp(Data.boundary,'NN') || ...
    strcmp(Data.boundary,'NI') || ...
    strcmp(Data.boundary,'IN'))
    p = K_nbc\b_nbc;
else
    [K,b,u_g] = BoundaryConditions(K_nbc,b_nbc,Femregion,Data);
    p = K\b;
    p = p + u_g;
end

%[p] = Snapshot(Femregion, p, 'Pressure');

% Computation of the velocity 
n = length(p);
vel(2:n-1,1) = (p(3:n,1)-p(1:n-2,1))/(2*Femregion.h);
vel(1) = Data.gN1(Data.omega,Data.ro,Data.vel);
vel(n) = Data.gN2(Data.omega,Data.ro,Data.vel);
% dis = 1/(Data.ro*Data.omega).^2 * vel;

% Calcola soluzione esatta sui nodi FEM
x_nodes = Femregion.coord(:,1);
u_ex_nodes = Data.uex(x_nodes, Data.omega, Data.ro, Data.vel);

% Snapshot pressione
%if PlotOpts.show_snapshot
%    snap_title = PlotOpts.name_pressure;
%    if isempty(snap_title), snap_title = 'Pressure'; end
%    [p] = Snapshot(Femregion, p, snap_title, PlotOpts, u_ex_nodes);
%end

% Snapshot velocità — qui non abbiamo la velocità esatta in generale,
% quindi non passiamo u_ex (a meno che non la definisci in DataTest)
%if PlotOpts.show_velocity
%    snap_title = PlotOpts.name_velocity;
%    if isempty(snap_title), snap_title = 'Velocity'; end
%    [vel] = Snapshot(Femregion, vel, snap_title, PlotOpts);
%end

% Impedance in x=0
% Z_imp = p(1)/Data.gN1(Data.omega,Data.ro,Data.vel);

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[Solutions] = PostProcessing(Data, Femregion, p);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
Errors = [];
if (Data.calc_errors)
    [Errors] = ComputeErrors(Data,Femregion,Solutions);
end

%==========================================================================
% EIGENVALUE ANALYSIS
%==========================================================================

if Data.eig_analysis
    % eigenmodes(resonant frequencies) and eingevectors
    if (strcmp(Data.boundary,'NN') || ...
        strcmp(Data.boundary,'NI') || ...
        strcmp(Data.boundary,'IN'))
        M = M_nbc; A = A_nbc;
    else
        [M,~, ~] = BoundaryConditions(M_nbc, b_nbc, Femregion, Data);
        [A,~, ~] = BoundaryConditions(A_nbc, b_nbc, Femregion, Data);
    end
    
    
    % eigenvalues D -- eigenvector V
    % A x = k^2 M x
    [V,D] = eigs(A,M,6,0.1);
    
    k2 = diag(D);
    omega_h = sqrt(k2*Data.vel^2);
    modes = omega_h/(2*pi);
    kex = [0:5]'*pi/Data.domain(2);
    modes_ex = sqrt(kex.^2*Data.vel^2)/(2*pi);
    
    %% plot of 6 eigenvectors
    if Data.plot_eigvct
        for i = 1:6
            v = V(:,i);
            strTitle = sprintf('Eigenvector %d', i);
            [v] = Snapshot(Femregion, v, strTitle);
        end
    end
    
    disp([modes modes_ex]);
    
end

