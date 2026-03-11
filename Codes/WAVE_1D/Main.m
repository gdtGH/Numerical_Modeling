function [Errors,Solutions,Femregion,Data] = Main(Data, nEl)
%%
%    INPUT:
%          Data    : (struct) Data struct
%          nEl     : (int)    Number of mesh elements
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) finite element space
%
%          Data        : (struct)  Data struct
%

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

[A_no_bc,M_no_bc] = Matrix1D(Data,Femregion);
if strcmp(Data.boundary,'RR')
    R = 0*A_no_bc;
    R(1,1) = Data.a;
    R(end,end) = Data.b;
    A_no_bc = A_no_bc + R;

elseif strcmp(Data.boundary,'AA')
    S = 0*A_no_bc;
    S(1,1) = Data.mu/Data.c * Data.alfa;
    S(end,end) = Data.mu/Data.c * Data.alfa;
end


%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================

[b_nbc] = Rhs1D(Data,Femregion);

%==========================================================================
% BUILD INITIAL CONDITIONS
%==========================================================================

x = Femregion.coord;
u0 = Data.u0(x);
v0 = Data.v0(x);

% [u0] = Snapshot(Femregion, u0, Data, 0);

u_snap = u0;

%% First step of leapfrog

rhs = 0.5*Data.dt^2*b_nbc(:,1) - (0.5*Data.dt^2*A_no_bc - M_no_bc)*u0 ...
    + Data.dt*M_no_bc*v0;

if (strcmp(Data.boundary,'AA'))
    rhs = rhs - 0.5*Data.dt^2*S*v0;
end

u_g = 0;

if (strcmp(Data.boundary,'DD') || strcmp(Data.boundary(1),'D') ...
        || strcmp(Data.boundary(2),'D'))
    [M,rhs, u_g] = BoundaryConditions(M_no_bc, rhs, Femregion, Data, 0);
else
    M = M_no_bc;
end

if strcmp(Data.boundary,'PP')
    M(1,:) = M(1,:) + M(end,:);
    M(end,:) = 0; M(end,1) = 1; M(end,end) = -1;
    rhs(1) = rhs(1) + rhs(end);
    rhs(end) = 0;
end

u1 = M\rhs;
u1 = u1 + u_g;



%% End of first step leapfrog

fprintf('Starting time-loop ... \n');

u_snap(:,2) = u1;
k = 2;

% Visualization of computational progress
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

for t = Data.dt : Data.dt : Data.T - Data.dt

    % Visualization of computational progress
    prog = ( 100*(t/Data.T) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);

    % General time step
    rhs = Data.dt^2*b_nbc(:,k) - (Data.dt^2*A_no_bc - 2*M_no_bc)*u1 - M_no_bc*u0;

    if (strcmp(Data.boundary,'AA'))
        rhs = rhs + 0.5*Data.dt*S*u0;
    end

    if (strcmp(Data.boundary,'DD') || strcmp(Data.boundary(1),'D') ...
            || strcmp(Data.boundary(2),'D'))
        [M,rhs, u_g] = BoundaryConditions(M_no_bc, rhs, Femregion, Data, t);
    elseif (strcmp(Data.boundary,'AA'))
        M = M_no_bc + 0.5*Data.dt*S;
    else
        M = M_no_bc;
    end

    if strcmp(Data.boundary,'PP')
        M(1,:) = M(1,:) + M(end,:);
        M(end,:) = 0; M(end,1) = 1; M(end,end) = -1;
        rhs(1) = rhs(1) + rhs(end);
        rhs(end) = 0;
    end

    u2 = M\rhs;
    u2 = u2 + u_g;

    u_snap(:,k+1) = u2;
    k = k +1;

    % Update the solution
    u0 = u1;
    u1 = u2;

end


%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

%Movie
if Data.snapshot
    t = 0;
    for i = 1 : size(u_snap,2)
        [~] = Snapshot(Femregion, u_snap(:,i), Data, t);
        % Put a pause between one step and the other to see the plot
        % pause(0.01);
        t = t + Data.dt;
    end
end

%Solution at final time T
[Solutions] = PostProcessing(Data,Femregion,u2);


%Surface plot
if Data.surf
    [xx,tt] = meshgrid(0:Data.dt:Data.T, Femregion.coord);
    figure(3);
    surf(xx,tt,u_snap,'EdgeColor','None');
    xlabel('time'); ylabel('x'); zlabel('u(x,t)');
end


%==========================================================================
% ERROR ANALYSIS
%==========================================================================
Errors = [];
if (Data.calc_errors)
    [Errors] = ComputeErrors(Data,Femregion,Solutions);
end



