clear; close all;

%% Path
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing
addpath SemLib


%% Data for Test
Data = DataTest('Test2');

%% Options
Data.surf = 1; % solution in space-time domain
Data.snapshot = 0; % movie of the solution 
Data.visual_graph = 1; % solution at final time T
Data.calc_errors = 0;  % error calculation


%% Main routine
% [err1,sol1,fem1,D1] = Main(Data,16);
% [err2,sol2,fem2,D2] = Main(Data,32);
% [err3,sol3,fem3,D3] = Main(Data,64);
[err4,sol4,fem4,D4] = Main(Data,128);
% 
%% Plot Errors
hVec   = [fem1.h, fem2.h, fem3.h, fem4.h]; 
eVecL2 = [err1.L2, err2.L2, err3.L2, err4.L2];
eVecH1 = [err1.H1, err2.H1, err3.H1, err4.H1];

hs = subplot(2,1,1);
loglog(hVec,hVec.^2,'-+b','Linewidth',2); hold on;
loglog(hVec,eVecL2,'-or','Linewidth',2);
legend(sprintf('h^%i',2),'||u-u_h||_{L^2}');
ylabel('L^2-error');
xlabel('h');
hs.FontSize = 12;

hs = subplot(2,1,2);
loglog(hVec,hVec,'-+b','Linewidth',2); hold on;
loglog(hVec,eVecH1,'-or','Linewidth',2);
legend(sprintf('h^%i',1),'||u-u_h||_{H^1}');
ylabel('H^1-error')
xlabel('h');
hs.FontSize = 12;


