function [femregion] = CreateFemregion(Data,Region) 
%% [femregion] = CreateFemregion(Data,Region)
%==========================================================================
% Creates finite element space
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Data        : (struct)  see DataTest.m
%          Region      : (struct)  see CreateMesh.m
%
%    OUTPUT:
%          femregion    : (struct) finite element space

fprintf('Creating finite element space ... \n');


nln = Data.p + 1;
        
bound_pts = ones(length(Region.coord),1);
bound_pts = find((abs(Region.coord(:,1) - Data.domain(1,1)) < 1.e-5) ...
    | (abs(Region.coord(:,1) - Data.domain(1,2)) < 1.e-5));


%==========================================================================
% COSTRUZIONE STRUTTURA FEMREGION
%==========================================================================
femregion=struct('fem',Data.p,...
                'domain',Region.domain,...
                'h', Region.h,...
                'nln',nln,...
                'ndof',length(Region.coord),...
                'ne',Region.ne,...
                'dof',Region.coord,...
                'nqn_1D',2,...
                'coord',Region.coord,...
                'connectivity',Region.connectivity,...
                'boundary_points',bound_pts);
            
            