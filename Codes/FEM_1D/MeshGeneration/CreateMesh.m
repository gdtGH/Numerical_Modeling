function [Region] = CreateMesh(Data,nEl)
%% [Region] = CreateMesh(Data, nEl)
%=========================================================================
% Creates regular mesh
%==========================================================================
%    INPUT:
%          Data    : (struct)  see DataTest.m
%          nEl     : (int)     Number of mesh elements  
%
%    OUTPUT:
%          Region  : (struct)  mesh structure


x0 = Data.domain(1);
xL = Data.domain(2);

%================================================
% Geometrical info
 nVert = nEl + 1;
 p = linspace(x0,xL,nVert);
 t = [(1:nVert-1)' (2:nVert)']';
 MeshSize = (xL-x0)./nEl;
%================================================

% Mesh data structure
Region = struct('dim',1,...
               'domain',Data.domain,...
               'h',MeshSize,...
               'nvert',nVert,...
               'ne',nEl,...
               'coord',p',...
               'boundary_points',[x0,xL],...
               'connectivity',t);
           
           
