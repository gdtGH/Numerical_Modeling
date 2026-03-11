function [Region] = CreateMesh(Data,nEl)
%% [Region] = C_create_mesh(Data, nEl)
%=========================================================================
% Creates regular mesh
%==========================================================================
%    called in Main.m
%
%    INPUT:
%          Data    : (struct)  see DataTest.m
%          nEl     : (int)    Number of mesh elements  
%
%    OUTPUT:
%          Region  : (struct) having fields: dimension
%                                            domain 
%                                            mesh size
%                                            number of vertices
%                                            number of elements
%                                            coordinates
%                                            boundary points
%                                            connectivity


x0 = Data.domain(1);
xL = Data.domain(2);

if Data.p == 1
    npdx = 2;
elseif Data.p == 2
    npdx = 3;
elseif Data.p == 3
    npdx = 4;
elseif Data.p == 4
    npdx = 5;
elseif Data.p == 5
    npdx = 6;
else
    disp('case not implemented')
end


%================================================
% Geometrical info
 MeshSize = (xL-x0)./nEl;
%================================================

i = 0; 
for ie = 1 : nEl
    xb_ie = x0 + ie*MeshSize;
    xa_ie = xb_ie - MeshSize;
    
    [xp,~] = xwlgl(npdx,xa_ie,xb_ie); 
    p(i+1:npdx+i) = xp;
    i = i + npdx ;
end

p = uniquetol(p,1.e-5);
nVert = size(p,2); 

if (Data.p == 1)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]']';
elseif (Data.p == 2)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]']';
elseif (Data.p == 3)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]' [4:npdx-1:nVert-npdx+4]']';
elseif (Data.p == 4)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]' [4:npdx-1:nVert-npdx+4]' [5:npdx-1:nVert-npdx+5]']';
elseif (Data.p == 5)
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]' [4:npdx-1:nVert-npdx+4]' [5:npdx-1:nVert-npdx+5]' [6:npdx-1:nVert-npdx+6]']';
else
    disp('case not implemented')
end

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
           
           
