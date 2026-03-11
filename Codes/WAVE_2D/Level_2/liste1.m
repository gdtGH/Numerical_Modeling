function [lbor,lint,lintint,lgamma]=liste1(ifro);
% LISTE1 Assembles lists of internal, boundary, interface nodes (similar to liste)
%
% [lbor,lint,lintint,lgamma]=liste1(ifro);
%
% It computes lists of internal, boundary, interface nodes.
% All lists are referred to global ordering on Omega
%
% Input : ifro = column array of length noe=nov(npdx*npdy,ne):
%            if (x_i,y_i) is internal to Omega then ifro(i)=0,
%            if (x_i,y_i) is on \partial\Omega then ifro(i)=1,
%
% Output:
% lbor: list of boundary nodes
% lint: list of internal nodes
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

noe=length(ifro);


lint=[];lbor=[];lgamma=[]; lintint=[];
for i=1:noe
if(ifro(i)==1)
lbor=[lbor;i];
else
lint=[lint;i];
end
end

for i=1:length(lint)
if(ifro(i)==-1)
lgamma=[lgamma;i];
else
lintint=[lintint;i];
end
end

return
