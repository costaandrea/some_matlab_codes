function [islands, neigh, allisl, allc] = islands_finder(a,varargin)
%islands_finder finds the islands in a multiply connected domain.
% 
% [islands, neigh] = islands_finder(a,varargin)
% Feed in a mask with ones (on land) and zeros (ocean).
% An optiona flag (==4 or ==8) specifies how many neighbors of a cell you
% want. Default is ==8.
%
% Outputs: ISLANDS is a struct with ones on each island (one struct element
%                  per island).
%          NEIGH has ones on the closest cells to each island (one struct 
%                element per island).
%
%          ALLISL has ones on all the islands.
%
%          ALLC has ones on all cells next to coastlines.
%
% This requires the bfs_bord.m (breadth first search) and convolve.m scripts

%Example 1:
% a=zeros(6,3);
% a(3,1)=1;
% a(4,1)=1;
% a(4,2)=1;
% a(5,1)=1;
% a(5,2)=1;
% a(6,1)=1;
% a(1,2)=1;
% a
%
%Example 2:
% defiles
% detpath('netcdf_files')
% depthM=ncread('grid.nc','Depth');%depth matrix [m]
% 
% a = depthM;
% a(a==0)=1;
% a(a~=1)=0;



if nargin==2
    if varargin{1}==4
        mask = [0 1 0; 1 0 1; 0 1 0];
    elseif varargin{1}==8
        mask = ones(3);
    else
        error('Illegal flag value.')
    end
    
elseif nargin==1
    mask = ones(3);
end



b = bfs_bord(a); %find connected components (i.e., island)
clear a

kkk = setdiff(unique(b),0); %how many connected components

for ii=1:numel(kkk)
    c=zeros(size(b));
    
    ind=find(b == kkk(ii)); %indeces of the current connected component
    c(ind) = kkk(ii)*ones(1,numel(ind));

    islands{ii} = c/kkk(ii);
    
    c(convolve(c,mask)~=0) = 1;
    neigh{ii} = c - islands{ii}; %contours of the differen islands
end

allisl=b;
allisl(allisl~=0)=1;  %ones on all islands

tmp=convolve(allisl,mask);
tmp(tmp~=0)=1;
allc = tmp - allisl;   %ones next to all coastlines

%%%Exercise for the future user:
%%%what if islands have common neighbors?? (as in the first example above)