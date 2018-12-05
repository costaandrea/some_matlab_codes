function cinds = conschunks(inds)
%conschunks(inds)
% divide an array INDS in chunks of consecutive elements.
% The output is a cell with the elements of the different 
% chunks in different cells.

[v,x] = find(diff(inds)>1); %find "jumps"

xx=[0 x length(inds)]; 

for ii=1:length(xx)-1
    
   cinds{ii} = inds(xx(ii)+1:xx(ii+1)); %output struct array
   
end