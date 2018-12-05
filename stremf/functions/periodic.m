function newM = periodic(mat,direc,varargin)
%periodic
% Use: newM = periodic(mat,direc,varargin)
% where varsgin can specify how many columns append/prepend (default is 1)


if nargin==2
    hm=1;
elseif nargin==3
    hm=varargin{1};
else
    error('wrong input number')
end


if hm==0 || hm==1 %%append one column/rows
if direc==1
    newM = [mat(end,:); mat; mat(1,:)];

elseif direc==2
    
    newM = [mat(:,end) mat mat(:,1)];
    
else
    error(' invalid dimension value')
end

else  %%append hm column/rows
    
    if direc==1
    newM = [mat(end-hm+1:end,:); mat; mat(1:hm,:)];

elseif direc==2
    
    newM = [mat(:,end-hm+1:end) mat mat(:,1:hm)];
    
else
    error(' invalid dimension value')
    end

end