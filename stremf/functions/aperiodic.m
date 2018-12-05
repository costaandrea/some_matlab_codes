function newM = aperiodic(mat,direc,varargin)

if nargin==2
    hm=1;
elseif nargin==3
    hm=varargin{1};
else
    error('wrong input number')
end

if direc==1
    newM = mat(hm+1:end-hm,:);

elseif direc==2
    
    newM = mat(:,hm+1:end-hm);
    
else
    error(' invalid dimension value')
end