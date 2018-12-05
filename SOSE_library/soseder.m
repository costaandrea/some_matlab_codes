function [dfdx, dx, dfdy, dy] = soseder(x,y,f,symm,ord,varargin)


if nargin==5

    [dx,dy] = ll2dist(x,y);

elseif nargin>5
    
    dx = varargin{1};
    dy = varargin{2};
    
end


[M,N] = size(f);

dfx = NaN(M,N);
dfy = NaN(M,N);


if ord==1 %FIRST DERIVATIVE

ii = 2:M-1;
jj = 2:N-1;

%%Take second order central derivative on interior points    
    dfx(ii,:) = (f(ii+1,:) - f(ii-1,:));
    dfy(:,jj) = (f(:,jj+1) - f(:,jj-1));
    
    %%% Derivative %%%
    dfdx(ii,:) = dfx(ii,:) ./ (2*dx(ii,:));
    dfdy(:,jj) = dfy(:,jj) ./ (2*dy(:,jj));   

    
  if strcmp(symm,' ') || strcmp(symm,'-')
%%Take second order forward derivative on first point    
    dfx(1,:) = (-f(1+2,:) + 4*f(1+1,:) - 3*f(1,:));
    dfy(:,1) = (-f(:,1+2) + 4*f(:,1+1) - 3*f(:,1));
    
    %%% Derivative %%%
    dfdx(1,:) = dfx(1,:) ./ (2*dx(1,:));
    dfdy(:,1) = dfy(:,1) ./ (2*dy(:,1));   
    
%%Take second order backward derivative on last point
    dfx(M,:) = (3*f(M,:) - 4*f(M-1,:) + f(M-2,:));
    dfy(:,N) = (3*f(:,N) - 4*f(:,N-1) + f(:,N-2));

%%% Derivative %%%
    dfdx(M,:) = dfx(M,:) ./ (2*dx(M,:));
    dfdy(:,N) = dfy(:,N) ./ (2*dy(:,N));   

    
  elseif strcmp(symm,'x')
%%x boundary (second order central first derivative)
    dfx(1,:) = (f(2,:) - f(M,:));
    dfx(M,:) = (f(1,:) - f(M-1,:));
%%Take second order forward derivative on first y point    
    dfy(:,1) = (-f(:,1+2) + 4*f(:,1+1) - 3*f(:,1));
%%Take second order backward derivative on last y point
    dfy(:,N) = (3*f(:,N) - 4*f(:,N-1) + f(:,N-2));
     
%%% Derivative %%%
    dfdx(1,:) = dfx(1,:) ./ (2*dx(1,:));
    dfdx(M,:) = dfx(M,:) ./ (2*dx(M,:));
    dfdy(:,1) = dfy(:,1) ./ (2*dy(:,1));       
    dfdy(:,N) = dfy(:,N) ./ (2*dy(:,N));       
    
  end
 
  
  
elseif ord==2 %SECOND DERIVATIVE

ii = 2:M-1;
jj = 2:N-1;
    
%%Take second order central derivative on interior points    
    dfx(ii,:) = (f(ii+1,:) - 2*f(ii,:) + f(ii-1,:));
    dfy(:,jj) = (f(:,jj+1) - 2*f(:,jj) + f(:,jj-1));
    
    %%% Derivative %%%
    dfdx(ii,:) = dfx(ii,:) ./ (dx(ii,:).^2);
    dfdy(:,jj) = dfy(:,jj) ./ (dy(:,jj).^2);   
    
    
  if strcmp(symm,' ') || strcmp(symm,'-')
%%Take second order forward derivative on first point    
    dfx(1,:) = (-f(1+3,:) + 4*f(1+2,:) - 5*f(1+1,:) + 2*f(1,:));
    dfy(:,1) = (-f(:,1+3) + 4*f(:,1+2) - 5*f(:,1+1) + 2*f(:,1));
    
    %%% Derivative %%%
    dfdx(1,:) = dfx(1,:) ./ (dx(1,:).^2);
    dfdy(:,1) = dfy(:,1) ./ (dy(:,1).^2);   
    
%%Take second order backward derivative on last point
    dfx(M,:) = (2*f(M,:) - 5*f(M-1,:) + 4*f(M-2,:) - f(M-3,:));
    dfy(:,N) = (2*f(:,N) - 5*f(:,N-1) + 4*f(:,N-2) - f(:,N-3));

    %%% Derivative %%%
    dfdx(M,:) = dfx(M,:) ./ (dx(M,:).^2);
    dfdy(:,N) = dfy(:,N) ./ (dy(:,N).^2);   

    
  elseif strcmp(symm,'x')
%%x boundary (second order central second derivative)
    dfx(1,:) = (f(1+1,:) - 2*f(1,:) + f(M,:)); 
    dfx(M,:) = (f(1,:) - 2*f(M,:) + f(M-1,:));
%%Take second order forward derivative on first y point    
    dfy(:,1) = (-f(:,1+3) + 4*f(:,1+2) - 5*f(:,1+1) + 2*f(:,1));
%%Take second order backward derivative on last y point
    dfy(:,N) = (2*f(:,N) - 5*f(:,N-1) + 4*f(:,N-2) - f(:,N-3));
     
    %%% Derivative %%%
    dfdx(1,:) = dfx(1,:) ./ (dx(1,:).^2);
    dfdx(M,:) = dfx(M,:) ./ (dx(M,:).^2);
    dfdy(:,1) = dfy(:,1) ./ (dy(:,1).^2);       
    dfdy(:,N) = dfy(:,N) ./ (dy(:,N).^2);       
    
  end

 
   
elseif ord==3 %THIRD DERIVATIVE

ii = 3:M-2;
jj = 3:N-2;

%%Take second order central derivative on interior points    
    dfx(ii,:) = (f(ii+2,:) - 2*f(ii+1,:) + 2*f(ii-1,:) - f(ii-2,:));
    dfy(:,jj) = (f(:,jj+2) - 2*f(:,jj+1) + 2*f(:,jj-1) - f(:,jj-2));
    
    %%% Derivative %%%
    dfdx(ii,:) = dfx(ii,:) ./ (2*dx(ii,:).^3);
    dfdy(:,jj) = dfy(:,jj) ./ (2*dy(:,jj).^3);   

    
  if strcmp(symm,' ') || strcmp(symm,'-')
%%Take second order forward derivative on first 2 points    
    dfx(1,:) = (-3*f(1+4,:) + 14*f(1+3,:) - 24*f(1+2,:) +18*f(1+1,:) - 5*f(1,:)); 
    dfy(:,1) = (-3*f(:,1+4) + 14*f(:,1+3) - 24*f(:,1+2) +18*f(:,1+1) - 5*f(:,1));
    dfx(2,:) = (-3*f(2+4,:) + 14*f(2+3,:) - 24*f(2+2,:) +18*f(2+1,:) - 5*f(2,:)); 
    dfy(:,2) = (-3*f(:,2+4) + 14*f(:,2+3) - 24*f(:,2+2) +18*f(:,2+1) - 5*f(:,2));
    
    %%% Derivative %%%
    dfdx(1,:) = dfx(1,:) ./ (2*dx(1,:).^3);
    dfdy(:,1) = dfy(:,1) ./ (2*dy(:,1).^3);   
    dfdx(2,:) = dfx(2,:) ./ (2*dx(2,:).^3);
    dfdy(:,2) = dfy(:,2) ./ (2*dy(:,2).^3);   
 
%%Take second order backward derivative on last 2 points
    dfx(M,:) = (5*f(M,:) - 18*f(M-1,:) + 24*f(M-2,:) - 14*f(M-3,:) + 3*f(M-4,:));
    dfy(:,N) = (5*f(:,N) - 18*f(:,N-1) + 24*f(:,N-2) - 14*f(:,N-3) + 3*f(:,N-4));
    dfx(M-1,:) = (5*f(M-1,:) - 18*f(M-2,:) + 24*f(M-3,:) - 14*f(M-4,:) + 3*f(M-5,:));
    dfy(:,N-1) = (5*f(:,N-1) - 18*f(:,N-2) + 24*f(:,N-3) - 14*f(:,N-4) + 3*f(:,N-5));
    
    %%% Derivative %%%
    dfdx(M,:) = dfx(M,:) ./ (2*dx(M,:).^3);
    dfdy(:,N) = dfy(:,N) ./ (2*dy(:,N).^3);   
    dfdx(M-1,:) = dfx(M-1,:) ./ (2*dx(M-1,:).^3);
    dfdy(:,N-1) = dfy(:,N-1) ./ (2*dy(:,N-1).^3);   

    
  elseif strcmp(symm,'x')
%%x boundary (second order central third derivative)
    dfx(1,:) = (f(1+2,:) - 2*f(1+1,:) + 2*f(M,:) - f(M-1,:));
    dfx(M,:) = (f(2,:) - 2*f(1,:) + 2*f(M-1,:) - f(M-2,:)); 
    dfx(2,:) = (f(2+2,:) - 2*f(2+1,:) + 2*f(1,:) - f(M,:)); 
    dfx(M-1,:) = (f(1,:) - 2*f(M,:) + 2*f(M-2,:) - f(M-3,:)); 
%%Take second order forward derivative on first 2 y point    
    dfy(:,1) = (f(:,1+3) - 3*f(:,1+2) + 3*f(:,1+1) - f(:,1));
    dfy(:,2) = (f(:,2+3) - 3*f(:,2+2) + 3*f(:,2+1) - f(:,2));
%%Take second order backward derivative on last 2 y point
    dfy(:,N) = (5*f(:,N) - 18*f(:,N-1) + 24*f(:,N-2) - 14*f(:,N-3) + 3*f(:,N-4));
    dfy(:,N-1) = (5*f(:,N-1) - 18*f(:,N-2) + 24*f(:,N-3) - 14*f(:,N-4) + 3*f(:,N-5));
    
    %%% Derivative %%%
    dfdx(1,:) = dfx(1,:) ./ (2*dx(1,:).^3); 
    dfdx(M,:) = dfx(M,:) ./ (2*dx(M,:).^3); 
    dfdx(2,:) = dfx(2,:) ./ (2*dx(2,:).^3); 
    dfdx(M-1,:) = dfx(M-1,:) ./ (2*dx(M-1,:).^3); 
    dfdy(:,1) = dfy(:,1) ./ (2*dy(:,1).^3);       
    dfdy(:,2) = dfy(:,2) ./ (2*dy(:,2).^3);       
    dfdy(:,N) = dfy(:,N) ./ (2*dy(:,N).^3);       
    dfdy(:,N-1) = dfy(:,N-1) ./ (2*dy(:,N-1).^3);       
    
  end
  
  
  
elseif ord==4 %FOURTH DERIVATIVE
   
ii = 5:M-4;
jj = 5:N-4;

%%Take second order central derivative on interior points    
    dfx(ii,:) = (f(ii+2,:) - 4*f(ii+1,:) + 6*f(ii,:) - 4*f(ii-1,:) + f(ii-2,:));
    dfy(:,jj) = (f(:,jj+2) - 4*f(:,jj+1) + 6*f(:,jj) - 4*f(:,jj-1,:) + f(:,jj-2));
    
    %%% Derivative %%%
    dfdx(ii,:) = dfx(ii,:) ./ (dx(ii,:).^4);
    dfdy(:,jj) = dfy(:,jj) ./ (dy(:,jj).^4);   
    
    
  if strcmp(symm,' ') || strcmp(symm,'-')
%%Take second order forward derivative on first 4 points    
    dfx((1:4),:) = (-2*f((1:4)+5,:) + 11*f((1:4)+4,:) - 24*f((1:4)+3,:) + 26*f((1:4)+2,:) - 14*f((1:4)+1,:) + 3*f((1:4),:)); 
    dfy(:,(1:4)) = (-2*f(:,(1:4)+5) + 11*f(:,(1:4)+4) - 24*f(:,(1:4)+3) + 26*f(:,(1:4)+2) - 14*f(:,(1:4)+1) + 3*f(:,(1:4)));
    
    %%% Derivative %%%
    dfdx((1:4),:) = dfx((1:4),:) ./ (dx((1:4),:).^4);
    dfdy(:,(1:4)) = dfy(:,(1:4)) ./ (dy(:,(1:4)).^4);   
  
%%Take second order backward derivative on last 4 points
    dfx((M-3:M),:) = (3*f((M-3:M),:) - 14*f((M-3:M)-1,:) + 26*f((M-3:M)-2,:) - 24*f((M-3:M)-3,:) + 11*f((M-3:M)-4,:) - 2*f((M-3:M)-5,:));
    dfy(:,(N-3:N)) = (3*f(:,(N-3:N)) - 14*f(:,(N-3:N)-1) + 26*f(:,(N-3:N)-2) - 24*f(:,(N-3:N)-3) + 11*f(:,(N-3:N)-4) - 2*f(:,(N-3:N)-5));
    
    %%% Derivative %%%
    dfdx((M-3:M),:) = dfx((M-3:M),:) ./ (dx((M-3:M),:).^4);
    dfdy(:,(N-3:N)) = dfy(:,(N-3:N)) ./ (dy(:,(N-3:N)).^4);   
  
  elseif strcmp(symm,'x')
%%x boundary (second order central fourth derivative)
    dfx(1,:) = (f(1+2,:) - 4*f(1+1,:) + 6*f(1,:) - 4*f(M,:) + f(M-1,:)); 
    dfx(2,:) = (f(2+2,:) - 4*f(2+1,:) + 6*f(2,:) - 4*f(2-1,:) + f(M,:)); 
    dfx((3:4),:) = (f((3:4)+2,:) - 4*f((3:4)+1,:) + 6*f((3:4),:) - 4*f((3:4)-1,:) + f((3:4)-2,:)); 
    dfx(M,:) = (f(2,:) - 4*f(1,:) + 6*f(M,:) - 4*f(M-1,:) + f(M-2,:)); 
    dfx(M-1,:) = (f(1,:) - 4*f(M,:) + 6*f(M-1,:) - 4*f(M-2,:) + f(M-3,:)); 
    dfx((M-3:M-2),:) = (f((M-3:M-2)+2,:) - 4*f((M-3:M-2)+1,:) + 6*f((M-3:M-2),:) - 4*f((M-3:M-2)-1,:) + f((M-3:M-2)-2,:)); 

%%Take second order forward derivative on first 4 y point    
    dfy(:,(1:4)) = (-2*f(:,(1:4)+5) + 11*f(:,(1:4)+4) - 24*f(:,(1:4)+3) + 26*f(:,(1:4)+2) - 14*f(:,(1:4)+1) + 3*f(:,(1:4)));

%%Take second order backward derivative on last 4 y point
    dfy(:,(N-3:N)) = (3*f(:,(N-3:N)) - 14*f(:,(N-3:N)-1) + 26*f(:,(N-3:N)-2) - 24*f(:,(N-3:N)-3) + 11*f(:,(N-3:N)-4) - 2*f(:,(N-3:N)-5));
    
    %%% Derivative %%%
    dfdx((1:4),:) = dfx((1:4),:) ./ (dx((1:4),:).^4); 
    dfdx((M-3:M),:) = dfx((M-3:M),:) ./ (dx((M-3:M),:).^4); 
    dfdy(:,(1:4)) = dfy(:,(1:4)) ./ (dy(:,(1:4)).^4);       
    dfdy(:,(N-3:N)) = dfy(:,(N-3:N)) ./ (dy(:,(N-3:N)).^4);       
    
  end

end 