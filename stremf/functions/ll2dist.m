function [dx,dy] = ll2dist(x,y)
%You input lon/lat matrices and I' ll give you the meters between the points

R = 6378137; %Earth radius [m]


[M,N] = size(x);

dx = NaN(M,N);
dy = NaN(M,N);

ii = 2:M-1;
jj = 2:N-1;


%%convert to radians
x=x*pi/180;
y=y*pi/180;

%%Haversine formulae 
%
%%interior points
dx(ii,:) = 2*R*asin(sqrt( sin((y(ii+1,:)-y(ii,:))/2).^2 + cos(y(ii,:)).*cos(y(ii+1,:)).*...
    sin((x(ii+1,:)-x(ii,:))/2).^2 )); %[m]
dy(:,jj) = 2*R*asin(sqrt( sin((y(:,jj+1)-y(:,jj))/2).^2 + cos(y(:,jj)).*cos(y(:,jj+1)).*...
    sin((x(:,jj+1)-x(:,jj))/2).^2 )); %[m]
%
%%boundaries
dx(1,:) = 2*R*asin(sqrt( sin((y(2,:)-y(1,:))/2).^2 + cos(y(1,:)).*cos(y(2,:)).*...
    sin((x(2,:)-x(1,:))/2).^2 )); %[m]
dx(M,:) = 2*R*asin(sqrt( sin((y(M,:)-y(M-1,:))/2).^2 + cos(y(M-1,:)).*cos(y(M,:)).*...
    sin((x(M,:)-x(M-1,:))/2).^2 )); %[m]
dy(:,1) = 2*R*asin(sqrt( sin((y(:,2)-y(:,1))/2).^2 + cos(y(:,1)).*cos(y(:,2)).*...
    sin((x(:,2)-x(:,1))/2).^2 )) ; %[m]  
dy(:,N) = 2*R*asin(sqrt( sin((y(:,N)-y(:,N-1))/2).^2 + cos(y(:,N-1)).*cos(y(:,N)).*...
    sin((x(:,N)-x(:,N-1))/2).^2 )) ; %[m]