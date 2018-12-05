clear vars
close all
clc

addpath('/Users/andreacosta/Desktop/SOSE')
%to make m_pcolor not screw up when put together with m_grid
set(0,'DefaultFigureRenderer','zbuffer')

projection = 'lambert';

detpath('netcdf_files')

%%%define colormap
ss=diverging_map(0:.001:1,[0.370000 0.620000 0.630000],[1.000000 0.130000 0.320000]);
colormap(ss)
cm=colormap;
cm(1,:)= [1 1 1];


%%%filenames
defiles

%%%grid stuff
lat=ncread(uvelnc,'lat');
lon=ncread(uvelnc,'lon');
lat2=ncread(vvelnc,'lat');
lon2=ncread(vvelnc,'lon');

%%%select domain
time=ncread(thetanc,'time');
T1=min(findnearest(time,datenum('01JAN2009')));
T2=min(findnearest(time,datenum('01FEB2009')));
time=time(T1:T2);

x1=findnearest(40,lon);
x2=findnearest(70,lon);
y1=findnearest(-72,lat);
y2=findnearest(-50,lat);

lonMAP=lon2(x1:x2); %useful for interpolating and mapping
latMAP=lat2(y1:y2);


%%%LOAD info on partial grid cells 
DYG   =ncread('grid.nc','DYG',[x1 y1],[x2-x1+1 y2-y1+1]);
DRF   =ncread('grid.nc','DRF');
hFaCW =ncread('grid.nc','hFacW',[x1 y1 1],[x2-x1+1 y2-y1+1 Inf]);
DXG   =ncread('grid.nc','DXG',[x1 y1],[x2-x1+1 y2-y1+1]);
hFacS =ncread('grid.nc','hFacS',[x1 y1 1],[x2-x1+1 y2-y1+1 Inf]);


%%%INTERPOLATE U at V points (where also T,S are)

%matrix of lat lon for U
depthU=ncread(uvelnc,'depth');

%matrix of lat lon for V
depthV=ncread(vvelnc,'depth');

%read U and V
Uvel=squeeze(ncread(uvelnc,'Uvel',[x1 y1 1 T1],[x2-x1+1 y2-y1+1 Inf T2-T1+1]));
Vvel=squeeze(ncread(vvelnc,'Vvel',[x1 y1 1 T1],[x2-x1+1 y2-y1+1 Inf T2-T1+1]));
Uvel(Uvel==-9999)=NaN;
Vvel(Vvel==-9999)=NaN;

%interpolation of U on V's grid-points
[Xu,Yu,Zu,Tu]=ndgrid(lon(x1:x2),lat(y1:y2),depthU,time);
[Xv,Yv,Zv,Tv]=ndgrid(lon2(x1:x2),lat2(y1:y2),depthV,time);

UvelINT=interpn(Xu,Yu,Zu,Tu,Uvel,Xv,Yv,Zv,Tv,'linear',NaN);


%%%CALCULATE TRANSPORTS
DYGr=repmat(DYG,1,1,numel(depthU),numel(time));
DRFr_temp=repmat(DRF,1,size(Uvel,1),size(Uvel,2),numel(time));
DRFr=permute(DRFr_temp,[2 3 1 4]);
hFaCWr=repmat(hFaCW,1,1,1,numel(time));

Utr = UvelINT.*DYGr.*DRFr.*hFaCWr;%%%


DXGr=repmat(DXG,1,1,numel(depthU),numel(time));
hFacSr=repmat(hFacS,1,1,1,numel(time));

Vtr = Vvel.*DXGr.*DRFr.*hFacSr;%%%


%column-integrated transports
intUtr = squeeze(nansum(Utr,3));
intVtr = squeeze(nansum(Vtr,3));

%time average of transports
TintUtr = squeeze(nansum(intUtr,3));
TintVtr = squeeze(nansum(intVtr,3));



%%%PLOTTING time series of transport through a section
% meridional section %
xs=min(findnearest(50,lonMAP)); %make sure that this is in the domain!!

figure
pcolor(time,latMAP,squeeze(intUtr(xs,:,:)))
shading flat;
colormap(cm);
xlabel('Time')
ylabel('Latitude')
colorbar

% zonal section %
ys=min(findnearest(50,latMAP)); %make sure that this is in the domain!!

figure
pcolor(lonMAP,time,squeeze(intUtr(:,ys,:))')
shading flat;
colormap(cm);
xlabel('Longitude')
ylabel('Time')
colorbar

% arbitrary section %
xs1=41;  %make sure that these are in the domain!!
xs2=69;
ys1=-79;
ys2=-62;
N=20; %how many points on the section

xq=linspace(xs1,xs2,N); %from origin to end of the section
yq=linspace(ys1,ys2,N);

%interpolate over the section
for ttt=1:numel(time)
    Utr_sec(:,ttt) = interpn(lonMAP,latMAP,intUtr(:,:,ttt),xq,yq,'linear',NaN);
end

dq=sqrt((xq-xs1).^2 + (yq-ys1).^2); %axis along section

figure
pcolor(dq,time,Utr_sec')
shading flat;
colormap(cm);
colorbar
xlabel('Along-section axis')
ylabel('Time')



%%%PLOTTING time-av col-int transport

figure;hold on
set(gcf, 'Position', [0 0 1000 1000])
colormap(cm);

m_proj(projection,'lon',[min(lonMAP) max(lonMAP)],'lat',[min(latMAP) max(latMAP)]);

m_pcolor(lonMAP,latMAP,TintUtr'); 
shading flat;

m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_coast('linewidth',0.5,'color','k');

m_grid('box','fancy','tickdir','in');

h=colorbar;
caxis([-4e6 4e6])

set(get(h,'ylabel'),'String','Time av. of column int. U transport', 'fontsize',14);
xlabel('Latitude', 'fontsize',14,'fontweight','bold')
ylabel('Longitude', 'fontsize',14,'fontweight','bold')
set(gca, 'fontsize',14)



%%%PLOTTING stacked horiz slices of transport

depth=ncread('grid.nc','Depth',[x1 y1],[x2-x1+1 y2-y1+1]);
X=repmat(lonMAP,1,numel(latMAP));
Y=repmat(latMAP',numel(lonMAP),1);

figure
set(gcf, 'Position', [0 0 700 900])

%transport slices
%1
Z=repmat(depthV(1),numel(lonMAP),numel(latMAP));
TR_t=squeeze(Utr(:,:,1,1));

surf(X,Y,Z,TR_t);hold on
shading flat
colormap(cm)
cb=colorbar;
ylabel(cb,'m^3/s','interpreter','latex')
caxis([-100000 10000])

%freezeColors

%2
Z=repmat(depthV(31),numel(lonMAP),numel(latMAP));
TR_t=squeeze(Utr(:,:,31,1));

surf(X,Y,Z,TR_t,'EdgeColor','none');hold on
%colormap(cm) %change it if you want. But then de-comment freezeColors here above

%freezeColors

%3
Z=repmat(depthV(41),numel(lonMAP),numel(latMAP));
TR_t=squeeze(Utr(:,:,41,1));

surf(X,Y,Z,TR_t,'EdgeColor','none');hold on
%colormap(cm)

%freezeColors

%4
Z=repmat(depthV(45),numel(lonMAP),numel(latMAP));
TR_t=squeeze(Utr(:,:,45,1));

surf(X,Y,Z,TR_t,'EdgeColor','none');hold on
%colormap(cm)

freezeColors

%bathymetry
depth_modif = depth; %cheating just a little bit 'cuz surf and patch do not work well together
depth_modif(1,:)=6000;
depth_modif(end,:)=6000;
depth_modif(:,1)=6000;
depth_modif(:,end)=6000;

p=surf(lonMAP,latMAP,depth_modif','FaceColor','interp','EdgeColor','none');hold on
set(gca,'dataaspectratio',[1 cos(nanmean(nanmean(latMAP'))*pi/180) 60]);

view(165,19)
set(gca,'ZDir','Reverse')
camlight

set(p,'edgecolor','none') %gotta do it now or shading fucks it up
set(p,'FaceColor',[.8 .8 .8]);


axis([30 80 -80 -60 0 6000])
xlabel('Longitude')
ylabel('Latitude')
zlabel('Depth [m]')
