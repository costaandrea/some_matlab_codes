clearvars
close all
clc


detpath('netcdf_files')

%%%flienames
defiles

%%%grid stuff
lat=ncread(vvelnc,'lat');
lon=ncread(vvelnc,'lon');
depth=ncread(vvelnc,'depth');
time=ncread(vvelnc,'time');%01JAN2008 to 31DEC2012

lonM=repmat(lon,1,numel(lat));
latM=repmat(lat',numel(lon),1);

y1=findnearest(-90,lat);
y2=findnearest(-60,lat);

latx = findnearest(-71.7,lat);
lo = findnearest(170.5,lat);
lo2 = max(findnearest(258,lon));

HFacW =ncread('grid.nc','hFacW'); %free dz pcntg
DRF=ncread('grid.nc','DRF'); %dz
DYG=ncread('grid.nc','DYG'); %dy [m]
DXG=ncread('grid.nc','DXG'); %dy [m]
depthM = ncread('grid.nc','Depth');%depth matrix
depthM = depthM(lo:lo2,latx);
indsd = find(depthM==0);   %land indeces periodic matrix


g = 9.806-.5*(9.832-9.780)*cos(2*lat*pi/180); % gravity [m/s2]

heat=nan(lo2-lo+1,numel(DRF),numel(time));
intZ=nan(numel(lon(lo:lo2)),numel(time));
HT=nan(1,numel(time));
for ttt=254:numel(time)
%     fprintf(['\n','TIME ',num2str(ttt),'/',num2str(numel(time)),'\n'])
    
parfor la=1:numel(lat)
%     fprintf([num2str(la),'/',num2str(numel(lat)),'...\n'])
    
    detpath('netcdf_files')
%     time=ncread(vvelnc,'time');
    %fprintf('Loading Vvel... \n')
    Uvel = ncread(uvelnc,'Uvel',[lo la 1 ttt],[lo2-lo+1 1 Inf 1]);
    Uvel(Uvel==-9999)=NaN;

    
    %fprintf('Loading density... \n')
    rho = ncread(densnc,'GAMMA',[lo la 1 ttt],[lo2-lo+1 1 Inf 1]); 
    rho(rho==-9999)=NaN;
    rho=rho+1000;
    %%calculate pression
    press = gsw_p_from_z(-depth,lat(la));
    %rho(450,1)*9.81*depth(1) *1e-5 *10 %[dbar]
    
    
    %fprintf('Loading salinity... \n')
    Salt = ncread(saltnc,'Salt',[lo la 1 ttt],[lo2-lo+1 1 Inf 1]); 
    Salt=squeeze(Salt);
    
    %absolute salinity
    SA = gsw_SA_from_SP(Salt,press,lon(lo:lo2),lat(la));
%     clear Salt
    
    
%     fprintf('Loading temperature... \n')
    Temp = ncread(thetanc,'Theta',[lo la 1 ttt],[lo2-lo+1 1 Inf 1]);
    Temp(Temp==-9999)=NaN;
    Temp=squeeze(Temp);
    

    %fprintf('Calculating Cp... \n')
    Cp = gsw_cp_t_exact(SA,Temp,press); %[ J/(kg*K) ]
%     clear press SA
    
    
    %fprintf('Calculating heat flux... \n')
    heat = Cp.*squeeze(rho).*squeeze(Uvel).* (Temp+273.15).* squeeze(HFacW(lo:lo2,la,:))...
        .*repmat(DRF',numel(lon(lo:lo2)),1).*repmat(DXG(lo:lo2,la),1,size(Temp,3)); %[W]
    
    
    intZ_cums = nancumsum(heat,2);
    
%     intZ(:,ttt) = nansum(heat(:,:,ttt) ,2); 
%     clear Vvel rho Cp
    
%     HT(ttt) = nansum(intZ(:,ttt),1); 



    cd  /Volumes/AC_Thunder_2/SOSE_HT/
    ncfile = ['SOSE_HT_u_',num2str(ttt,'%03d'),'_',num2str(la,'%03d'),'.nc']; 

 if exist(ncfile,'file')>0
     delete(ncfile)
 end
 
 %%dimension
 nccreate(ncfile,'lon','Dimensions',{'lon',numel(lon)},'DeflateLevel',5);
 ncwrite(ncfile,'lon',lon,1); 
 nccreate(ncfile,'lat','Dimensions',{'lat',1},'DeflateLevel',5);
 ncwrite(ncfile,'lat',lat(la),1); 
 nccreate(ncfile,'z','Dimensions',{'z',numel(depth)},'DeflateLevel',5);
 ncwrite(ncfile,'z',depth,1); 
 nccreate(ncfile,'time','Dimensions',{'time',1},'DeflateLevel',5);
 ncwrite(ncfile,'time',time(ttt),1);
 
 %%variables
 nccreate(ncfile,'heat','Dimensions',{'lon',numel(lon), 'z',numel(depth), 'lat',1, 'time',1},'DeflateLevel',5); %create nc file and variable 
 ncwriteatt(ncfile,'heat','Units','[W]');
 ncwriteatt(ncfile,'heat','Name','Zonal Heat Transport section');
 ncwriteatt(ncfile,'heat','Dimensions','[lon,z,lat,t]');
 ncwrite(ncfile,'heat',heat,[1, 1, 1, 1]); 
%
 nccreate(ncfile,'intZ_cums','Dimensions',{'lon',numel(lon), 'z',numel(depth), 'lat',1, 'time',1},'DeflateLevel',5); %create nc file and variable 
 ncwriteatt(ncfile,'intZ_cums','Units','[W]');
 ncwriteatt(ncfile,'intZ_cums','Name','Zonal Heat Transport Cumulative Vertical Integral');
 ncwriteatt(ncfile,'intZ_cums','Dimensions','[lon,z,lat,t]');
 ncwrite(ncfile,'intZ_cums',intZ_cums,[1, 1, 1, 1]); 
%
 
 
end%lat

end

 