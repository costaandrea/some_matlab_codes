close all
clearvars
clc

detpath('netcdf_files')
defiles

latu=ncread(uvelnc,'lat');
lonu=ncread(uvelnc,'lon');
latv=ncread(vvelnc,'lat');
lonv=ncread(vvelnc,'lon');

lonM=repmat(lonu,1,numel(latu));
latM=repmat(latu',numel(lonu),1);


U = ncread('SOSE_velocity_surface_highp_3mo.nc','Uvel_highp');
V = ncread('SOSE_velocity_surface_highp_3mo.nc','Vvel_highp');
time = ncread('SOSE_velocity_surface_highp_3mo.nc','time');


[Xv,Yv,tv] = ndgrid(lonv,latv,time,'x'); clear lonv latv
[Xu,Yu,tu] = ndgrid(lonu,latu,time,'x');
V = interpn(Xv,Yv,tv, V, Xu,Yu,tu, 'linear');
clear Xu Yu Xv Yv tu tv 


%%calc vorticity
fprintf('Calculating vorticity... \n')
omega = NaN(size(U));
parfor tt=1:numel(time)
    [~, ~, dudy, ~] = soseder(lonM,latM,U(:,:,tt),'x',1); 
    [dvdx, ~, ~, ~] = soseder(lonM,latM,V(:,:,tt),'x',1); 
    omega(:,:,tt) = dvdx - dudy;
end

 
%%calc EKE
fprintf('Calculating  EKE... \n')
EKE = .5*(U.^2+V.^2);
clear U V


%%%seasonal averages
seasons{1}='Spring';
seasons{2}='Summer';
seasons{3}='Autumn';
seasons{4}='Winter';

omega_s = nan(size(EKE,1),size(EKE,2),4);
EKE_s = nan(size(EKE,1),size(EKE,2),4);
parfor  sss = 1:4
    season = seasons{sss};
    
    inds=indseason(season,time); %seasonal indexes

    omega_s(:,:,sss) = nanmean(omega(:,:,inds),3); 
    EKE_s(:,:,sss) = nanmean(EKE(:,:,inds),3); 
end
clear omega
clear EKE

% %%enstrophy
% enstr_s = omega_s.^2;

figure; imagesc(omega_s(:,:,2)./omega_s(:,:,1))
caxis([.5 2])

figure; imagesc(EKE_s(:,:,2)./EKE_s(:,:,1))
caxis([.5 2])


%%%stratification
fprintf('Stratification... \n')
N2_int_s = nan(size(omega_s,1),size(omega_s,2),4);
g=9.8;
rho0=1028;
parfor sss=1:4
    
season=seasons{sss};

drdz = season_ave(stratnc,'Strat', season,time, [1 1 1],[Inf Inf Inf]);
drdz(drdz==-9999)=NaN;
tmp = -drdz * g / rho0; clear drdz
N2_int_s(:,:,sss) = nansum(tmp,3); clear tmp
end


