clearvars
close all
clc

defcms
defiles

addpath(genpath('E:\MATLAB_library\'))

detpath('netcdf_files')

fig_folder = '/Users/andreacosta/Desktop/QUESTIONSacoc/Q2_uvel/barotropic/';
file_folder = '/Users/andreacosta/Desktop/QUESTIONSacoc/Q2_uvel/barotropic/';

season='clima';
load([file_folder,'psiA_',season],'psiA');


latu=ncread(uvelnc,'lat');
lonu=ncread(uvelnc,'lon');
lonM=repmat(lonu,1,length(latu));
latM=repmat(latu',length(lonu),1);

depthM=ncread('grid.nc','Depth');%depth matrix [m]
indsd = find(depthM==0);


%%%
%%% BAROTROPIC EQUATION - steady and linearized
%%%
% psit = psiA*1e-6*1000;

psit = psiA*1e-6;

%%%
% $A_{H}\nabla^2\psi + \frac{f}{2\alpha H}\nabla^2\psi - \frac{f}{H}J(H,\psi) - 
% \beta\frac{\partial\psi}{\partial x} = 
% -\frac{f}{\rho_0g H}\nabla\times\mathbf{\tau} - \frac{g}{2\alpha\rho_0
% H}\int_{-H}^{0}z\nabla^2\rho dz + \frac{g}{\rho_0 H}\int_{-H}^0z
% J(H,\rho) dz$
%
% $\alpha=\sqrt{\frac{f}{2A_V}}$


Ah = 10;
Av = 1e-5;

om = 2*pi/(24*60*60); % Earth rotation angle velocity
% Set Coriolis force coefficients
f = 2*om*sin(latM(1,:)*pi/180); %[1/s]
f = repmat(f,numel(lonM(:,1)),1);
[~,~, beta,~] = soseder(lonM,latM,f,'x',1);

alpha = sqrt(abs(f)/2/Av);


rho0 = 9.998e2; %http://mitgcm.org/public/r2_manual/latest/online_documents/node101.html

g=9.806-.5*(9.832-9.780)*cos(2*latu*pi/180); % gravity [m/s2]
g = repmat(g',numel(lonu),1);


%%% ONE %%%
[nap4x,~,nap4y,~] = soseder(lonM,latM,psit,'x',4);
[nap2x,~,nap2y,~] = soseder(lonM,latM,psit,'x',2);
nap4 = nap4x+nap4y+2*nap2x.*nap2y;
one = Ah*nap4;

figure;pcolor(lonu,latu,one'*1e12);shading flat;title('one')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)
set(gcf, 'Position', [0 0 1200 600])
print('-f1',[fig_folder,'one_clima'],'-dpng')

%%% TWO %%%
% [nap2x,~,nap2y,~] = soseder(lonM,latM,psit,'x',2);
nap2 = nap2x+nap2y;
two = f./2./alpha./(-depthM).*nap2;

figure;pcolor(lonu,latu,two'*1e14);shading flat;title('two')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)
set(gcf, 'Position', [0 0 1200 600])
print('-f2',[fig_folder,'two_clima'],'-dpng')

%%% THREE %%%
[dpdx,~,dpdy,~] = soseder(lonM,latM,psit,'x',1);
%JEBAR
[dHdx,~,dHdy,~] = soseder(lonM,latM,-depthM,'x',1);
JB = dpdx.*dHdy - dHdx.*dpdy;
three = -f./(-depthM).*JB;

figure;pcolor(lonu,latu,three'*1e12);shading flat;title('three')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)
set(gcf, 'Position', [0 0 1200 600])
print('-f3',[fig_folder,'three_clima'],'-dpng')

%%% FOUR %%%
four = -beta.*dpdx;

figure;pcolor(lonu,latu,four'*1e13);shading flat;title('four')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)
set(gcf, 'Position', [0 0 1200 600])
print('-f4',[fig_folder,'four_clima'],'-dpng')

%%% FIVE %%%
%load tau
myFile = matfile([fig_folder,'tauX_clima']);
% myFile = matfile([fig_folder,'tauX_ave_',season]);
taux = myFile.taux(:, :); %[kg/m3]
myFile = matfile([fig_folder,'tauY_clima']);
% myFile = matfile([fig_folder,'tauY_ave_',season]);
tauy = myFile.tauy(:, :); %[kg/m3]

[dtydx,~,~,~] = soseder(lonM,latM,tauy,'x',1); %clear tauy
[~,~,dtxdy,~] = soseder(lonM,latM,taux,'x',1); %clear taux

rottau = dtydx - dtxdy;

% close all
% figure;hold on
% m_proj('stereographic','lat',-90,'long',0,'radius',60);
% m_pcolor(lonM,latM,taux);
% shading flat
% m_grid('xtick',10,'tickdir','in','ytick',-[80 70 60 50 40],'linest','--');
% m_coast('patch',[.7 .7 .7],'edgecolor','none');
% m_coast('linewidth',0.5,'color','k');
% caxis([-.2 .2])
% colorbar
% cmocean('curl', 'NLevels', 20, 'pivot', 0)
% set(gcf, 'Position', [0 0 1000 1000])


% figure;hold on
% m_proj('miller','lon',[min(min(lonu)) max(max(lonu))],'lat',[min(min(latu)) -60]);
% m_contourf(lonM,latM,taux,20);
% m_grid('xtick',10,'tickdir','in','ytick',-[80 70 60 50 40],'linest','--');
% m_coast('patch',[.7 .7 .7],'edgecolor','none');
% m_coast('linewidth',0.5,'color','k');
% caxis([-.2 .2])
% colorbar
% cmocean('curl', 'NLevels', 20, 'pivot', 0)
% set(gcf, 'Position', [0 0 2000 1400])
% % print('-f2',[fig_folder,'rottau_clima_zoom'],'-dpng')
% 
% figure;hold on
% contourf(lonM(:,1:161),latM(:,1:161),taux(:,1:161),10);
% colorbar
% caxis([-.2 .2])
% cmocean('curl', 'NLevels', 20, 'pivot', 0)
% set(gca,'fontsize',16)
% set(gcf, 'Position', [0 0 3000 400])
% % print('-f3',[fig_folder,'taux_clima_zoom'],'-dpng')
% print('-f3',[fig_folder,'taux_',season,'_zoom'],'-dpng')




five = -f/rho0./g./(-depthM).*rottau;

figure;pcolor(lonu,latu,five'*1e18);shading flat;title('five')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)
set(gcf, 'Position', [0 0 1200 600])
print('-f5',[fig_folder,'fived_clima'],'-dpng')

%%% SIX %%%
load('depth_mask.mat');
DRF=ncread('grid.nc','DRF'); %dz

% for ii=1:size(depthM,1)
%     for jj=1:size(depthM,2)
%     
%         na2r = squeeze(ncread([fig_folder,'rho_nabla2_clima.nc'],'nar2',[ii jj 1],[1 1 Inf]));
%         
% %         dint = nansum(squeeze(depth_mask(ii,jj,:)).*-depth.*na2r.*DRF);
%         dint = nansum(-depth.*na2r.*DRF);
%         
%     end
% end
% 
% six = -g/2./alpha/rho0./(-depthM).*dint;
% six(indsd)=NaN;
% 
% figure;pcolor(lonu,latu,six'*1e12);shading flat;title('six')
% caxis([-5 5]);colorbar
% cmocean('delta', 'NLevels', 20, 'pivot', 0)
% set(gcf, 'Position', [0 0 1200 600])


%%% SEVEN %%%
%load int z J(H,rho)
% for ii=1:size(depthM,1)
%     for jj=1:size(depthM,2)
%     
%         jac = squeeze(ncread('jac_depth_rho_clima.nc','Jacobian of depth and density',[ii jj 1],[1 1 Inf]));
%         
% %         dintj = nansum(squeeze(depth_mask(ii,jj,:)).*-depth.*jac.*DRF);
%         dintj = nansum(-depth.*jac.*DRF);
%         
%     end
% end
% load([fig_folder,'zjint_clima.mat']);
% seven = g/rho0./depthM.*dintj;
seven = g/rho0./depthM.*zjint;
figure;imagesc(seven*1e8);title('seven')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)
set(gcf, 'Position', [0 0 1200 600])



figure;hold on
subplot(5,2,1)
pcolor(lonu,latu,one'*1e12);shading flat;title('Dissip.')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)

subplot(5,2,2)
pcolor(lonu,latu,two'*1e14);shading flat;title('Vorticity cons.')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)
subplot(5,2,3)
pcolor(lonu,latu,three'*1e12);shading flat;title(' ')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)

subplot(5,2,4)
pcolor(lonu,latu,four'*1e12);shading flat;title('BARBE')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)

subplot(5,2,5)
pcolor(lonu,latu,five'*1e12);shading flat;title('Wind forc.')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)

subplot(5,2,6)
pcolor(lonu,latu,six'*1e12);shading flat;title('Bottom Fric.')
caxis([-5 5]);colorbar

subplot(5,2,7)
pcolor(lonu,latu,seven'*1e12);shading flat;title('JEBAR')
caxis([-5 5]);colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)

subplot(5,2,9:10)
pcolor(lonM,latM,one+two+three+four-five-seven);shading flat;title('SUM')
caxis([-10 10]);
colorbar
cmocean('delta', 'NLevels', 20, 'pivot', 0)





figure;hold on
plot(lonu,nanmean(one,2),'r')
plot(lonu,nanmean(two,2),'b')
plot(lonu,nanmean(three,2),'k')
plot(lonu,nanmean(four,2),'g')
plot(lonu,nanmean(five,2),'c')
plot(lonu,nanmean(seven,2),'m')

figure;hold on
plot(lonu,nanmean(one,2),'r')
plot(lonu,nanmean(two,2),'b')
plot(lonu,nanmean(three,2),'k')

plot(lonu,nanmean(five,2),'c')
plot(lonu,nanmean(seven,2),'m')


keyboard

figure;hold on
set(gcf, 'Position', [0 0 3000 400])
m_proj('miller','lon',[min(min(lonu)) max(max(lonu))],'lat',[min(min(latu)) -60]);
m_contourf(lonM,latM,zjint,20);
m_grid('xtick',10,'tickdir','in','ytick',-[80 70 60 50 40],'linest','--');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_coast('linewidth',0.5,'color','k');
caxis([-3 3]*1e-2)
colorbar
cmocean('balance', 'NLevels', 20, 'pivot', 0)

figure;hold on
set(gcf, 'Position', [0 0 1000 1000])
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(lonM,latM,psit);
shading flat
m_grid('xtick',10,'tickdir','in','ytick',-[80 70 60 50 40],'linest','--');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_coast('linewidth',0.5,'color','k');
caxis([-80 80])
colorbar
cmocean('balance', 'NLevels', 20, 'pivot', 0)

figure;hold on
set(gcf, 'Position', [0 0 1000 1000])
m_proj('stereographic','lat',-90,'long',0,'radius',60);
[c,h]=m_contourf(lonM,latM,psit,20);
m_contour(lonM,latM,f./depthM,[10.^(-7.45) 10.^(-7.45)],'linewidth',2,'linecolor','y');
m_grid('xtick',10,'tickdir','in','ytick',-[80 70 60 50 40],'linest','--');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_coast('linewidth',0.5,'color','k');
caxis([-80 80])
colorbar
cmocean('balance', 'NLevels', 20, 'pivot', 0)

figure;hold on
m_proj('miller','lon',[min(min(lonu)) max(max(lonu))],'lat',[min(min(latu)) -60]);
m_contourf(lonM,latM,psit,20);
m_grid('xtick',10,'tickdir','in','ytick',-[80 70 60 50 40],'linest','--');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_coast('linewidth',0.5,'color','k');
caxis([-10 10])
colorbar
cmocean('balance', 'NLevels', 20, 'pivot', 0)
set(gcf, 'Position', [0 0 1800 400])

figure;hold on
set(gcf, 'Position', [0 0 1800 400])
m_proj('miller','lon',[min(min(lonu)) max(max(lonu))],'lat',[min(min(latu)) -60]);
m_contourf(lonM,latM,psit,20);
m_grid('xtick',10,'tickdir','in','ytick',-[80 70 60 50 40],'linest','--');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_coast('linewidth',0.5,'color','k');
caxis([-10 10])
colorbar
cmocean('balance', 'NLevels', 20, 'pivot', 0)


figure;hold on
contourf(lonM(:,1:161),latM(:,1:161),psit(:,1:161),[-8 -5 -3 -1 0 1 3 5 8]);
colorbar
caxis([-8 8])
colormap(cmocean('balance', 'NLevels', 20, 'pivot', 0))
set(gca,'fontsize',16)
title([season,' barotropic stream function'],'interpreter','latex','fontsize',20)
set(gcf, 'Position', [0 0 3000 400])

% print('-f2',[fig_folder,'psi_',season,'_zoom'],'-dpng')
print('-f2',[fig_folder,'psi_clima_zoom'],'-dpng')
%%%various terms magnitude
%%%average


%%%median


%%%same but only in acoc zone




