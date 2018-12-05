function testsolution(psi,Uvel,Vvel, lonu,latu,lonM,latM,inds_l)

%%%%
%%% TEST SOLUTION
%%%%
psi = aperiodic(psi,1,2);
lonM = aperiodic(lonM,1,2);
latM = aperiodic(latM,1,2);


[vv,~,uv,~] = soseder(lonM,latM,psi,'x',1);
vv(inds_l)=NaN; %on land
uv(inds_l)=NaN;

diffu = Uvel-(-uv); %discrepancy between input and reconstructed zonal velocity
diffv = Vvel-(vv); %discrepancy between input and reconstructed zonal velocity


figure;hold on
subplot(1,2,1)
imagesc(lonu,latu, diffu);colorbar;%caxis([-300 700])
subplot(1,2,2)
imagesc(lonu,latu,diffv);colorbar;%caxis([-300 700])
colormap(cool(20))
set(gcf, 'Position', [0 0 800 1200])
%%the plots should show values near to zaro everywhere