close all
clearvars
clc

detpath('netcdf_files')

%%%filenames
defiles

%%%grid stuff
lon = ncread(stratnc,'lon');
lat = ncread(stratnc,'lat'); 
time=ncread(stratnc,'time');%01JAN2008 to 31DEC2012

lonM = repmat(lon,1,numel(lat));
latM = repmat(lat',numel(lon),1);

y2=findnearest(-60,lat);
lo=294;
lo2=775;

SSH = ncread(sshnc,'ETAN');
tmp= nanmean(SSH,3);
s_prime = SSH - tmp;

depthMall=ncread('grid.nc','Depth');


%%%make movie
f=figure;
set(gcf, 'Position', [0 0 1500 600])

subplot(2,2,[1 2]);hold on
title(datestr(time(2),'dd-mmm-yy'))
tmp=s_prime(:,:,2);
tmp(depthM(1:size(tmp,1),1:size(tmp,2))==0)=NaN;
pcolor(lonM(:,1:y2+40),latM(:,1:y2+40),tmp(:,1:y2+40));shading flat
xlabel('Lon')
ylabel('Lat')
caxis([-.2,.2])
colormap(bone(20))
cb=colorbar;
ylabel(cb,'$SSH$','interpreter','latex')
axis tight
contour(lonM(:,1:y2+40),latM(:,1:y2+40),depthM(:,1:y2+40),[1800 1800],'w','linewidth',1.5)

subplot(2,2,[3 4]);hold on
plot(time(2:end),HT)
p2=plot(time(2),HT(2),'.','color','blue','markersize',14);
datetick('x','mm-yy','keepticks')
axis tight
xlabel('Time')
ylabel('HT [W]')

set(gcf, 'color', 'white');
F = getframe(gcf);
im = frame2im(F);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,'HT_SSH_whole.gif','gif', 'Loopcount',1);

pause(.6)
%waitforbuttonpress
set(p2,'Visible','off')

for tt=3:numel(time)
    
    subplot(2,2,[1 2]);hold on
    title(datestr(time(tt),'dd-mmm-yy'))
    tmp=s_prime(:,:,tt);
    tmp(depthM(1:size(tmp,1),1:size(tmp,2))==0)=NaN;
    pcolor(lonM(:,1:y2+40),latM(:,1:y2+40),tmp(:,1:y2+40));shading flat
    xlabel('Lon')
    ylabel('Lat')
    caxis([-.2,.2])
    colormap(bone(20))
    cb=colorbar;
    ylabel(cb,'$SSH$','interpreter','latex')
    axis tight
    contour(lonM(:,1:y2+40),latM(:,1:y2+40),depthM(:,1:y2+40),[1800 1800],'w','linewidth',1.5)

    subplot(2,2,[3 4]);hold on
    p2=plot(time(tt),HT(tt),'.','color','blue','markersize',14);
    
    drawnow
    
    F = getframe(gcf);
    im = frame2im(F);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,'HT_SSH_whole.gif','gif','WriteMode','append');
    
    pause(.6)
    
    %if tt==29||tt==62||tt==143||tt==147||tt==195
    %waitforbuttonpress
    %end

    set(p2,'Visible','off')
    
end