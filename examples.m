close all
clearvars
clc

%%PLOT 1 :: different shadings in same figure
figure
surf(0:10,0:10,repmat(1,11,11),peaks(11),'EdgeColor','none');
hold on;     %You don't need to use hold on again and again
surf(0:10,0:10,repmat(3,11,11),peaks(11),'FaceColor','flat','EdgeColor','none');
surf(0:10,0:10,repmat(5,11,11),peaks(11),'FaceColor','interp','EdgeColor','none');
view(-15,32);



%%PLOT 2 :: change cm in zone of interest
[x, y, z] = peaks(1000);

figure;
p = surf(x,y,z);
shading flat

cm=colormap;

Xm=-7;
XM=8;
caxis([Xm XM])

x1=2;
x2=4;
l1=findnearest(x1,linspace(Xm,XM,length(cm)));
l2=findnearest(x2,linspace(Xm,XM,length(cm)));

%blurring
n = floor((l2-l1+1)/2);
fade = linspace(0, 1, n)';
cm(l1+n,:) = [1 0 0];
cm(l1:l1+n-1,:)=fade * [1 0 0] + (1-fade) * cm(l1-1, :);
cm(l2:-1:l2-n+1,:)=fade * [1 0 0] + (1-fade) * cm(l2+1, :);
colormap(cm)

Xm=-7;
XM=8;
caxis([Xm XM])
colorbar



%%PLOT3 :: two colormaps on different zones of scalar field
[X1,Y1,Z1] = peaks(250);

[X2,Y2,Z2] = peaks(250);
idx = find(X2 > 0);
Z2(idx) = NaN;

f = figure();
set(gcf,'Position',[0 0 800 800])
title('');

ax1 = axes();
pcolor(X1,Y1,Z1);
shading flat;
xlabel(ax1,'stuff')
ylabel(ax1,'other stuff')

view(2);
ax2 = axes();
pcolor(X2,Y2,Z2);
shading flat;
ax2.Visible = 'off';

colormap(ax1,bone());
colormap(ax2,jet(26));

set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'position',[.88 .11 .03 .815]);
set(get(cb1,'ylabel'),'String','whatever','interpreter','latex', 'fontsize',20);
cb2 = colorbar(ax2,'position',[.08 .11 .03 .815]);
set(get(cb2,'ylabel'),'String','whatever','interpreter','latex', 'fontsize',20);

caxis(ax1,[-5 5]);
caxis(ax2,[-5 5]);

set(ax1,'Tag','keep');
set(ax2,'Tag','keep');

delete(findall(f,'Type','Axes','-not','Tag','keep'));


%%PLOT 4 :: arbitrary transect on scalar field
figure
z=peaks(50);

% Create x,y coordinates of the data
[x,y]=meshgrid(1:50);

% Plot Data and the slicing curve
surf(z);
hold on

X=[1 21 35 47 29 25 8];
Y=[5 19 24 26 14 39 47];

plot3(X,Y,-10*ones(1,numel(X)),'r','linewidth',3);
plot3(X,Y,10*ones(1,numel(X)),'r','linewidth',3);

patch([X fliplr(X)],[Y fliplr(Y)],[-10*ones(1,numel(X)) 10*ones(1,numel(X))],...
    'r','FaceAlpha',0.21)

axis([0 50 0 50])

patch([X(1:end-1);X(2:end);X(2:end);X(1:end-1)],[Y(1:end-1);Y(2:end);Y(2:end);Y(1:end-1)],...
    [-10*ones(1,numel(X)-1);-10*ones(1,numel(X)-1);10*ones(1,numel(X)-1);10*ones(1,numel(X)-1)],...
    'r','FaceAlpha',0.21)
