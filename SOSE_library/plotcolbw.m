function plotcolbw(xx,yy,fun,inds,liml,limu,xlab,ylab,alphaVal)

defcms

Z1 = fun;

Z2 = fun;
Z2(~inds) = NaN;

f = figure();
set(gcf,'Position',[0 0 3500 500])
title('');

ax1 = axes();
p1=pcolor(xx,yy,Z1);
shading flat;
set(p1,'facealpha',alphaVal)
xlabel(ax1,xlab)
ylabel(ax1,ylab)
set(gca,'fontsize',18)

view(2);
ax2 = axes();
pcolor(xx,yy,Z2);
shading flat;
ax2.Visible = 'off';
set(gca,'fontsize',18)

colormap(ax1,bone(10));
colormap(ax2,VELcm);

% set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'position',[.08 .11 .03 .815]);
% set(get(cb1,'ylabel'),'String','whatever','interpreter','latex', 'fontsize',20);
cb2 = colorbar(ax2,'position',[.92 .11 .03 .815]);
% set(get(cb2,'ylabel'),'String','whatever','interpreter','latex', 'fontsize',20);
 
caxis(ax1,[liml limu]);
caxis(ax2,[liml limu]);

set(ax1,'Tag','keep');
set(ax2,'Tag','keep');

delete(findall(f,'Type','Axes','-not','Tag','keep'));

% Get the color data of the object that correponds to the colorbar
cdata = cb1.Face.Texture.CData;

% % Change the 4th channel (alpha channel) to (1-alphaVal)% of it's initial value (255)
% cdata(end,:) = uint8((1-alphaVal) * cdata(end,:));

% Ensure that the display respects the alpha channel
cb1.Face.Texture.ColorType = 'truecoloralpha';

% Update the color data with the new transparency information
cb1.Face.Texture.CData = cdata;