function hcb = AAplot(lon,lat,field,label,varargin)
%AAplot make a standard plot for a 2D field around Antarctica.
%   [hf,hcb] = AAplot(LON,LAT,FIELD,LABEL)
%
%   LON and LAT are 2D matrices.
%   LABEL is the colorbar label. For default the colorbar 
%        label interpreter is set to 'latex' and the fontsize to 16.
%   HF is the figure handle. HCB is the colorbar handle.
%
%   The default projection is 'lambert'. 
%   Shading is 'flat'.
%   I use m_coast and m_grid (no m_gssh).
%
%   N.B. I will mess with the FigureRenderer


%to make m_pcolor not screw up when put together with m_grid
set(0,'DefaultFigureRenderer','zbuffer') %'painters' screws up the colorbar
                                         %probably because of
                                         %diverging_map.m used to define
                                         %the colormap

if nargin>4
    projection = varargin{1};
else
    projection = 'lambert';
end
                                         
                                         
% hf = figure; hold on
hold on

set(gcf, 'Position', [0 0 1000 1000])

m_proj(projection,'lon',[min(min(lon)) max(max(lon))],'lat',[min(min(lat)) max(max(lat))]);


m_pcolor(lon,lat,field); 
shading flat

m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_coast('linewidth',0.5,'color','k');

m_grid('linewi',2,' tickdir','out');
% m_grid('box','fancy','tickdir','in');


hcb=colorbar;
set(get(hcb,'ylabel'),'String',label,'interpreter','latex', 'fontsize',24)
set(gca, 'fontsize',24)


%set(0,'DefaultFigureRenderer','auto')