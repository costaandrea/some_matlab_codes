function inds = indseason(season,time) 
    
inds=[];
    
 if strcmp(season,'Spring')
%%%spring
    for yy=8:12

      start=findnearest(datenum(['22SEP20',num2str(yy, '%02i')]),time);
      fin=findnearest(datenum(['21DEC20',num2str(yy, '%02i')]),time);
      inds=[inds start:fin];
    end

 elseif strcmp(season,'Summer')
%%%summer
    for yy=8:12

      start=findnearest(datenum(['22DEC20',num2str(yy-1, '%02i')]),time);
      fin=findnearest(datenum(['21MAR20',num2str(yy, '%02i')]),time);
      inds=[inds start:fin];
    end
    
 elseif strcmp(season,'Autumn')
%%autumn
    for yy=8:12

        start=findnearest(datenum(['22MAR20',num2str(yy, '%02i')]),time);
        fin=findnearest(datenum(['21JUN20',num2str(yy, '%02i')]),time);
        inds=[inds start:fin];
    end
    
 elseif strcmp(season,'Winter')
%%%winter
    for yy=8:12

        start=findnearest(datenum(['22JUN20',num2str(yy, '%02i')]),time);
        fin=findnearest(datenum(['21SEP20',num2str(yy, '%02i')]),time);
        inds=[inds start:fin];
    end
    
 else
     
     error('No corresponding season. \n Check cases.')

 end