function var_seas = season_ave(filename,varname,season,time, otherinds_in, otherinds_fin)
%season_ave calculates the seasonal avergae of a variable.
%
%   Inputs: FILENAME: the name of the netcdf file containing the variable
%           VARNAME: the name of the variable
%           SEASON:  the name of the season 'spring' 'summer' 'autumn'
%                                           'winter'
%           TIME: time vector
%           OTHERINDS_IN: starting indexes of coordinates other then time
%           OTHERINDS_FIN: final indexes of coordinates other then time
%
% This script needs indseason.m and conschunks.m


%%%select season
inds=indseason(season,time); %seasonal indexes
cinds = conschunks(inds); %consecutive indexes


if numel(otherinds_in)==3 %4D fields
%average
for kkk=1:length(cinds)
    if kkk == 1
        var =  ncread(filename,varname,[otherinds_in cinds{kkk}(1)],...
                          [otherinds_fin cinds{kkk}(end)-cinds{kkk}(1)+1]);
        var(var==-9999)=NaN;
        var = squeeze(nanmean(var,4)); %partial seasonal average
        
    else
        var_t = ncread(filename,varname,[otherinds_in cinds{kkk}(1)],...
                          [otherinds_fin cinds{kkk}(end)-cinds{kkk}(1)+1]);
        var_t(var_t==-9999)=NaN;
        var_t = squeeze(nanmean(var_t,4)); %partial seasonal average
        
        var = var+var_t;
        clear var_t
        
    end
end
var_seas = var/length(cinds); %final seasonal average


elseif numel(otherinds_in)==2 %3D fields
%average
for kkk=1:length(cinds)
    if kkk == 1
        var =  ncread(filename,varname,[otherinds_in cinds{kkk}(1)],...
                          [otherinds_fin cinds{kkk}(end)-cinds{kkk}(1)+1]);
        var(var==-9999)=NaN;
        var = squeeze(nanmean(var,3)); %partial seasonal average
        
    else
        var_t = ncread(filename,varname,[otherinds_in cinds{kkk}(1)],...
                          [otherinds_fin cinds{kkk}(end)-cinds{kkk}(1)+1]);
        var_t(var_t==-9999)=NaN;
        var_t = squeeze(nanmean(var_t,3)); %partial seasonal average
        
        var = var+var_t;
        clear var_t
        
    end
end
var_seas = var/length(cinds); %final seasonal average

end