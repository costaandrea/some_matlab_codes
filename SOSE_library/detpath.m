function detpath(folder)

match=0;

if strcmp(getenv('systemroot'),'C:\WINDOWS')
    flag=1;
    match=1;
else
    flag=0;
    match=1;
end

if strcmp(folder,'netcdf_files')&&flag==1
    cd E:\SOSE\netcdf_files\
    match=1;  
elseif strcmp(folder,'netcdf_files')&&flag==0
    cd /Volumes/AC_Thunder_2/netcdf_files/
    match=1;
end
    
if strcmp(folder,'FIGURES')&&flag==1
    cd E:\SOSE\FIGURES
    match=1;    
elseif strcmp(folder,'FIGURES')&&flag==0
    cd /Users/andreacosta/Desktop/SOSE/FIGURES
    match=1;
end

if strcmp(folder,'SOSE')&&flag==1
    cd E:\SOSE\
    match=1;        
elseif strcmp(folder,'SOSE')&&flag==0
    cd /Users/andreacosta/Desktop/SOSE/
    match=1;
end

if strcmp(folder,'ERA')&&flag==1
    cd E:\ERAinterim\
    match=1;  
elseif strcmp(folder,'ERA')&&flag==0
    cd /Users/andreacosta/Desktop/ERAinterim/
    match=1;
end


if match==0
    error('Unrecognized folder. Check the name.')
end