function ffolder = defpath(folder)

if strcmp(getenv('systemroot'),'C:\WINDOWS')
    flag=1;
elseif
    flag=0;
end

if strcmp(folder,'netcdf_files')&&flag==1
    ffolder = 'E:\SOSE\netcdf_files\';
    
elseif strcmp(folder,'netcdf_files')&&flag==0
    ffolder = '/Users/andreacosta/Desktop/SOSE/netcdf_files/';
end
    
if strcmp(folder,'FIGURES')&&flag==1
    ffolder = 'E:\SOSE\FIGURES';
        
elseif strcmp(folder,'netcdf_files')&&flag==0
    ffolder = '/Users/andreacosta/Desktop/SOSE/FIGURES';
end

if strcmp(folder,'SOSE')&&flag==1
    ffolder = 'E:\SOSE\';
        
elseif strcmp(folder,'netcdf_files')&&flag==0
    ffolder = '/Users/andreacosta/Desktop/SOSE/';
end
