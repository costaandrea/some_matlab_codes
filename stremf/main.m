clearvars
close all
clc

addpath('./functions')

%% load grid 
load('./data/grid_stuff.mat','latu');
load('./data/grid_stuff.mat','lonu');

lonM=repmat(lonu,1,length(latu));  %longitude matrix
latM=repmat(latu',length(lonu),1); %latitude matrix


load('./data/grid_stuff.mat','depthM');%bottom depth matrix [m]
inds_l = find(depthM==0);         %land indeces
not_inds_l = find(depthM>0);       %ocean indeces
%%% make it periodic: append the 2 first columns at the end  & 
                    % prepend the last 2 columns at the beginning
                    % USE periodic.m for this
               
depthM_pd = periodic(depthM,1,2); %make periodic depthM in the first direction
                                %by appending/prepending 2 columns
inds_l_pd = find(depthM_pd==0);   %land indeces periodic matrix
not_inds_l_pd = find(depthM_pd>0); %ocean indeces periodic matrix
depthM_pd(inds_l_pd)=NaN;         %NaN on land


%%% find the islands
a = depthM_pd; %to feed to the island finder
a(inds_l_pd)= 1; %on land
a(not_inds_l_pd)= 0; %on water

[islands, neigh, all_isl, allc] = islands_finder(a,4); %see  >>help islands_finder
                %Because of how it is built, it starts from Antarctica (index ==1)
inds_c_pd = find(allc); %indeces of the grid points nearby all the coasts
                            %%where circulations must be calculated

[Nx,Ny] = gradient(all_isl); %test this w/ a = [1 1 1 1 1 1 1; 1 0 0 0 0 0 1;  1 1 0 0 0 1 1; 1 1 0 0 0 1 1; 1 0 0 0 0 0 1; 1 0 0 0 0 0 1; 1 1 1 1 1 1 1]
Nx = sign(Nx); %x-component of normal vector 
Ny = sign(-Ny); %y-component of normal vector 

%% CALCULATION OF PSI

for ttt=1:1%2

    fprintf(['Step ',num2str(ttt),'\n\n'])
    
    %%load data
    myFile = matfile('./data/Uvel.mat');
    Uvel = myFile.Uvel(:,:,1, ttt); %zonal velocity (1st dimension) [m2/s]
    myFile = matfile('./data/Vvel.mat');
    Vvel = myFile.Vvel(:,:,1, ttt);  %meridional velocity (2nd dimension) [m2/s]
    
    Uvel(inds_l)=0; %zeros on land (before it was NaN)
    Vvel(inds_l)=0; %zeros on land
    
    %%load boundary condition at open north points
    myFile = matfile('./data/bc_north.mat');
    psiN = myFile.psiN(:,end-1:end, 1, ttt); %BC at north boundary
    psiN =  periodic(psiN,1,2); %make it periodic
    
    
    %%%calculate vorticity
    fprintf('Calculating vorticity \n')
    [~, dx, dudy, dy] = soseder(lonM,latM, Uvel, 'x',1);%[m/s] %clear Uvel 
    [dvdx, ~, ~, ~] = soseder(lonM,latM, Vvel, 'x',1,dx,dy);%[m/s] %clear Vvel 

    omega=dvdx - dudy; %[1/s] 
    omega = periodic(omega,1,2); %make it periodic
    
    
    %%%
    %%%calculate psi
    %%%
    
    %%% start with calculating PSI_0 (Eq 18)
    [N,M] = size(omega);
    if ttt==1
        %%initial guess
        psi=zeros(N,M); %initial guess (it matches b.c. at coast) 
        
        %%parameters for the solver
        beta=dx./dy; %grid's aspect ratio
        beta = periodic(beta,1,2); %make it periodic
        dx = periodic(dx,1,2); %make it periodic
        
        Merr = 1e-4; %error for convergence 
        Miter = 1000; %max no. of interations
 
        w = 1; %relaxation factor
        
    else 
        %psi=psi; %use psi from previous step
        
        %%parameters for the solver
        beta=dx./dy; %grid's aspect ratio
        beta = periodic(beta,1,2); %make it periodic
        dx = periodic(dx,1,2); %make it periodic
 
        %other parameters unchanged
    end
    
    fprintf('SOR begun for PSI_0... \n') 

       psi = sorramelo(psi,omega, psiN, inds_c_pd, inds_l_pd, Merr,Miter,dx,beta,w);
                       %psi,omega, north_bc, coast indeces, land indeces, parameteres
    fprintf('SOR for PSI_0 ended \n\n') 
    
    psi0 = psi; %%%PSI_0 [m3/s]
    
    %%% ISLAND RULE
    %%%
    %1) solve Eq 17 for all islands: it results in N PSI_k
      %2) calculate circulations of the PSI_k around relative islands
          %2.1) calculate normal derivatives
            %2.1.1)categorize contour pixels, define normal vector components
                   %DONE at beginning
            %2.1.2)calculate gradient and do dot product with normal vector
          %2.2) calculate line integral
        %3) solve system in Eq 19: it results in N mu_k values
          %4) calculate PSI as PSI_0+sum(mu_k*PSI_k)

          
  psikst = zeros(N,M); %for final sum of PSI_ks
  %%1) calculate PSI_k
   for isl=1:numel(islands)
       
       psis{isl} = psi0; %initialize a struct containing the solutions (use psi0 as first guess)
                  %maybe should use zeros(N,M);?
%        psis{isl}(find([neigh{isl}]>0)) = 1; %set PSI_k ==1 around k-th island
   
       fprintf(['SOR begun for PSI_',num2str(isl),'... \n'])
       
           psi = sorramelo(psi,zeros(N,M), psiN, {inds_c_pd,neigh{isl}}, inds_l_pd, Merr,Miter,dx,beta,w);
                           %psi,0, north_bc, coast indeces, land indeces, parameteres
           
           psiks.(strcat('island_',num2str(isl))) = psi;                
       fprintf(['SOR ended for PSI_',num2str(isl),' \n\n'])
       
   
   end%islands
   
   lonM=periodic(lonM,1,2);
   latM=periodic(latM,1,2);
   dy=periodic(dy,1,2);
   %%2) calculate circulations of the PSI_k around relative islands
for isl=1:numel(islands)
    fprintf(['Calculating circulations for island ',num2str(isl),'... \n'])
    %2.1) calculate normal derivatives
        %2.1.1)categorize contour pixels, define normal vector components
                      %DONE at beginning
        %2.1.2)calculate gradient and do dot product with normal vector
        %derivatives
        [psi0dx,~,psi0dy,~] = soseder(lonM,latM,psi0,'x',1);
        %vector m
        m = -psi0dy+psi0dx; %%%%%%%%%%%%% CHECK THIS!!!!
        
        %normal derivatives
        psi0dx = psi0dx.*Nx;
        psi0dy = psi0dy.*Ny;
 
        inds_c_current = find([neigh{isl}]>0);
         A = dx(inds_c_current).*dy(inds_c_current); %area elements
         circ_psi0(isl) = sum((psi0dx(inds_c_current)+psi0dy(inds_c_current)).*A); 
         circ_m(isl) = sum(m(inds_c_current).*A); clear m
         
        for jj=1:numel(islands)
             [psidx,~,psidy,~] = soseder(lonM,latM,psiks.(strcat('island_',num2str(jj))),'x',1);
              psidx = psidx.*Nx;
              psidy = psidy.*Ny;
        
              circ_.(strcat('psi_',num2str(isl))) = sum((psidx(inds_c_current)+psidy(inds_c_current)).*A); clear psidx psidy      
              
        end
        
         all_circ_k = 0;
          for jj=1:numel(islands) 
             all_circ_k = all_circ_k + circ_.(strcat('psi_',num2str(isl)));
          end
          clear circ
          all_circ(isl) = all_circ_k; clear all_circ_k
end%2)

  %%3) solve system in Eq 19: it results in N mu_k values
  fprintf(['\n Solving the system... \n'])
        %%mu_j*circ_jk+circ_0k=circ_mk   
            
        MU  = (circ_m - circ_psi0)\all_circ;   

  %%4) calculate PSI as PSI_0+sum(mu_k*PSI_k)
  fprintf(['\n Calculataing PSI... \n'])
    for isl=1:numel(islands)
    
      psikst = psikst + MU(isl)*psiks.(strcat('island_',num2str(isl)));  
      
    end%islands
%       clear psiks

     psi = psi0 + psikst; 
%      clear psikst  
    
    %%%%
    %%% TEST SOLUTION
    %%%%
    testsolution(psi0,Uvel,Vvel, lonu,latu,lonM,latM,inds_l)
    title(['$\psi_0$ Discrepancies at step ',num2str(ttt)],'interpreter','latex')
%     testsolution(psiks.island_1(:,:),Uvel,Vvel, lonu,latu,lonM,latM,inds_l)
%     title(['$\psi_1$ Discrepancies at step ',num2str(ttt)],'interpreter','latex')
    testsolution(psi,Uvel,Vvel, lonu,latu,lonM,latM,inds_l)
    title(['$\psi$ discrepancies','interpreter','latex')

    drawnow
end%time