function psi = sorramelo(psi,omega, bc_open, cinds,inds_l_pd, Merr,Miter,dx,beta,w)


if iscell(cinds) && numel(cinds)==2 %PSI_k
    inds_c_current = find([cinds{2}]>0); %indeces of current coast
    inds_c_others = setdiff([cinds{1}], inds_c_current); %indeces of other coasts (all - current)
    valc = 1; %values at coast for current islands
    
elseif ~iscell(cinds) %PSI_0
    inds_c_current = cinds;
    inds_c_others = [];   %indeces of other coasts not needed
    valc = 0; %values at coast for current islands
end
clear cinds


[N,M] = size(psi);

for iter=1:Miter  %sor/Jacobi solving the problem \nabla2 f = Q
    
    Err=0.0; %error
    psi_old=psi; %to use for convergence test (psi at previous time step)
    
    for jj=3:M-2
        for ii=3:N-2
            %%fourth order formula
            psi(ii,jj) = (1-w)*psi(ii,jj) + w*...
                ( 16*( psi(ii+1,jj) + psi(ii-1,jj) ) - psi(ii+2,jj) - psi(ii-2,jj) +...
                beta(ii,jj)^2*( -psi(ii,jj+2) + 16*(psi(ii,jj+1) + psi(ii,jj-1)) -psi(ii,jj-2) ) -...
                12*dx(ii,jj)^2*omega(ii,jj) )...
                /(30*(1+(beta(ii,jj)^2)));
            
            %%second order formula  (to use this gotta modify for cycle extremes and periodic inputs
            %%% (use 2:M-1 and 2:N-1 and periodic ,1,1)
            %          psi(ii,jj) = (1-w)*psi(ii,jj) + ...
            %              w*( (psi(ii+1,jj)+beta(ii,jj)^2*psi(ii,jj+1)+psi(ii-1,jj)+beta(ii,jj)^2*psi(ii,jj-1)) -...
            %              2*dx(ii,jj)^2*omega(ii,jj) )...
            %              /(2*(1+(beta(ii,jj)^2)));
            
        end
    end
    
    %%%OPEN BOUNDARY BC
    psi(:,end-1:end) = bc_open;
    
    %%%BC AT COASTS
    psi(inds_c_current)= valc; %constant at current coasts
    psi(inds_c_others) = 0.0; %constant at other coasts
    
    psi(inds_l_pd) = 0.0; %zero on all continents
    
    
    %%%check for convergence
    Err(end+1) = sqrt( sum(sum( (psi - psi_old).^2 )) ); %psi - psi_p1
    
    if Err(end) <= Merr
%         Err(iter)
%         iter
%         figure;plot(Err)
          break; %stop if iteration has converged
          
    elseif Err(end)>Merr && iter==Miter
        warning('convergence error')
        
    end
    
end%iter