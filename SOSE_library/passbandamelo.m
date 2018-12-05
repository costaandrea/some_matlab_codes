function [bpassed,FTb,FT,f] = passbandamelo(field,fs,fhigh,flow)
%PASSBANDAMELO passbands a signal.
%
%   [BPASSED, FTb, FT, F] = passbandamelo(FIELD,FS,FHIGH,FLOW)
%
%   FIELD must be a 3D matrix with time as third dimension.
%
%   FS is the sampling frequency [Hz]. 
%      Ex: fs=1/(3*86400) is a 3-days sampling frequency.
%
%   Set FHIGH or FLOW are the high and low cut-off requencies 
%   respectively. Freqs higher than FHIGH are suppressed and vice-versa.
%   Set either of them to NaN if they are not needed.
%
%   
%   BPASSED is the band-passed signal. It has the same mean of 
%   the original signal.
%   FTb is the (one-sided) power spectrum of the band-passed signal.
%   FT is the (one-sided) power spectrum of the original signal.
%   F is the frequency vector.

%   Shame on Andrea Costa version 2017.09.28 in case of bugs.



%ensure that N is multiple of 2
L = size(field,3);
newL = 2*floor(L/2);
field = field(:,:,1:newL);


FT=fft(field-nanmean(field,3)  ,[],3) /size(field,3); %note the normalizzation

FT=FT(:,:,1:size(FT,3)/2+1);    %one-sided spectrum
FT=2*FT;   %conserve energy


f = fs*(1:size(FT,3)) /size(field,3); %normalized frequencies

if isnan(flow)*isnan(fhigh)~=1 %filtering
%bandpassing
    FTb=FT;
    if ~isnan(fhigh)
        FTb(:,:, f>fhigh)=0; %rectangular window
    end

    if ~isnan(flow)
         FTb(:,:, f<flow)=0;
    end


    bpassed=ifft(cat(3,FTb,flip(FTb(:,:,1:end-2),3)) *size(field,3)/2  ,[],3,'symmetric');
    %note the size(field,3)/2 rescaling. Needed for the reconstructed signal to
    %have the right amplitude

    bpassed = bpassed+nanmean(field,3); %re-adding mean of original signal


elseif isnan(flow)*isnan(fhigh)==1 %NO filtering !!
    
    bpassed = -9999;
    FTb = -999;
 
end