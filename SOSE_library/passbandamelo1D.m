function [bpassed,FTb,FT,f] = passbandamelo1D(field,fs,fhigh,flow)
%PASSBANDAMELO1D passbands a signal.
%
%   [BPASSED, FT, F] = passbandamelo(FIELD,FS,FHIGH,FLOW)
%
%   FIELD must be a 1D array.
%
%   FS is the sampling frequency [Hz]. 
%      Ex: fs=1/(3*86400) is a 3-days sampling frequency.
%
%   Set FHIGH or FLOW are the high and low cut-off requencies 
%   respectively.
%   Set either of them to NaN if they are not needed.
%
%   
%   BPASSED is the band-passed signal. It has the same mean of 
%   the original signal.
%   FTb is the (one-sided) power spectrum of the band-passed signal.
%   FT is the (one-sided) power spectrum of the original signal.
%   F is the frequency vector.

%   Shame on Andrea Costa version 2017.10.09 in case of bugs.



%ensure that N is multiple of 2
L = numel(field);
newL = 2*floor(L/2);
field = field(1:newL);


FT=fft(field-nanmean(field)) /numel(field); %note the normalizzation

FT=FT(1:numel(FT)/2+1);    %one-sided spectrum
FT=2*FT;   %conserve energy


f = fs*(1:numel(FT)) /numel(field); %normalized frequencies

if isnan(flow)*isnan(fhigh)~=1 %filtering
%bandpassing
    FTb=FT;
    if ~isnan(fhigh)
        FTb(f>fhigh)=0; %rectangular window
    end

    if ~isnan(flow)
         FTb(f<flow)=0;
    end


    bpassed=ifft([FTb; flip(FTb(1:end-2))] *numel(field)/2  ,'symmetric');
    %note the size(field,3)/2 rescaling. Needed for the reconstructed signal to
    %have the right amplitude

    bpassed = bpassed+nanmean(field); %re-adding mean of original signal


elseif isnan(flow)*isnan(fhigh)==1 %NO filtering !!
    
    bpassed = -9999;
    FTb = -999;
 
end