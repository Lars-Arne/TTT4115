function fFir = Freq2Time(F)
% F is the vector of frequency values between (0:N-1)/N which correspond to
% 0:(1-1/N)*Fs Hz with a gap of Fs/N Hz sampling distance in frequency or
% 1/N in the normalized frequency range.

% N should be an odd number
if size(F,2)>size(F,1) % Make sure if F is a column vector
    F = F';
end

if mod(length(F),2)==0
    error('Use a filter with odd number of taps to get better results');
end

NFir = length(F);
fFir = ifft(F, NFir); % Estimated unit pulse response
fFir=real([fFir((NFir+1)/2+1:NFir,:);fFir(1:(NFir+1)/2,:)]);% adjusting the 
% correct position of samples
end

