function f = fouriersolver(F)
% f = fouriersolver(F)
%
% Brian Goodwin 2014-04-05
%
% Provide the frequency domain solution and the result will be the time
% domain solution. The function assumes that F is single sided. The length
% of the input must be N/2+1, where N is length of the DFT. i.e. the DFT is
% symmetric and the result will be real.
%
% This function assumes the double sided DFT to be of length 2^n where n is
% an integer.
%
% i.e. if F is 65pts long, it returns a 128 point ifft (time domain signal)
%
% INPUTS:
% F: single sided DFT (MxN), where the DFT is assumed to have a length of
%      2*(N-1). Each DFT in the matrix must run along the second dimension.
%      In other words, they must be arrays, not column vectors.
% 
% OUTPUTS:
% f: time series where ifft is calculated along the second dimension of the
%     F vector.

n = size(F,2);
N = 2*(n-1);
f = cat(2,F(:,1),F(:,2:n-1)/2,real(F(:,n)),conj(fliplr(F(:,2:n-1))/2))*N;
f = ifft(f,[],2);