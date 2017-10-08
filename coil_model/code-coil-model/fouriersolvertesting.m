% fouriersolvertesting.m
%
% 2013-12-26, Brian Goodwin
%
% This just tests the validity of the fourier solver to ensure that I am
% performing the correct operations. For example, when the input wvform
% is a square wv and the electrode sees the square wv, the electrode
% records a scaled version of it. So when we do the ifft(), we should still
% see a square wv at the electrode. This code tests this for verification
% purposes.

% Scaling value:
k = 0.5;
n = 128;
wv = zeros(n,1);
wv(20:30) = 1;

N = n/2+1;
H = fft(wv,n)/n;
H = H(1:N).*[1;2*ones(N-2,1);1];

v = H.*k;

v(1) = real(v(1));
v(end) = real(v(end));
v(2:end-1) = v(2:end-1)/2;
v = [v;flipud(conj(v(2:end-1)))]*n;
v_t = ifft(v);

stem([wv,v_t],'.','markersize',15)