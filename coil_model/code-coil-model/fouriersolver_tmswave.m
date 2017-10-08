% fouriersolver_tmswave.m
%
% 2013-12-26, Brian Goodwin
%
% Uses the Fourier Solver to solve for the time dependent waveform of the
% voltage at a node in our FEM.
%
% Assumes that the solution frequency of the FEM was 1kHz.
%
% This prevents us from having to do the extra work of converting
% the voltage at each node to the time domain. The solution can be
% sent to NRNIV directly since the TMS wave. The following two files
% have been generated using this m-script:
% 
% */Project_NEURON_MODEL/AmatrudoWeaver_GoodwinMod/modified_amatrudo_goodwin/tms_wave10ms_MagstimRapid.txt
% */Project_NEURON_MODEL/AmatrudoWeaver_GoodwinMod/modified_amatrudo_goodwin/tms_wave10ms_Magstim200.txt

load ./tms_ms200_wave.mat % currentwave
% load ./tms_msrapid_wave.mat % currentwave

dt = 1e-5;
fs = 1/dt;

wv = currentwave.';
N = length(wv);
n = nextpow2(N);
if n<512, n = 512; end

N = n/2+1;
H = fft(wv,n)/n;
H = H(1:N).*[1;2*ones(N-2,1);1]; % single sided.

% Result in mesh that corresponds to normalized "currentwave".
% This way the resulting waveform will accurately scale the voltage
% solution as obtained from the FEM results.
v = 1i;

% Assume solution generated at 1000Hz
fsol = 1000; % Hz - frequency of voltage solution
phi = angle(v)+angle(H); % new phase 

finc = fs/n; % Freq increment
f = (0:n-1)*finc; % frequencies

F = f(1:n/2+1)'; % positive frequencies
v = abs(v).*F/fsol; % frequency scaling
v = v.*abs(H);

v = complex(v.*cos(phi),v.*sin(phi));
v(1) = real(v(1));
v(end) = real(v(end));
v_fft = v;
v(2:end-1) = v(2:end-1)/2;
v = [v;flipud(conj(v(2:end-1)))]*n;
v_t = ifft(v);

v_t = v_t(1:length(wv));

figure
axes1 = axes('FontSize',12,'FontName','Arial');
box(axes1,'off');
hold(axes1,'all');
stem((0:length(wv)-1).*dt,wv,'.k','markersize',8)
axis([0,.5e-3,0,1.2]);

figure
axes1 = axes('FontSize',12,'FontName','Arial');
box(axes1,'off');
hold(axes1,'all');
stem((0:length(v_t)-1).*dt,v_t,'.k','markersize',8)
axis([0,.5e-3,-1,3.25]);

figure
axes1 = axes('FontSize',12,'FontName','Arial');
box(axes1,'off');
hold(axes1,'all');
stem((0:length(v_fft)-1).*finc,abs(v_fft),'.k','markersize',8)
axis([0,20e3,0,.1]);

figure
axes1 = axes('FontSize',12,'FontName','Arial');
box(axes1,'off');
hold(axes1,'all');
stem((0:length(H)-1).*finc,abs(H),'.k','markersize',8)
axis([0,20e3,0,.14]);