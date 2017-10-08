% fourierreconstructiondescription.m
%
% 20110713, Goodwin, Brian
%
% Performs frequency analysis for 20110711_SphereModel.mph to obtain the
% time transient with respect to the TMS waveform.

clc
clear

load waveform % waveform

dt=1e-5; % time step for 'waveform.mat' file.
n=pow2(nextpow2(length(waveform))); % number of fft points

fs=1/dt;

%% Spectral analysis
y=fft(waveform,n)/n;
Y=y(1:n/2+1).*[1,2*ones(1,n/2-1),1]; % single-sided fft
amp=abs(Y);
phase=angle(Y);
power=abs(Y).^2/n; % Power of the dft
finc=fs/n; % Freq increment
f=(0:n-1)*finc; % frequencies

% V scaled for frequencies
V=-5.516976426297033E-8-0.012208966059577678i; % Solution from point in mesh
fsol=1000; % Hz - frequency of voltage solution
Vphi=angle(V); % phase of V solution
F=f(1:n/2+1); % positive frequencies
V=V*F/1e3; % frequency scaling
V=V.*amp; % amplitude scaling from tms waveform
mag=abs(V); % Magnitude of V
phi=Vphi+phase; % new phase 

sol2=complex(mag.*cos(phi),mag.*sin(phi));
sol3=[sol2 fliplr(conj(sol2(2:n/2)))]*n/2;
sol3(n/2+1)=abs(sol3(n/2+1));
new=-ifft(sol3);

new=new/max(new);


%% Plotting results
load reduceddata % current dt t tta
t=(0:n-1)*dt;
recorded=tta(6,:);
recorded=recorded(1:n)/max(recorded);

subplot(4,1,1)
plot(t*1e3,waveform,'k','Linewidth',2)
xlabel 'Time (ms)'
ylabel({'Current Waveform','Normaized Units'})
axis([0 1.3 -.4 1.2])
subplot(4,1,2)
plot(F/1e3,amp/max(amp),'k','Linewidth',2)
xlabel 'Frequency (kHz)'
ylabel({'Current Waveform';'Frequency Domain'})
axis([0 30 0 1.2])
subplot(4,1,3)
plot(F/1e3,mag/max(mag),'k','Linewidth',2)
xlabel 'Frequency (kHz)'
ylabel({'Voltage Solution';'Frequency Response'})
axis([0 30 0 1.2])
subplot(4,1,4)
plot(t*1e3,new,'k','Linewidth',2)
xlabel 'Time (ms)'
ylabel({'Voltage Solution';'Time Response'})
axis([0 1.3 -.4 1.2])

%% plotting for Grant

o=zeros(1,floor(.01/dt/6));
waveform=[o,waveform,o];
new=[o,new,o];
t=[(-length(o):-1)*dt,(0:length(new)-length(o)-1)*dt];

figure(2)
subplot(2,1,1)
plot(t*1e3,waveform,'k','Linewidth',1)
xlabel 'Time (ms)'
ylabel 'normalized units'
% axis([0 1.3 -.4 1.2])

subplot(2,1,2)
plot(t*1e3,new,'k','Linewidth',1)
xlabel 'Time (ms)'
% axis([0 1.3 -.4 1.2])