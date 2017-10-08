% strdur_make_bstimstepfun.m
%
% Brian Goodwin, 2014-08-31
%
% Make step functions for strength-duration relationship.
%
% Minimum step function is 2*dt with maximum size of 200ms.

dt = 1e-5; % [s]
dt = dt*1e3; % [ms]

t = logspace(log10(.1),log10(100),20); % [ms]
n = numel(t);

fidpar = fopen('bstim_stepfunctionLegend.txt','w');
fprintf(fidpar,'%d\t %f\n',cat(1,1:n,t));
fclose(fidpar);

for k = 1:n
    N = ceil(t(k)/dt);
    
    wv = [linspace(0,1,N),ones(1,512)];
    
    % scaling factor (for magnetic stimulation):
    wv = fouriersolver_magnetic2efield(wv.',dt/1e3,1i,1000);
    scaling = mean(wv(2:N-1));
    
    stepfun = cat(1,zeros(floor(2/dt),1),ones(N,1)*scaling,0);
    
    fid = fopen(['bstimstepfun_',num2str(k),'.txt'],'w');
    fprintf(fid,'%f\n',stepfun);
    fclose(fid);
end