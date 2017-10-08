% strdur_make_stepfun.m
%
% Brian Goodwin, 2014-08-31
%
% Make step functions for strength-duration relationship.
%
% Minimum step function is 0.1ms with maximum size of 100ms.

dt = 1e-5; % [s]
dt = dt*1e3; % [ms]

t = logspace(log10(.1),log10(100),20); % [ms]
n = numel(t);

fidpar = fopen('./stepfunctionLegend.txt','w');
fprintf(fidpar,'%d\t %f\n',cat(1,1:n,t));
fclose(fidpar);

for k = 1:n
    N = ceil(t(k)/dt);
    
    stepfun = cat(1,zeros(floor(10/dt),1),ones(N,1),0);
    
    fid = fopen(['./input/stepfun_',num2str(k),'.txt'],'w');
    fprintf(fid,'%f\n',stepfun);
    fclose(fid);
end