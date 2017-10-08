% strdur_analysis.m
%
% Brian Goodwin 2014-09-02
%
% Analysis of strength-duration results.

% All variables have *.dur *.str structures.
load strdur_ecs_result.mat % ecs
load strdur_bstim_result.mat % bstim
load strdur_currinj_result.mat % currinj

% Solution for bstim in FEM was 1000Hz and 1000A
% the mean E-field magnitude was 14.3401 V/m with a std of 0.0416 V/m

% Solution for ECS FEM was at 1V
% the mean E-field magnitude was 16.7795 V/m with a std of 4.405 V/m

it = length(bstim.dur);

% % The following would be to find the stimulating current through the
% % coil
% dt = 1e-5;
% for k = 1:it
%     n = ceil(bstim.dur(k)/1e3/dt);
%     wv = fouriersolver_magnetic2efield(linspace(0,1,n).',dt,1i,1000);
%     scaling = mean(wv(2:n-5));
%     bstim.str(k) = bstim.str(k)./scaling;
% end

figure
semilogx(bstim.dur,bstim.str,'ok-',ecs.dur,ecs.str,'or-')

bstim.str = bstim.str*14.3401; % [V/m]
bstim.std = bstim.str*0.0416./14.3401; % [V/m]
chr = min(bstim.str)*2;
tmp1 = find(bstim.str>chr,1,'last');
tmp2 = find(bstim.str<chr,1,'first');
tmp = (chr-bstim.str(tmp2))/(bstim.str(tmp1)-bstim.str(tmp2));
bstim.tau = (bstim.dur(tmp1)-bstim.dur(tmp2))*tmp+bstim.dur(tmp2);

disp(bstim.tau)

ecs.str = ecs.str*16.7795; % [V/m]
ecs.std = ecs.str*4.405./16.7795; % [V/m]
chr = min(ecs.str)*2;
tmp1 = find(ecs.str>chr,1,'last');
tmp2 = find(ecs.str<chr,1,'first');
tmp = (chr-ecs.str(tmp2))/(ecs.str(tmp1)-ecs.str(tmp2));
ecs.tau = (ecs.dur(tmp1)-ecs.dur(tmp2))*tmp+ecs.dur(tmp2);
disp(ecs.tau)

chr = min(currinj.str)*2;
tmp1 = find(currinj.str>chr,1,'last');
tmp2 = find(currinj.str<chr,1,'first');
tmp = (chr-currinj.str(tmp2))/(currinj.str(tmp1)-currinj.str(tmp2));
currinj.tau = (currinj.dur(tmp1)-currinj.dur(tmp2))*tmp+currinj.dur(tmp2);
disp(currinj.tau)

figure
semilogx(bstim.dur,bstim.str./min(bstim.str),'ok-',ecs.dur,ecs.str./min(ecs.str),'xk-',currinj.dur,currinj.str./min(currinj.str),'^k-')
legend('COIL','ECS','Somatic INJ')

figure
semilogx(bstim.dur,bstim.str,'ok-',ecs.dur,ecs.str,'or-')

