% analyze_somaV_vs_thresh.m
%
% Brian Goodwin 2014-09-03
%
% Analyzes data from soma voltage vs threshold experiment.

do.gatherdata = 0;
do.plotstuff1 = 0;
do.plotstuff2 = 0;
do.analyzebi72m1_neurons = 1;

load nrn_withaxon_20mm_z_oriented.mat % neuron; index *.axon *.cellbody *.hill *.soma
load matrix_map_Vintra2Vextra.mat % s
t = 7:.05:17.5;
ncoords = 690;
ntime = 212;
if do.gatherdata
    %%% FYI 0.01 is actually .005 and
%     amps = [.01,.05,.1,.12,.13,.14,.15,.2,.25,.3,.35];
    amps = .1:.05:.35;
    
    n = length(amps);
    v_mono = cell(n,1);
    v_bi = cell(n,1);
    soma_mono = zeros(n,ntime-1);
    soma_bi = zeros(n,ntime-1);
    for k = 1:n
        tmp = read_nrniv_bin_vout(['./output/somaV_vs_thresh_monophase_',num2str(amps(k),'%.2f'),'.bin'],ncoords,ntime);
        v_mono{k} = s*tmp;
        soma_mono(k,:) = v_mono{k}(index.soma,:);
        
        tmp = read_nrniv_bin_vout(['./output/somaV_vs_thresh_biphase_',num2str(amps(k),'%.2f'),'.bin'],ncoords,ntime);
        v_bi{k} = s*tmp;
        soma_bi(k,:) = v_bi{k}(index.soma,:);
    end
end
if do.plotstuff1
    
    % at t = 9.95, index = 60;
    tindex = 60;
    
    data = cat_txt_files('./output/somaV*txt',3); % num mono bi
    amps(1) = .005;
    
    figure
    plot(t,soma_mono)
    title 'mono'
    xlim([9,17])
    
    
    figure
    plot(t,soma_bi)
    title 'bi'
    xlim([9,17])
    
    monothresh = data(6:end,2);
    monothresh = monothresh./max(monothresh);
    bithresh = data(6:end,3);
    bithresh = bithresh./max(bithresh);
    
    figure
    plot(soma_mono(6:end,tindex),monothresh,'^k-')
    hold on
    plot(soma_bi(6:end,tindex),bithresh,'ok-')
    hold off
    legend('mono','bi')
    
end

if do.plotstuff2
    tindex = 60;
    
    nlen = size(soma_bi,1);
    bilat = zeros(nlen,1);
    monolat = zeros(nlen,1);
    for k = 1:nlen
        bilat(k) = t(find(soma_bi(k,:)>0,1,'first'))-10;
        monolat(k) = t(find(soma_mono(k,:)>0,1,'first'))-10;
    end
    bilat(:,2) = soma_bi(:,tindex);
    monolat(:,2) = soma_mono(:,tindex);
    
    figure
    scatter(bilat(:,2),bilat(:,1),40,'ok')
    hold on
    scatter(monolat(:,2),monolat(:,1),40,'^k')
    hold off
    ylim([2,5])
end

if do.analyzebi72m1_neurons
    load nu357_somav.mat % nu357
    load nu3236_somav.mat % nu3236
    
    n = 690; % number of recorded compartments
    
    fempref = {'357';'3236'};
    
    % Monophase 357
    fname = dir(['./output/v',fempref{1},'*monophase*bin']);
    
    vname = {fname.name}.';
    
    vlen = length(vname);
    v357mono = zeros(vlen,1);
    for k = 1:vlen
        tmp = read_nrniv_bin_vout(['./output/',vname{k}],n);
        v357mono(k) = tmp(1,55);
    end
    
    % Biphase 357
    fname = dir(['./output/v',fempref{1},'*biphase*bin']);
    
    vname = {fname.name}.';
    
    vlen = length(vname);
    v357bi = zeros(vlen,1);
    for k = 1:vlen
        tmp = read_nrniv_bin_vout(['./output/',vname{k}],n);
        v357bi(k) = tmp(1,55);
    end
    
    % Monophase 3236
    fname = dir(['./output/v',fempref{2},'*monophase*bin']);
    
    vname = {fname.name}.';
    
    vlen = length(vname);
    v3236mono = zeros(vlen,1);
    for k = 1:vlen
        tmp = read_nrniv_bin_vout(['./output/',vname{k}],n);
        v3236mono(k) = tmp(1,55);
    end
    
    % Biphase 3236
    fname = dir(['./output/v',fempref{2},'*biphase*bin']);
    
    vname = {fname.name}.';
    
    vlen = length(vname);
    v3236bi = zeros(vlen,1);
    for k = 1:vlen
        tmp = read_nrniv_bin_vout(['./output/',vname{k}],n);
        v3236bi(k) = tmp(1,55);
    end
    
    % Plotting
    figure
    plot(v357mono,nu357(:,2)./max(nu357(:,2)),'-xr')
    hold on
    plot(v357bi,nu357(:,3)./max(nu357(:,3)),'-xb')
    hold off
    figure
    plot(v3236mono,nu3236(:,2)./max(nu3236(:,2)),'-^r')
    hold on
    plot(v3236bi,nu3236(:,3)./max(nu3236(:,3)),'-^b')
    hold off
    legend(fempref{1},fempref{1},fempref{2},fempref{2})
    
end