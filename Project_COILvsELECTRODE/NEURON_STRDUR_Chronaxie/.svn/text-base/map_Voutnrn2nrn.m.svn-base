% map_Voutnrn2nrn.m
%
% Brian Goodwin, 2014-08-04
%
% This was designed to generate a matrix to transform intracellular
% Voltages to the extracullular neuron compartments. The *hoc code does
% not record the intracellular voltage at every compartment that the
% extraceulluar voltage is imposed upon. This matrix (s) would be used
% primarily for visualization purposes.

nnrns = 6111; % this is the number that "run.orient_cell_body" gave.

suff = '';

% Neuron files
% File extensions
ext.i2m = '_i2m';
ext.srnfld = '_scirunfld';
ext.nrns = ['_',num2str(nnrns)];
ext.merged = '_merged';
ext.unmerged = '_unmerged';
ext.mat = '.mat';
ext.txt = '.txt';
ext.pts = '.pts';

nrnf.neuronpath = './';
nrnf.corticalpath = '../../../Project_CORTICAL_MODEL/JL_CorticalModel/';
nrnf.nrnivoutpath = '/Volumes/Macintosh HD 2/Brian/NRNIV_Results/';
nrnf.pop_nrns_cellstruct = [nrnf.corticalpath,'pop_withaxons_20mm_export_comsol_sol_cellstruct',ext.nrns,ext.mat];
nrnf.pop_nrns_scirunviz = [nrnf.corticalpath,'pop_withaxons_20mm_export_comsol_sol',ext.srnfld,ext.nrns,ext.mat];
nrnf.edge_connections = [nrnf.neuronpath,'nrn_edgeconn_withaxon_20mm',ext.mat];
nrnf.placements = [nrnf.corticalpath,'nrn_placements_',num2str(nnrns),suff,ext.mat];
nrnf.nrnpialindices = [nrnf.corticalpath,'nrnplace_GMindices_',num2str(nnrns),suff,ext.mat];
nrnf.coords_raw = [nrnf.neuronpath,'coords_withaxon_ax20mm',ext.pts];
nrnf.extracellularsecs = [nrnf.neuronpath,'coords_section_orders',ext.txt];
nrnf.recordedsecs = [nrnf.neuronpath,'VoutSections_ax20mm',ext.txt];
nrnf.recordedpts = [nrnf.neuronpath,'VoutLocations_ax20mm',ext.txt];
nrnf.nrn_z_oriented = ['nrn_withaxon_20mm_z_oriented',ext.mat];
nrnf.Vout_nrn_z_oriented = ['VoutNeuron_withaxon20mm_z_oriented',ext.mat];
nrnf.efield_mult_matrix = ['efield_mult_matrix_ax20mm',ext.mat];
nrnf.threshold_results = [nrnf.nrnivoutpath,'JL_thresholdsV_1thru79443_all',ext.mat];

%% CODE %%
load(nrnf.Vout_nrn_z_oriented) % voutneuron *.node; vindex *.soma *.iseg *.hill *.axon *.cellbody
load(nrnf.nrn_z_oriented) % neuron *.node *.edge; index *.axon *.cellbody *.hill *.soma
load(nrnf.edge_connections) % edge; {secnames}
fid = fopen(nrnf.recordedsecs,'r');
recsecs = textscan(fid,'%s');
fclose(fid);

% Matrix needs to be length(neuron.node)-by-length(voutneuron.node)
nrnlen = size(neuron.node,1);
vnrnlen = size(voutneuron.node,1);
I = 1:nrnlen; % I is the corresponding neuron node
J = zeros(1,nrnlen);
for k = 1:nrnlen
    loc = strfind(secnames{1}{k},'_');
    loc = loc(end);
    
    J(k) = find(strncmpi(secnames{1}{k},recsecs{1},length(secnames{1}{k}(1:loc-1))));
    
end
s = sparse(I,J,1,nrnlen,vnrnlen);

README = 'Map intracellular voltage to extracellular node points by v_nrn = s*v_recnrn.';
save([nrnf.corticalpath,'matrix_map_Vintra2Vextra',ext.mat],'s','README');
save([nrnf.neuronpath,'matrix_map_Vintra2Vextra',ext.mat],'s','README');

