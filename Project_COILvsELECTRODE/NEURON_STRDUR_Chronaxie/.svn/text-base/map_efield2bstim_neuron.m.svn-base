% map_efield2bstim_neuron.m
%
% Brian Goodwin
%
% Maps the results from single circular coil FEM to the neuron. The
% solution was obtained at 1A and 1000Hz. The FEM is scaled so that it is
% as if the solution was obtained at 1000A for accuracy purposes in NEURON
% HOC.

load nrn_withaxon_20mm_z_oriented.mat % neuron *.node *.edge

coord = neuron.node./1e3; % [um] --> [mm]
cent = mean([max(coord);min(coord)]);
m = makehgtform('translate',[45,0,-10],'axisrotate',[0,1,1],pi,'translate',-cent);

nrn = xfm3d(coord,m,1);

[th,r,z] = cart2pol(nrn(:,1),nrn(:,2),nrn(:,3));

if ~exist('model','var')
    model = mphload('/Users/1773goodwib/Documents/Project_STRENGTHDURATION/MultiturnCoil_axisymmetric_ch1.mph');
end

e = mphinterp(model,'mef.Ephi','phase','90','coord',cat(2,r,z).'); % [V/m]

ex = e.*-sin(th).';
ey = e.*cos(th).';

V = map_efield2V_nrn(nrn,neuron.edge,ex,ey,zeros(1,length(ex)));

% scatter(nrn(:,1),nrn(:,2),1,get_color_for_colorbar((V-min(V))./max(V-min(V))));

fid = fopen('bstimfem.txt','w');
fprintf(fid,'%.9f\n',V*1000);
fclose(fid);