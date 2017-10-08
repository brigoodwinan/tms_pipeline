% map_V2estim_neuron.m
%
% Brian Goodwin
%
% Maps the results from cylindrical ECS electrode FEM to the neuron. The
% solution was obtained at 1V. 
%
% The neuron is placed 1cm below the ECS electrode in the middle. The axon
% hillock is located directly below the centroid of the ECS electrode.

load nrn_withaxon_20mm_z_oriented.mat % neuron *.node *.edge

coord = neuron.node./1e3; % [um] --> [mm]
cent = mean([max(coord);min(coord)]);
m = makehgtform('translate',[0,0,-10],'axisrotate',[0,1,1],pi);

nrn = xfm3d(coord,m,1);

[th,r,z] = cart2pol(nrn(:,1),nrn(:,2),nrn(:,3));

if ~exist('model','var')
    model = mphload('/Users/1773goodwib/Documents/Project_STRENGTHDURATION/Electrode_axisymmetric_ch1.mph');
end

V = mphinterp(model,'V','coord',cat(2,r,z).','recover','ppr'); % [V]
normE = mphinterp(model,'ec.normE','coord',cat(2,r,z).','recover','ppr'); % [V]

% % Uncomment the following to generate FEM text file for NEURON.
% scatter(nrn(:,1),nrn(:,2),5,get_color_for_colorbar((V-min(V))./max(V-min(V))));
% axis equal
% 
% fid = fopen('estimfem.txt','w');
% fprintf(fid,'%.9f\n',V);
% fclose(fid);