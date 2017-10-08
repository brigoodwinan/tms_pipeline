function A = calcMagneticVectorPotential(p,m)
% A = calcMagneticVectorPotential(p,m)
%
% Goodwin, Brian 2012-02-24
%
% Versions:
% 2013-08-15 - v2
% 2014-08-26 - v3
% 2017-06-12 - v4
%
% REQUIRES: https://github.com/brigoodwinan/matlab-general
%
% Calculates the A-field (magnetic vector potential) using the figure-8 TMS
% coil model ('./coilmat') at given points in 3-D space and a given
% transformation matrix (4-by-4)
%
% INPUTS:
% p: n-by-3 column vector of points ([x,y,z])
%
% m: (optional) 4-by-4 transformation matrix having the form that would 
%      come from executing (e.g.):
%      >> makehgtform('translate',[1,3,5],'axisrotate',[4,6,1],pi/8)
%
% OUTPUT:
% A: n-by-3 column vector of the magnetic vector potential field calculated
%      at points, p; i.e., 
%      >> a_field = struct('node',p,'data',A); % scirun format
%
% E-field equation, which requires waveform (after A-field is calculated):
% E = -nabla(phi)-d(A)/dt

% Coil located at [0,0,0] with the primary stimulation direction pointing in the y-direction
load ./coil.mat % coil; *.node *.data

if nargin>1
    coil.node = xfm3d(coil.node,m,1);
    coil.data = xfm3d(coil.data,m,0);
end

% Magnetic permeability of free space
mu_0 = 4*pi*1E-7; % [N/A^2]  = [H/m]  = [J/(A^2m)]
mu_r = 1;

n = length(coil.node);

Qp = mat2cell(coil.node,n,ones(1,3));
Q = mat2cell(coil.data,n,ones(1,3));

% What is this?
N = size(p,1);
N = cat(1,ones(floor(N/200),1)*200,rem(N,200));

p = mat2cell(p,N,ones(1,3));

N = num2cell(N);
Nlen = length(N);
A = cell(Nlen,1);

% This for-loop approach conserves RAM
for k = 1:Nlen
    % Magnetic vector potential calculation (A-field)
    M = mu_0*mu_r/4/pi./sqrt(...
        bsxfun(@minus,p{k,1},repmat(Qp{1}',N{k},1)).^2+...
        bsxfun(@minus,p{k,2},repmat(Qp{2}',N{k},1)).^2+...
        bsxfun(@minus,p{k,3},repmat(Qp{3}',N{k},1)).^2);
    
    A{k} = M*cell2mat(Q); % [V-s/m] [Wb/m] [H-A/m]
end

A = cell2mat(A);
return