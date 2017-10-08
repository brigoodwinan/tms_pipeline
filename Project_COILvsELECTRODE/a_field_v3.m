% a_field.m
%
% Goodwin, Brian 2012-02-24
%
% Edits:
% 2013-08-15 - v2
% 2014-08-26 - v3
%
% Reshapes the Efield into column vectors.
%
% Magnetic vector potential calculation for testing. This calculation would
% ultimately be utilized in the FEM.
% This program specifically calculates the maximum mangetic vector
% potential (A) that would be observed in free space.
%
% VERSION NOTES:
% v2
% gives the option of calculating the A-field along a line on the x-axis
% underneath the coil for comparison purposes.
%
% v3
% general edits.

% Coil located at [0,0,0] with the primary stimulation direction pointing in the y-direction
load /Users/1773goodwib/Documents/Project_HEAD_MODEL/coil.mat % coil; *.node *.data
% Magnetic permeability of free space
mu_0 = 4*pi*1E-7; % [N/A^2]  = [H/m]  = [J/(A^2m)]
mu_r = 1;

% Line. These need to be in [m]
x = (-.1:.001:.1).';
y = zeros(length(x),1);
z = ones(length(y),1)*-.0035;
N = numel(x);

% % Create Plane
% [X,Y] = meshgrid(-140:3:140,-100:3:100);
% X = X./1e3;
% Y = Y./1e3;
% [xm,xn] = size(X);
% x = reshape(X,[],1);
% y = reshape(Y,[],1);
% z = ones(length(y),1)*-.01;

% % Create Box (for curl calc)
% [X,Y,Z] = meshgrid(-5:.1:5,-5:.1:5,-10:.1:0);
% X = X./1e3;
% Y = Y./1e3;
% Z = Z./1e3;
% [xm1,xm2,xm3] = size(X);
% N = numel(X);
% x = reshape(X,N,1);
% y = reshape(Y,N,1);
% z = reshape(Z,N,1);

% % Single point evaluation
% x = 0;
% y = 0;
% z = -.01;

n = length(coil.node);

Qp = mat2cell(coil.node,n,[1,1,1]);
Q = mat2cell(coil.data,n,[1,1,1]);

N = cat(1,ones(floor(N/200),1)*200,rem(N,200));
x = mat2cell(x,N,1);
y = mat2cell(y,N,1);
z = mat2cell(z,N,1);
N = num2cell(N);
Nlen = length(N);
A = cell(Nlen,1);

if isPoolOpen
    parfor k = 1:Nlen
        % Magnetic vector potential calculation (A-field)
        M = mu_0*mu_r/4/pi./sqrt(...
            bsxfun(@minus,x{k},repmat(Qp{1}',N{k},1)).^2+...
            bsxfun(@minus,y{k},repmat(Qp{2}',N{k},1)).^2+...
            bsxfun(@minus,z{k},repmat(Qp{3}',N{k},1)).^2);
        
        A{k} = cat(2,M*Q{1},M*Q{2},M*Q{3}); % [V-s/m] [Wb/m] [H-A/m]
    end
else
    for k = 1:Nlen
        % Magnetic vector potential calculation (A-field)
        M = mu_0*mu_r/4/pi./sqrt(...
            bsxfun(@minus,x{k},repmat(Qp{1}',N{k},1)).^2+...
            bsxfun(@minus,y{k},repmat(Qp{2}',N{k},1)).^2+...
            bsxfun(@minus,z{k},repmat(Qp{3}',N{k},1)).^2);
        
        A{k} = cat(2,M*Q{1},M*Q{2},M*Q{3}); % [V-s/m] [Wb/m] [H-A/m]
    end
end
A = cell2mat(A);
x = cell2mat(x);

% Ax = reshape(A(:,1),xm1,xm2,xm3);
% Ay = reshape(A(:,2),xm1,xm2,xm3);
% Az = reshape(A(:,3),xm1,xm2,xm3);
% 
% [bx,by,bz,bav] = curl(X,Y,Z,Ax,Ay,Az);

a = vectormag(A);
e = a./100e-6;% For Salinas et al. (2007) comparison. For realistic E-field use: a*16193; % constant has units of s^-1: max(diff(currentwave))./1e-5

% 1D Plot
figure
plot(-x*1e3,e,'k')
ylim([0,1.1*340])

% % 2D Vector
% figure
% quiver(1e3*x,1e3*y,-A(:,1),-A(:,2),'k')
% axis equal

% % 2D Contour Plot
% a = a./max(a);
% a = reshape(a,xm,xn);
% figure
% contour(X*1e3,Y*1e3,a,.1:.1:1)
% axis equal

% fig = figure;
% ax = axes('Parent',fig,'fontname','arial','fontsize',12);
% box(ax,'off');
% plot(x*1e3,a./max(a)*350,'k','linewidth',3)
% xlim(ax,[-100,100]);
% ylim(ax,[0,400]);
% xlabel(ax,'X Direction [mm]','fontname','arial','fontsize',12);
% ylabel(ax,'Electric Field [V/m]','fontname','arial','fontsize',12)

%%% Writing files
% Creates text file for COMSOL interpolation function
if 0
    matrix = [x,y,z,A];
    fid = fopen(['comsolfieldinput_',datestr(now,30),'.txt'],'w');
    fprintf(fid,'%8.8f %8.8f %8.8f %8.8f %8.8f %8.8f\n',matrix');
    fclose(fid);
end

% Write Field points and data for SCIRun
if 0
    fid = fopen('fieldpts.pts','w');
    fprintf(fid,'%12.12f %12.12f %12.12f\n',fieldpts');
    fclose(fid);
    fid = fopen('fielddata_bf.txt','w');
    fprintf(fid,'%12.12f %12.12f %12.12f\n',bf'); % writes magnetic field (not vector potential field).
    fclose(fid);
    fid = fopen('fielddata_af.txt','w');
    fprintf(fid,'%12.12f %12.12f %12.12f\n',A'); % writes magnetic vector potential field).
    fclose(fid);
    
    % Write coil to file for SCIRun
    fid = fopen('coil_with_height.pts','w');
    fprintf(fid,'%8.8f %8.8f %8.8f\n',coil.node);
    fclose(fid);
    fid = fopen('coil_with_height.txt','w');
    fprintf(fid,'%8.8f %8.8f %8.8f\n',coil.data);
    fclose(fid);
end
