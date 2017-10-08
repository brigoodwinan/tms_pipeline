% a_field.m
%
% Goodwin, Brian 2012-02-24
%
% Versions:
% 2013-08-15 - v2
% 2014-08-26 - v3
% 2017-06-12 - v4
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
%
% v4
% 2017-05-09
% Some general edits to how the A-field is calculated (nothing changed
% though). This was an update from the original A-field script. The
% a_field.m file from MagneticVectorPotentialFieldGenCode contains a
% deprecated version for computing the A-field at a minimum amount of
% points, e.g., a line, say.

% E-field equation:
% E = -nabla(phi)-d(A)/dt

% Coil located at [0,0,0] with the primary stimulation direction pointing in the y-direction
load coil.mat % coil; *.node *.data
% Magnetic permeability of free space
mu_0 = 4*pi*1E-7; % [N/A^2]  = [H/m]  = [J/(A^2m)]
mu_r = 1;

%% For Calculating Sensitivity of E-field in X and Y
% x = (-.1:.001:.1).';
% y = zeros(numel(x),1);
% z = ones(numel(x),1)*-.005;
% N = numel(x);
% 
% gps = ones(N,1);
% 
% x = cat(1,x,zeros(numel(x),1));
% y = cat(1,y,(-.1:.001:.1).');
% z = repmat(z,2,1);
% 
% gps = cat(1,gps,ones(N,1)*2);
% N = numel(x);

%% For Calculating Sensitivity of E-field in Z
% z = linspace(-.005,-.1,200)';
% N = numel(z);
% x = zeros(N,1);
% y = x;

%% For Calculating Sensitivity of E-field in pitch and yaw
% th = linspace(270-30,270+30,200)*pi/180;
% R = 0.05; % assume target is 1cm away... technically could use a variety of these.
% [y,z] = pol2cart(th,R);
% y = y';
% z = z';
% N = numel(y);
% x = zeros(N,1);
% gps = ones(N,1);
% 
% [tmpx,tmpz] = pol2cart(th,R);
% x = cat(1,x,tmpx');
% y = cat(1,y,zeros(N,1));
% z = cat(1,z,tmpz');
% gps = cat(1,gps,ones(N,1)*2);
% N = numel(x);

%% Create Plane
% [X,Y] = meshgrid(-140:3:140,-100:3:100);
% X = X./1e3;
% Y = Y./1e3;
% [xm,xn] = size(X);
% x = reshape(X,[],1);
% y = reshape(Y,[],1);
% z = ones(length(y),1)*-.01;

%% Create Box (for curl calc)
% [X,Y,Z] = meshgrid(-5:.1:5,-5:.1:5,-10:.1:0);
% X = X./1e3;
% Y = Y./1e3;
% Z = Z./1e3;
% [xm1,xm2,xm3] = size(X);
% N = numel(X);
% x = reshape(X,N,1);
% y = reshape(Y,N,1);
% z = reshape(Z,N,1);

%% Single point evaluation
% x = 0;
% y = 0;
% z = -.01;

n = length(coil.node);

Qp = mat2cell(coil.node,n,[1,1,1]);
Q = mat2cell(coil.data,n,[1,1,1]);

% What is this?
N = cat(1,ones(floor(N/200),1)*200,rem(N,200));

x = mat2cell(x,N,1);
y = mat2cell(y,N,1);
z = mat2cell(z,N,1);
N = num2cell(N);
Nlen = length(N);
A = cell(Nlen,1);

% This for-loop approach conserves RAM
for k = 1:Nlen
    % Magnetic vector potential calculation (A-field)
    M = mu_0*mu_r/4/pi./sqrt(...
        bsxfun(@minus,x{k},repmat(Qp{1}',N{k},1)).^2+...
        bsxfun(@minus,y{k},repmat(Qp{2}',N{k},1)).^2+...
        bsxfun(@minus,z{k},repmat(Qp{3}',N{k},1)).^2);
    
    A{k} = M*cell2mat(Q); % [V-s/m] [Wb/m] [H-A/m]
end

A = cell2mat(A);

a = vectormag(A);
% a = A(:,[1,2]);
x = cell2mat(x);
y = cell2mat(y);
z = cell2mat(z);


%% X and Y changes

pcx = bsxfun(@(a,ax) abs(a-ax)./ax,a,a(find(x==0,1,'first'),:)); % percent change

% 1D Plot
figure
plot(x(gps==1)*1e3,pcx(gps==1)*100,'k','linewidth',3)
hold on
plot(y(gps==2)*1e3,pcx(gps==2)*100,'r','linewidth',3)
hold off
xlim([-40,40])
grid on
ylabel('Percent Change')
xlabel('Coil Movement in X or Y Direction (mm)')
legend('Movement in X','Movement in Y')
saveFigureEps(gcf,6.4,3.2,'./EfieldSensitivity_XandY.png','png','-r200')

%% Z changes
% 1D Plot
figure
plot(z*1e3+5,pcx*100,'k','linewidth',3)
xlim([-40,0])
grid on
ylabel('Percent Change')
xlabel('Coil Movement in Z Direction (mm)')
legend('Movement in Z')
saveFigureEps(gcf,6.4,3.2,'./EfieldSensitivity_Z.png','png','-r200')

%% X Rotation Changes

[~,I] = min(z);
I = I(1);
pcxrot = bsxfun(@(a,ax) abs(a-ax)./ax,a,a(I,:)); % percent change

% 1D Plot
figure
plot(th*180/pi-270,pcxrot(gps==1)*100,'k','linewidth',3)
hold on
plot(th*180/pi-270,pcxrot(gps==2)*100,'r','linewidth',3)
hold off
xlim([-30,30])
grid on
ylabel('Percent Change')
xlabel('Coil Rotation About X or Y (deg)')
legend('Rotation about X','Rotation about Y')
title('E-field Change due to Pitch and Yaw')
% saveFigureEps(gcf,6.4,3.2,'./EfieldSensitivity_PitchYaw.png','png','-r200')