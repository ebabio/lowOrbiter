%% Plot phase portrait for the eccentricity vector

close all
clear

% Add path
path = [pwd,'/library']; %change to: added_path = '/path' for your required path
addpath(path);

clc

%% Function and domain definition

% Define ode
unstable=1;
alpha = 3;
if(unstable)
    alpha= 1/alpha;
end
odeEcc =@(~,x) (odeEccentricity([], x, alpha));
stabOmega = 1/2 * acos(alpha);

eDomain = linspace(0,.5,6);
omegaDomain = linspace(-pi, pi, 30+1);
eT = [.5 .5 .5 .5 .1 .1 .1 .1];
omegaT = [pi/2+stabOmega+0.01    pi/2+stabOmega-0.01  -pi/2+stabOmega+0.01   -pi/2+stabOmega-0.01 ...
          -pi/2-stabOmega+0.1     -pi/2-stabOmega-0.1   pi/2-stabOmega+0.1     pi/2-stabOmega-0.1];
tspan = 10;

%% Plot Vector field

% Evaluate
[e, omega] = meshgrid(eDomain(2:end), omegaDomain(2:end));
ecc = zeros([size(e) 2]);
eccDot = zeros([size(e) 2]);
for i = 1:size(e,1)
    for j = 1:size(e,2)
        ecc(i,j,:) = e(i,j) * [cos(omega(i,j)); sin(omega(i,j))];
        xDot = odeEcc(0,[e(i,j); omega(i,j)]);
        eccDot(i,j,:) =.05 * (xDot(1) * [cos(omega(i,j)); sin(omega(i,j))] + xDot(2) * [-sin(omega(i,j)); cos(omega(i,j))]);
    end
end

% Plot
figure(1)
title(['Phase portrait for the eccentricity vector for \alpha=' num2str(alpha)])
set(gcf,'units','inches','position',[0 0 8.5 6]);
set(gcf,'units','inches','PaperPosition',[0 0 8.5 6]);

xlabel('e cos \omega')
ylabel('e sin \omega')
hold on

q = quiver(ecc(:,:,1), ecc(:,:,2), eccDot(:,:,1), eccDot(:,:,2));
q.AutoScale = 'off';
axis equal

%% Plot trajectories
% Define integration termination conditions4
event = @(~,x) (integrationEvents([], x, max(eDomain)));

odeOptions = odeset('Reltol', 1e-6, 'AbsTol',1e-9, 'Events',event);
% Integrate trajectories
for i = 1:numel(eT)
    z0 = [eT(i); omegaT(i)];
    [t,z] = ode45(odeEcc,[0 tspan],z0, odeOptions);
    Z{i} = z';
    plot(Z{i}(1,:).*cos(Z{i}(2,:)), Z{i}(1,:).*sin(Z{i}(2,:)))
end

saveas(gcf, 'results/EccentricityPhasePortrait.png')

