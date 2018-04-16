%% Test 1:
% Obtaining Lyapunov family of orbits

%% Setup workspace
clear
% Add path
path1 = [pwd,'/library'];
path2 = [pwd,'/library/hillProblem'];
addpath(path1, path2);

Constants

clc

%% Define systems
JupiterSat = HillSystem(Jupiter, Ganymede);

%% Define initial conditions
a = 0.1668;
e = 0.0;
eccA = 0;
lonAscNode = 0;
inclination = deg2rad(87);
lonPeriapsis = deg2rad(-39.441); %-47.8847 for Case 1, -39.441 for Case 3

tf = 18*2*pi;
globalize = 1;

%% Integrate one orbit and plot it
% find initial state
elem0 = [a e eccA lonAscNode inclination lonPeriapsis];
x0 = kepler2state(elem0);
x0 = x0' - [0; 0; 0; cross([0;0;1], x0(1:3)')];

% integrate: standard
if(~globalize)
    orbit = HillOrbit(x0, JupiterSat);
    orbit.setStoppingCondition(@(t,x) norm(x(1:3))<(JupiterSat.centralRadius), 0);
    orbit.integrateX0([0, tf]);
    dC = orbit.jacobiAnalysis();
end

% integrate: global manifolds
if(globalize)
    orbit = HillOrbit(x0, JupiterSat);
    orbit.setStoppingCondition(@(t,x) norm(x(1:3))<(JupiterSat.centralRadius), 0);
    orbit1 = orbit.copy();
    orbit1.integrateX0([0, -tf]);
    orbit2 = orbit.copy();
    orbit2.integrateX0([0, +tf]);
    dC = orbit1.jacobiAnalysis();
    % merge orbits
    orbit.t = [fliplr(orbit1.t) orbit2.t] - orbit1.tf; %1st orbit is progated backwards
    orbit.x = [fliplr(orbit1.x) orbit2.x];
    orbit.tf = max(orbit.t);
    orbit.x0 = orbit.x(:,1);
    orbit.xf = orbit.x(:,end);
end

% inertial analysis
xInertial = orbit.inertialMotion(orbit.t, orbit.x);
keplerianElements = state2kepler(xInertial');

%% Display results
[timeVal, units] = reasonableTime(orbit.tf*orbit.System.t);

figure(1)
clf reset
set(gcf,'units','inches','position',[0 0 7 5]);
set(gcf,'units','inches','PaperPosition',[0 0 7 5]);
JupiterSat.plotSystem3d();
orbit.plot();
view(30,45)
saveas(gcf, 'results/motionRotatingFrame.png')

figure(2)
set(gcf,'units','inches','position',[0 0 7 5]);
set(gcf,'units','inches','PaperPosition',[0 0 7 5]);
clf reset
hold on
JupiterSat.plotSystem3d();
plot3(xInertial(1,:), xInertial(2,:), xInertial(3,:))
title(['Motion in the non-rotating frame for ' num2str(orbit.tf) units])
title(['Motion in the non-rotating frame for ' num2str(timeVal) units])
xlabel('x [non-dimensional units]')
ylabel('y [non-dimensional units]')
zlabel('z [non-dimensional units]')
view(30,15)
saveas(gcf, 'results/motionFixedFrame.png')

figure(3)
set(gcf,'units','inches','position',[0 0 7 5]);
set(gcf,'units','inches','PaperPosition',[0 0 7 5]);
clf reset
subplot(2,2,1)
plot(orbit.t, keplerianElements(:,1)')
title('semi-major axis')
xlabel('t[nd units]')
ylabel('a[nd units]')

subplot(2,2,2)
plot(orbit.t, keplerianElements(:,2)')
title('eccentricity')
xlabel('t[nd units]')
ylabel('e')

subplot(2,2,3)
plot(orbit.t, rad2deg(keplerianElements(:,5)'))
title('inclination')
xlabel('t[nd units]')
ylabel('i[deg]')

subplot(2,2,4)
plot(orbit.t, rad2deg(keplerianElements(:,6)'))
title('lonPeriapsis')
xlabel('t[nd units]')
ylabel('\omega[deg]')
saveas(gcf, 'results/elements.png')

figure(4)
set(gcf,'units','inches','position',[0 0 7 5]);
set(gcf,'units','inches','PaperPosition',[0 0 7 5]);
clf reset
hold on
h = keplerianElements(:,2)'.* cos(keplerianElements(:,6)');
k = keplerianElements(:,2)'.* sin(keplerianElements(:,6)');
plot(h, k)
title('motion of the eccentricity vector')
xlabel('e cos \omega')
ylabel('e sin \omega')
axis equal
saveas(gcf, 'results/eccentricity.png')

if(globalize)
    figure(5)
    set(gcf,'units','inches','position',[0 0 7 5]);
    set(gcf,'units','inches','PaperPosition',[0 0 7 5]);
    clf reset
    hold on
    rp = (keplerianElements(:,1)'.* (1-keplerianElements(:,2)')-JupiterSat.centralRadius) * JupiterSat.l;
    ra = (keplerianElements(:,1)'.* (1+keplerianElements(:,2)')-JupiterSat.centralRadius) * JupiterSat.l;
    days = orbit.t*(JupiterSat.t/86400);
    plot(days, rp)
    plot(days, ra)
    title('Altitude profile for JUICE')
    xlabel('t [days]')
    ylabel('altitude [km]')
    legend('periapsis','apoapsis')
    saveas(gcf, 'results/missionProfileDim.png')
    
    figure(6)
    set(gcf,'units','inches','position',[0 0 7 5]);
    set(gcf,'units','inches','PaperPosition',[0 0 7 5]);
    clf reset
    hold on
    rpNd = (keplerianElements(:,1)'.* (1-keplerianElements(:,2)'));
    raNd = (keplerianElements(:,1)'.* (1+keplerianElements(:,2)'));
    days = orbit.t*(JupiterSat.t/86400);
    plot(orbit.t, rpNd)
    plot(orbit.t, raNd)
    title('Altitude profile for JUICE')
    xlabel('t [nd units]')
    ylabel('distance [nd units]')
    legend('periapsis','apoapsis')
    saveas(gcf, 'results/missionProfileNd.png')
end

display(['nondimensional time to impact: ' num2str(orbit.t(end))])
