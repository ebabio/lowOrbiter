function elem = state2kepler(x)

% position and velocity 
rVec = x(:,1:3);
r = sqrt(rVec(:,1).^2 + rVec(:,2).^2 + rVec(:,3).^2);
rHat = rVec ./ (r * ones(1,3));
vVec = x(:,4:6);
v = sqrt(vVec(:,1).^2 + vVec(:,2).^2 + vVec(:,3).^2);

% angular momentum
hVec = cross(rVec, vVec, 2);
h = sqrt(hVec(:,1).^2 + hVec(:,2).^2 + hVec(:,3).^2);
hHat = hVec ./ (h * ones(1,3));
thetaHat = cross(hHat, rHat, 2);

% elements
a = -.5./(.5 * v.^2 - 1./r);
e = sqrt(1 - h.^2 ./ a);
inclination = acos(hHat(:,3));
lonAscNode = atan2(hHat(:,1),-hHat(:,2));

trueA = acos( (a.*(1-e.^2) ./ r - 1)./e );
ascending = -1 + 2 * ((rVec(:,1).*vVec(:,1) + rVec(:,2).*vVec(:,2) + rVec(:,3).*vVec(:,3)) > 0);
trueA = ascending .* trueA;
thetaA = atan2(rHat(:,3), thetaHat(:,3));
eccA =  ascending .* acos( (1 - r./a) ./ e );
omega = wrapToPi(thetaA - trueA);

elem = [a e eccA lonAscNode inclination omega];
