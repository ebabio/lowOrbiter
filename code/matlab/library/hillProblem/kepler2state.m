function x = kepler2state(elem)
% orbital parameters
a = elem(:,1);
e = elem(:,2);
eccA = elem(:,3);
lonAscNode = elem(:,4);
inclination = elem(:,5);
lonPeriapsis = elem(:,6);

% intermediate quantities
energy = -1./(2*a);
h = sqrt(a .*(1-e.^2));
trueA = 2 * atan( sqrt( (1+e)./(1-e) ) .* tan(eccA/2) );
r = a * (1 - e .* cos(eccA));
v = sqrt(2 *energy + 2./r );
gamma = real(acos( (h./r) ./ v ));


% positition and velocity
rVec = zeros(size(elem,1), 3);
vVec = zeros(size(elem,1), 3);

for i=1:size(elem,1)
    DCM = rotZXZ(lonAscNode(i), inclination(i), lonPeriapsis(i)+trueA(i));
    rVec(i,:) = DCM * [r(i); 0; 0];
    
    if(trueA(i) <0)
        gamma(i) = -1 * gamma(i);
    end
    DCM = rotZXZ(lonAscNode(i), inclination(i), lonPeriapsis(i)+trueA(i)+pi/2-gamma(i));
    vVec(i,:) = DCM * [v(i); 0,;0];
end

x = [rVec vVec];


