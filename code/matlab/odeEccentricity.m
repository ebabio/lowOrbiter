function xDot = odeEccentricity(~, x, alpha)

e = x(1);
omega = x(2);

eDot = sin(2*omega) * e;
omegaDot = cos(2*omega) + alpha;

xDot = [eDot; omegaDot];
end