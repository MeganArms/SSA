% Thompson error - this could be converted to consider empircal s
a = 160;
s = 561/2/1.49;
gainMultiplier = 150;
b = ib/gainMultiplier;
% Get total photons collected for each particle as the sum of integrated
% intensities over the entire trajectory, each divided by the gain
% multiplier
Blong = B(cellfun(@length,B)>2); counts = cellfun(@sum,Blong);
N = counts/gainMultiplier;

delx = (s^2 + a^2/12)./N + (8*pi*s^4*b^2)/a^2./N.^2;
rmsx = sqrt(delx);

Nt = 4*sqrt(pi)*s^3*b^2/(a*(s^2 + a^2/12));