function k = koff(x,f,mu)

A = diff(f);
A = [A; A(end)];
dT = diff(x);
dT = [dT; dT(end)];

k = -(A./dT + f/mu)./f;
