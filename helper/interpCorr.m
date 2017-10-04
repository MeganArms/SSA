function [xmax, ymax] = interpCorr(C,n)

	a = max(max(C));
	i0 = find(c==a,1);
	[x0,y0] = ind2sub([n n],i0);
	b = 0.5*(C(x0+1,y0) - C(x0-1,y0));
	c = 0.5*(C(x0,y0+1) - C(x0,y0-1));
	d = 0.5*(C(x0+1,y0) + C(x0-1,y0)) - C(x0,y0);
	e = 0.5*(C(x0,y0+1) + C(x0,y0-1)) - C(x0,y0);

	xmax = b/(2*d); ymax = c/(2*e);
end