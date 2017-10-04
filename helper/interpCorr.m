function [xmax, ymax] = interpCorr(C,n)

	a = max(max(C));
	i0 = find(C==a);
	[x0,y0] = ind2sub([n n],i0);
    if length(i0)>1
        r0 = [n/2+0.5,n/2+0.5];
        dx0 = abs(x0-r0); dy0 = abs(y0-r0);
        x0 = x0(dx0==min(dx0)); y0 = y0(dy0==min(dy0));
    end
	b = 0.5*(C(x0+1,y0) - C(x0-1,y0));
	c = 0.5*(C(x0,y0+1) - C(x0,y0-1));
	d = 0.5*(C(x0+1,y0) + C(x0-1,y0)) - C(x0,y0);
	e = 0.5*(C(x0,y0+1) + C(x0,y0-1)) - C(x0,y0);

	xmax = b/(2*d); ymax = c/(2*e);
end