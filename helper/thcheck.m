function b = thcheck(x,th)
winsize = length(x);
% th = 0.25*(max(x) - min(x));
lut = zeros(1,winsize);

lut1 = lut;
lut1(1:floor(winsize/2)) = 1;
lut2 = lut;
lut2(ceil(winsize/2):end) = 1;

if (sum((lut1.*x)/sum(lut1)) - median(x)) >= th && (median(x) - sum((lut2.*x)/sum(lut2))) >= th
	b = 1;
else
	b = 0;
end

end