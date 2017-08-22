function bp = cold(subimg)

% COLSP (COLumn Density) is the function designed to work with COLFILT
% and a 'sliding' type 3x3 neighborhood. The purpose is to normalize the
% sum of the intensities of the region to the size of the region to get the
% density. This can be performed on grayscale or binary images.

% Generalizable with padding of the lut's to the size of nhood and changing
% the "9"s to be m*n product of size of nhood.

% subimg(subimg > 0) = 1;
[~,n] = size(subimg);
% dispm = floor(m/3):floor(2*m/3);
% dispn = 1:floor(n/3);
% padm = heaviside(m - 9)*(m - 9);
% padn = heaviside(n - 1)*(n - 1);

lut = [1 1 1; 1 1 1; 1 1 1];
bp = zeros(1,n);

lutcurrent = reshape(lut,[9,1]);
lutmatrix = repmat(lutcurrent,1,n);
outs = sum(subimg.*lutmatrix,1);
bp(1,:) = outs/(9);

end
