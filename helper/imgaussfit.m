function f = imgaussfit(img)

% Fit PSF
R = size(img,1);
% x = 1:R; y = 1:R;
[x,y] = meshgrid(1:R,1:R);
X = x(:); Y = y(:);
Z = img(:);


% PSF model
ft = fittype('A*exp((-(x-x0)^2-(y-y0)^2)/(2*sigma^2))+z0',...
    'independent',{'x','y'},'dependent','z');

% fit options
opts = fitoptions(ft);
opts.Display = 'off';
z_0 = 0;
% z_0 = min(min(img));
% if min(min(img - z_0)) < 0
%     error('Background setting too high');
% end
A_0 = img(R,R);
x_0 = R/2;
y_0 = R/2;
sigma_0 = 0.25*R;
opts.StartPoint = [A_0 sigma_0 x_0 y_0 z_0];
opts.Lower = [0 0.05*R 0.25*R 0.25*R z_0];
opts.Upper = [Inf 0.5*R 0.75*R 0.75*R z_0+1000];
% for i=1:5
%     if opts.Lower(i)>opts.Upper(i)
%         f=1;
%     end
% end

[f, ~] = fit([X,Y],Z,ft,opts);
