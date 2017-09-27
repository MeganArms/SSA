function [IMG,CCD_IMG] = imagesimulator3000(filename,snr,NA,lambda,ps,gs,iter,varargin)
% % Define the integral function for later substitution
% syms x y I mux muy a ib
% % int_fx = ib*x - (2^(1/2)*I*pi^(1/2)*exp(-muy^2/(2*a^2))*exp(-y^2/(2*a^2))...
% %     *exp((muy*y)/a^2)*erf((2^(1/2)*(mux - x)*(1/a^2)^(1/2))/2))/...
% %     (2*(1/a^2)^(1/2));
% int_f = ib*x*y - (I*pi*erf(-(2^(1/2)*((muy*1i)/a^2 - (y*1i)/a^2))...
%     /(2*(-1/a^2)^(1/2)))*erf((2^(1/2)*(mux - x)*(1/a^2)^(1/2))/2)*1i)...
%     /(2*(1/a^2)^(1/2)*(-1/a^2)^(1/2));

% Static noise function - make this it's own function.
% This is gonna be so darn slow. But it only has to happen once, not for
% every simulation iteration, or even every position on the HRG. Only run
% the first time for the first position.
iminf = imfinfo(filename); nFrames = length(iminf);
ipn = zeros(nFrames,1); ips = ipn;
for k = 1:nFrames-1
    if k == 1
        im1 = double(imread(filename,k));
        im2 = double(imread(filename,k+1));
        sumBG = im1;
    else
        im1 = im2;
        im2 = double(imread(filename,k+1));
        sumBG = sumBG + im2;
    end
    ipn(k) = std(im2(:)-im1(:))/sqrt(2);
    ips(k) = mean(im1(:));
    if k == nFrames - 1
        ipn(k+1) = std(im1(:)-im2(:))/sqrt(2);
        ips(k+1) = mean(im2(:));
    end
end
stnoise = fit(ips,ipn,'poly1'); m = stnoise.p1; b = stnoise.p2;

% Determine the background level
% v may be poorly affected by the number of bins, but as long as it's not
% too few, should be okay. Consider replacing with the rank-sum method in
% the future.
[counts,edges] = histcounts(ips,100,'normalization','pdf'); 
x = edges(1:end-1);
f = fit(x,counts,'gauss1');
ib = f.b1;

% Simulate particle on high-resolution grid (HRG). A pixel is 160 nm, so
% the grid will be 1/160 pixels in resolution. 
% Maybe worth considering the pros and cons of using units of pixels or
% nanometers.
% ps = 160; gs = 3;
xx = 0:1/ps:ps*gs; 
[xg,yg] = meshgrid(xx,xx);

% Determine the apparent particle radius
a = 2*pi*NA/lambda; % NA and lambda are user inputs

% The subpixel location will be take values ranging between from the center
% of the HRG minus 1.5 pixels to the center plus 1.5 pixels. Start with the
% center position, but in future it will be a user input. Grid goes from
% 0-gridsize (gs), so the center is at ceil(gs/2)
if isempty(varargin)
    mux = ceil(ps*gs/2); muy = ceil(ps*gs/2);
else
    mux = varargin{1}; muy = varargin{2};
end

% Run simulation
IMG = zeros(gs*ps,gs*ps,iter);
CCD_IMG = zeros(gs,gs,iter);

for k = 1:iter
    % Determine the peak signal intensity. SNR is user input.
    I = poissrnd(snr*var(sumBG/nFrames));

    % Simulate particle on an HRG
    ip = I*exp(-((xg-mux)^2 + (yg-muy)^2)/(2*a^2)) + ib;

    % Add noise to HRG
    IMG = ip + (m*ip+b)*normrnd(0,1);

    % Integrate to obtain pixel-resolution grid
%     F = matlabFunction(subs(int_f));
%     x0 = 1:ps:gs*ps; xf = ps:ps:ps*gs; y0 = x0; yf = xf;
%     [x0g,y0g] = meshgrid(x0,y0); [xfg,yfg] = meshgrid(xf,yf);
%     CCD_IMG = F(xfg,yfg) - F(x0g,y0g);
    
    % No anonymous or symbolic function method
    F0 = ib*x0g*y0g - (I*pi*erf(-(2^(1/2)*((muy*1i)/a^2 - (y0g*1i)...
        /a^2))/(2*(-1/a^2)^(1/2)))*erf((2^(1/2)*(mux - x0g)*...
        (1/a^2)^(1/2))/2)*1i)/(2*(1/a^2)^(1/2)*(-1/a^2)^(1/2));
    Ff = ib*xfg*yfg - (I*pi*erf(-(2^(1/2)*((muy*1i)/a^2 - (yfg*1i)...
        /a^2))/(2*(-1/a^2)^(1/2)))*erf((2^(1/2)*(mux - xfg)*...
        (1/a^2)^(1/2))/2)*1i)/(2*(1/a^2)^(1/2)*(-1/a^2)^(1/2));
    CCD_IMG = Ff - F0;
end



% Fx_subs = subs(Fx_definite,x,[x0,xf]);
% Fx = Fx_subs(2,:) - Fx_subs(1,:);
% intf
% 
% int_fx = matlabFunction(subs(int_fx));
% 
% int_f = matlabFunction(subs(int_f)); % use this to integrate to a low res grid
% 
% 
% B = im2col(IMG,[ps ps],'distinct');
% % Either need to rewrite the function to take columns or to take block
% % structs
