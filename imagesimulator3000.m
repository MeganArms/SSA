function [IMG,CCD_IMG] = imagesimulator3000(m,b,imstd,ib,snr,NA,lambda,...
    ps,gs,sf,iter,varargin)
% Simulate particle on high-resolution grid (HRG). A pixel is 160 nm, so
% the grid will be 1/160 pixels in resolution. 
% Maybe worth considering the pros and cons of using units of pixels or
% nanometers.
% ps = 160; gs = 30;
xx = sf/ps:sf/ps:ps*gs/sf; 
[xg,yg] = meshgrid(xx,xx);

% Determine the apparent particle radius
a = lambda/2/NA/sf; % NA and lambda are user inputs

% The subpixel location will be take values ranging between from the center
% of the HRG minus 1.5 pixels to the center plus 1.5 pixels. Start with the
% center position, but in future it will be a user input. Grid goes from
% 0-gridsize (gs), so the center is at ceil(gs/2)
if isempty(varargin)
    mux = ceil(ps/sf*gs/2); muy = ceil(ps/sf*gs/2);
else
    mux = varargin{1}; muy = varargin{2};
end

% Run simulation
IMG = zeros(length(xx),length(xx),iter);
CCD_IMG = zeros(gs,gs,iter); 
% f = @(x,y,A,C,a,mux,muy) A*exp(-((x-mux).^2+(y-muy).^2)/(2*a^2)) + C;
for k = 1:iter
    % Determine the peak signal intensity. SNR is user input.
    I = poissrnd(snr*imstd);

    % Simulate particle on an HRG
    % constant bg doesn't seem right
    ip = I*exp(-((xg-mux).^2 + (yg-muy).^2)/(2*a^2)) + ib; 
    
    % Add noise to HRG
    IMG(:,:,k) = ip + (m*ip+b).*normrnd(0,1,size(ip));
    
    % Sum to obtain pixel-resolution grid
    CCD_IMG = reshape(mean(im2col(IMG,[ps^2/sf^2 ps^2/sf^2],'distinct'),...
        1),[gs,gs]);
end