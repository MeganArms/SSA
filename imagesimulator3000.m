function IMG = imagesimulator3000(m,b,ib,B,S,NA,lambda,ps,gs,sf,iter,varargin)
% Generate high-resolution grid (HRG). 
xx = sf/ps:sf/ps:ps*gs/sf; [xg,yg] = meshgrid(xx,xx);

% Determine the apparent particle radius
a = lambda/2/NA/sf; 

% Determine the average signal
Blong = B(cellfun(@length,B)>2); avgB = mean(cellfun(@mean,Blong));
Slong = S(cellfun(@length,S)>2); avgS = mean(cellfun(@mean,Slong));
avgPeakSignal = avgB*0.341/avgS;

% Define true location of particle. Default is centered.
if isempty(varargin)
    mux = ceil(ps/sf*gs/2); muy = ceil(ps/sf*gs/2);
else
    mux = varargin{1}; muy = varargin{2};
end

% Run simulation
IMG = zeros(gs,gs,iter); 
% f = @(x,y,A,C,a,mux,muy) A*exp(-((x-mux).^2+(y-muy).^2)/(2*a^2)) + C;
for k = 1:iter
    % Determine the peak signal intensity with mean signal intensity from
    % the image file of interest
    I = poissrnd(avgPeakSignal);
    % Simulate particle on an HRG
    % constant bg doesn't seem right. but this is the noiseless image.
    ip = I*exp(-((xg-mux).^2 + (yg-muy).^2)/(2*a^2)) + ib; 
    % Sum to obtain pixel-resolution grid
    im = reshape(mean(im2col(ip,[ps^2/sf^2 ps^2/sf^2],'distinct'),1),[gs,gs]);
    % Add noise to pixel-resolution grid
    IMG(:,:,k) = im + (m*im+b).*normrnd(0,1,size(im));
end