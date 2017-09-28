function [IMG,CCD_IMG,F0,Ff] = imagesimulator3000(m,b,imstd,ib,snr,NA,lambda,ps,gs,sf,iter,varargin)
% % Define the integral function for later substitution
% syms x y I mux muy a ib
% % int_fx = ib*x - (2^(1/2)*I*pi^(1/2)*exp(-muy^2/(2*a^2))*exp(-y^2/(2*a^2))...
% %     *exp((muy*y)/a^2)*erf((2^(1/2)*(mux - x)*(1/a^2)^(1/2))/2))/...
% %     (2*(1/a^2)^(1/2));
% int_f = ib*x*y - (I*pi*erf(-(2^(1/2)*((muy*1i)/a^2 - (y*1i)/a^2))...
%     /(2*(-1/a^2)^(1/2)))*erf((2^(1/2)*(mux - x)*(1/a^2)^(1/2))/2)*1i)...
%     /(2*(1/a^2)^(1/2)*(-1/a^2)^(1/2));



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
    ip = I*exp(-((xg-mux).^2 + (yg-muy).^2)/(2*a^2)) + ib; % constant bg doesn't seem right
    % ip = f(xg,yg,I,ib,a,mux,muy);
    
    % Add noise to HRG
    IMG(:,:,k) = ip + (m*ip+b).*normrnd(0,1,size(ip));
    
    % Sum to obtain pixel-resolution grid
    CCD_IMG = reshape(mean(im2col(IMG,[ps^2/sf^2 ps^2/sf^2],'distinct'),1),[gs,gs]);
   
    % Integrate to obtain pixel-resolution grid
%     F = matlabFunction(subs(int_f));
%     x0 = 1:ps:gs*ps; xf = ps:ps:ps*gs; y0 = x0; yf = xf;
%     [x0g,y0g] = meshgrid(x0,y0); [xfg,yfg] = meshgrid(xf,yf);
%     CCD_IMG = F(xfg,yfg) - F(x0g,y0g);
    
    % No anonymous or symbolic function method - hold on. There's no noise
    % in this function.
%     x0 = 1:ps:gs*ps; xf = ps:ps:ps*gs; y0 = x0; yf = xf;
%     [x0g,y0g] = meshgrid(x0,y0); [xfg,yfg] = meshgrid(xf,yf);
%     F0 = ib.*x0g.*y0g - (I*pi*erf(-(2^(1/2)*((muy*1i)/a^2 - (y0g*1i)...
%         /a^2))/(2*(-1/a^2)^(1/2))).*erf((2^(1/2)*(mux - x0g)*...
%         (1/a^2)^(1/2))/2)*1i)/(2*(1/a^2)^(1/2)*(-1/a^2)^(1/2));
%     Ff = ib.*xfg.*yfg - (I*pi*erf(-(2^(1/2)*((muy*1i)/a^2 - (yfg*1i)...
%         /a^2))/(2*(-1/a^2)^(1/2))).*erf((2^(1/2)*(mux - xfg)*...
%         (1/a^2)^(1/2))/2)*1i)/(2*(1/a^2)^(1/2)*(-1/a^2)^(1/2));
%     CCD_IMG(:,:,k) = Ff - F0;
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
