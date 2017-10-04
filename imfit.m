function [static_error, err_x, err_y, xp, yp] = imfit(IMG,mux,muy,sf,ps)

% Convert mux into "units" of pixels
muX = mux*sf/ps; muY= muy*sf/ps;

% Fit the simulated images
[xc,yc,~,~] = radialcenter_stk(IMG);

% Find difference in position
err_x = (muX - xc)*ps;
err_y = (muY - yc)*ps;

% Find static error
static_error = std(err_x) + std(err_y);

% Calculate "mass" (intensity) in each neighborhood
% subtract background, estimated from mean of non-neighborhood px.
nhood = getnhood(strel('disk', floor(size(IMG,1)/2),0)); 
meanbkg = zeros(1, length(IMG));
savemass = zeros(1, length(IMG));
for k = 1:length(IMG)
    tempreg = IMG(:,:,k);
    cropreg = tempreg(nhood);
    bkgreg  = tempreg(~nhood);
    meanbkg(k) = mean(bkgreg(:));
    savemass(k) = sum(cropreg(:)) - sum(nhood(:))*meanbkg(k);
end

% Find precision
xp = std(xc)/sum(savemass);
yp = std(yc)/sum(savemass);