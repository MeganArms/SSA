function [static_error, err_x, err_y] = imfit(IMG,mux,muy,sf,ps)

% Convert mux into "units" of pixels
muX = mux*sf/ps; muY= muy*sf/ps;

% Fit the simulated images
[xc,yc,~,~] = radialcenter_stk(IMG);

% Find difference in position
err_x = (muX - xc)*ps;
err_y = (muY - yc)*ps;

% Find static error
static_error = std(err_x) + std(err_y);