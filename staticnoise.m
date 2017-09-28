function [m, b, ib, imstd] = staticnoise(filename)

iminf = imfinfo(filename); nFrames = length(iminf);
M = iminf.Width; diffBG = zeros(M,M,nFrames-1);
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
    diffBG(:,:,k) = im2 - im1;
    if k == nFrames - 1
        ipn(k+1) = std(im1(:)-im2(:))/sqrt(2);
        ips(k+1) = mean(im2(:));
    end
end
stnoise = fit(ips,ipn,'poly1'); m = stnoise.p1; b = stnoise.p2;
imstd = mean(mean(std(diffBG,0,3)));

% Determine the background level
% v may be poorly affected by the number of bins, but as long as it's not
% too few, should be okay. Consider replacing with the rank-sum method in
% the future.
[counts,edges] = histcounts(ips,100,'normalization','pdf'); 
x = edges(1:end-1);
f = fit(x',counts','gauss1');
ib = f.b1;