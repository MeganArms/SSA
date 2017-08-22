% Visualize overlap between long lifetimes and high frequency spots
inds = 1:numel(countsIm);
freqPlot = inds(countsIm ~= 0);
[x,y] = ind2sub([512 512],freqPlot);
r = iqr(ltIm(ltIm~=0));
lim75 = r/2+median(ltIm(ltIm~=0));
top75 = find(ltIm > lim75);
[p,q] = ind2sub([512 512],top75);
figure,plot(p,q,'k.')
hold on, plot(x,y,'b^'), hold off
legend('Top 75th Percentile of Lifetimes','Event Frequency','Location','best')
h = gca;
h.YLim = [0 512];
h.XLim = [0 512];
xlabel('X','FontSize',14),ylabel('Y','FontSize',14),title('overlap top75% w frequency')
clear h