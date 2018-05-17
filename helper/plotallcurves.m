function plotallcurves(x,Y)

% c = [0 0 0;1 0 0;0 0 1];
% figure,semilogy(x{1},Y{1},'s','MarkerFaceColor',c(1,:),'MarkerEdgeColor',c(1,:))
figure,semilogy(x{1},Y{1},'s')
hold on
for k = 2:length(Y)
    % plot(x{k},Y{k},'s','MarkerFaceColor',c(k,:),'MarkerEdgeColor',c(k,:))
    plot(x{k},Y{k},'s')
end

% figure,loglog(x,Y(1,:),'s-')
% hold on
% for k = 2:size((Y),1)
%     plot(x,Y(k,:),'s-')
% end

% M = max(cellfun(@max,Y));
% edges = linspace(0,M,100);
% figure
% for k = 1:length(Y)
%     c = histcounts(Y{k},edges);
%     % semilogy(edges(2:end),c,'s-'),hold on
%     semilogy(x{k}*edges(2:end),c,'s-'),hold on
% end

% legend('32 ms','64 ms','80 ms','100 ms','120 ms','160 ms','Location','best')
% legend('64 ms','80 ms','160 ms','Location','best')
xlabel('Connected residence time, \tau_c (s)'),ylabel('S(\tau_c)')
% xlabel('Light dosage, E (mJ)'),ylabel('S(E)')
h = gca;
h.FontSize = 14; h.FontWeight = 'bold';
clear h