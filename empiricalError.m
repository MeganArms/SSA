% 1
% lengths = cellfun(@length,C);
% empties = cellfun(@isempty,C);
% Clong = C(lengths>2&~empties);
% Blong = B(lengths>2&~empties);
% semx = cellfun(@(x)std(x(1,:)),Clong)./cellfun(@sum,Blong);
% semy = cellfun(@(x)std(x(2,:)),Clong)./cellfun(@sum,Blong);
% loc_err = mean(semx,'omitnan') + mean(semy,'omitnan');
% 
% lengths = cellfun(@(x)size(x,1),Ccorr);
% Clong = Ccorr(lengths>2&~empties);
% semx_corr = cellfun(@(x)std(x(:,1)),Clong)./cellfun(@sum,Blong);
% semy_corr = cellfun(@(x)std(x(:,2)),Clong)./cellfun(@sum,Blong);
% loc_err_corr = mean(semx_corr,'omitnan') + mean(semy_corr,'omitnan');

% 2
% empties = cellfun(@isempty,C); Clong = C(~empties);
% Blong = B(~empties);
% semx = cellfun(@(x)std(x(1,:)),Clong);%./cellfun(@sum,Blong);
% semy = cellfun(@(x)std(x(2,:)),Clong);%./cellfun(@sum,Blong);
% loc_err = mean(semx,'omitnan') + mean(semy,'omitnan');
% 
% empties = cellfun(@isempty,Ccorr); Clong = Ccorr(~empties); Blong = B(~empties);
% semx_corr = cellfun(@(x)std(x(:,1)),Clong);%./cellfun(@sum,Blong);
% semy_corr = cellfun(@(x)std(x(:,2)),Clong);%./cellfun(@sum,Blong);
% loc_err_corr = mean(semx_corr,'omitnan') + mean(semy_corr,'omitnan');

% 3
empties = cellfun(@isempty,C); Clong = C(~empties);
% Blong = B(~empties);
% semx = cellfun(@(x)sqrt(sum((x(1,:)-mean(x(1,:))).^2)),Clong)...
%     ./cellfun(@(x)sqrt(sum(x)),Blong);
% semy = cellfun(@(x)sqrt(sum((x(2,:)-mean(x(2,:))).^2)),Clong)...
%     ./cellfun(@(x)sqrt(sum(x)),Blong);
% loc_err = mean(semx,'omitnan') + mean(semy,'omitnan');

empties = cellfun(@isempty,Ccorr); Clongc = Ccorr(~empties); Blong = B(~empties);
% semx = cellfun(@(x)sqrt(sum((x(:,1)-mean(x(:,1))).^2)),Clong)...
%     ./cellfun(@(x)sqrt(sum(x)),Blong);
% semy = cellfun(@(x)sqrt(sum((x(:,2)-mean(x(:,2))).^2)),Clong)...
%     ./cellfun(@(x)sqrt(sum(x)),Blong);
% loc_err_corr = mean(semx_corr,'omitnan') + mean(semy_corr,'omitnan');

% 4
ssex = cell(length(Clong),1); ssey = ssex; ssex_c = ssex; ssey_c = ssey;
meansx = zeros(length(Clong),1); meansy = meansx; meansxc = meansx; meansyc = meansy;
for k = 1:length(Clong)
    currtraj = Clong{k};
    meansx(k) = mean(currtraj(1,:)); meansy(k) = mean(currtraj(2,:));
    ssex{k} = (currtraj(1,:)-mean(currtraj(1,:))).^2;
    ssey{k} = (currtraj(2,:)-mean(currtraj(2,:))).^2;
    currtraj = Clongc{k};
    meansxc(k) = mean(currtraj(:,1)); meansyc(k) = mean(currtraj(:,2));
    ssex_c{k} = (currtraj(:,1)-mean(currtraj(:,1))).^2;
    ssey_c{k} = (currtraj(:,2)-mean(currtraj(:,2))).^2;
end
        
