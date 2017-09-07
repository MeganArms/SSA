% Exponential Tests - need to run visualizeEvents to get "counts"
% lifetimefreqC = cell(length(counts),1); l = 1;
% for k = 1:length(counts)
%     lifetimestotake = counts(k);
%     if lifetimestotake ~= 0
%         lifetimefreqC{k} = indLengths_allTraj(l:l+lifetimestotake-1);
%         l = l+lifetimestotake;
%     end
% end
% 
% testlife = lifetimefreqC(cellfun(@(x) length(x)>3,lifetimefreqC));
testlife = cellfun(@(x) 2*(x-3),testlife,'UniformOutput',0);
lillian = cellfun(@(x) lillietest(x,'Distr','exp'),testlife);
ad = cellfun(@(x) adtest(x,'Distr','exp'),testlife);
agreement = sum(lillian == ad)/length(testlife);
Lpct = (length(testlife)-sum(lillian))/length(testlife);
ADpct = (length(testlife)-sum(ad))/length(testlife);

figure
for  k = 1:length(testlife)
    if ad(k) == 1
        text1 = 'A-D: Not Exponential';
    else
        text1 = 'A-D: Exponential';
    end
    if lillian(k) == 1
        text2 = 'Lilliefors: Not Exponential';
    else
        text2 = 'Lilliefors: Exponential';
    end
    histogram(testlife{k}),title([text1,', ',text2,', Spot No.',num2str(k),'/',num2str(length(testlife))])
    pause(2)
end


