% Normal tests on brightness
brightnessD = cell(length(counts),1); l = 1;
for k = 1:length(counts)
    lifetimestotake = counts(k);
    if lifetimestotake ~= 0
        brightnessD{k} = brightnesses_allTraj(l:l+lifetimestotake-1);
        l = l+lifetimestotake;
    end
end

testbright = brightnessD(cellfun(@(x) length(x)>3,brightnessD));
lillian = cellfun(@lillietest,testbright);
ad = cellfun(@adtest,testbright);
agreement = sum(lillian == ad)/length(testbright);
Lpct = (length(testbright)-sum(lillian))/length(testbright);
ADpct = (length(testbright)-sum(ad))/length(testbright);

figure
for  k = 1:length(testbright)
    if ad(k) == 1
        text1 = 'A-D: Not Normal';
    else
        text1 = 'A-D: Normal';
    end
    if lillian(k) == 1
        text2 = 'Lilliefors: Not Normal';
    else
        text2 = 'Lilliefors: Normal';
    end
    histogram(testbright{k},'NumBins',length(testbright{k})+1),title([text1,', ',text2,', Spot No.',num2str(k),'/',num2str(length(testbright))])
    pause(2)
end
