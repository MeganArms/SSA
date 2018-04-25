% Find the lifetimes of the individual events within multiple events

r = zeros(length(spotEvents)*max(eventfreq),1);
l = 1;
for k = 1:length(spotEvents)
    lts = regionprops(~isnan(spotEvents(k).trajectory),'Area');
    ltsvec = [lts.Area]';
    if ~isnan(spotEvents(k).trajectory(1))
        ltsvec(1) = [];
    end
    censtmp = zeros(length(ltsvec),1);
    if spotEvents(k).trajectory(end) == Nframes
        censtmp(end) = 1;
    end
    r(l:length(ltsvec)+l-1) = ltsvec;
    c(l:length(ltsvec)+l-1) = censtmp;
    l = l + length(lts);
end
r = r(r>0);
c = c(r>0);