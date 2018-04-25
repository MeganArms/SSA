function coords = viewblinks(spotEvents,objs_link,Nframes)

coords = nan(Nframes*length(spotEvents),2);
l = 1;
for k = 1:length(spotEvents)
    molecnums = spotEvents(k).trajectory;
    molecnums = molecnums(~isnan(molecnums));
    inds = l:l+length(molecnums)-1;
    coords(inds,:) = objs_link(1:2,molecnums)';
    l = l + length(molecnums);
end
figure,plot(coords(:,1),coords(:,2),'k.')