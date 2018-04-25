function viewtraces(traces,Nframes,et)
figure; 
t = et:et:Nframes*et;
for i = 1:size(traces,1) % i = 1298:size(traces,1)
    plot(t,traces(i,:));title(['Trace ',num2str(i)]);
    set(gca,'YLim',[900 5000])
    pause(1)
end