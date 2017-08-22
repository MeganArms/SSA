function z = landingrate(F,Nframes,lag,plotting)

% Landing rate provides the number of landing events per second. Landing
% events are events in which a bright spot appears for at least three
% frames.
% INPUT:
%   F - output from COLLATE with the frame numbers of molecule trajectories 
%   NFRAMES - number of frames in video
%   LAG - interval between the start of each frame
%   PLOTTING - set as 'on' to view the histogram of new events per frame,
%   labelled with the average landing rate.
% OUTPUT:
%   Z - landing events per second

longTrajs = F(cellfun(@length,F) > 2);
y = cellfun(@(x) x(1), longTrajs);

N = zeros(Nframes,1);
for t = 1:Nframes
    N(t) = sum(y == t);
end

% time = lag:lag:lag*Nframes;
% landingmat = [time',N];
% landingmat = landingmat(N ~= 0,:);

% f = fit(landingmat(:,1),landingmat(:,2),'poly1');
% f = fit(time(2:end)',N(2:end),'poly1');
% z = f.p2;
z = mean(N(2:end))/lag;

if strcmp(plotting,'on')
%     figure, plot(landingmat(:,1),landingmat(:,2),'ks')
%     figure, plot(time',N,'ks')
%     hold on, plot(f), hold off
%     xlabel('Time, t (s)'), ylabel('N(t)')
%     h = gca; h.YLim = [0 max(N)]; h.Children(2).MarkerFaceColor = [0 0 0];
    figure, histogram(N(2:end));
    xlabel('Number of new events per frame, N'), ylabel('f(N)')
    title({'Average landing rate: ',num2str(z),' s^{-1}'})
end

end