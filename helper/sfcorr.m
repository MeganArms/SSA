function s = sfcorr(x,f,mu)
ko = koff(x,f,mu);
delta = 1;
dt = [diff(x);2];
s = zeros(length(x),1);
% ko(ko>1) = 0.4999;
for k = 1:length(x)
    if k == 1
        s(k) = 1;
    elseif ~isinf(ko(k))
        Dt = sum(dt(k-delta:k-1));
        s(k) = sum([s(k-delta),-s(k-delta)*ko(k-delta)*Dt],'omitnan');
        if s(k) < 0
            s(k) = 0;
            delta = delta + 1;
        else
            delta = 1;
        end
    elseif isinf(ko(k))
        delta = delta + 1;
    end
end
% s = sfcorr(x,f,mu)
% dt = diff(x); df = diff(f);
% s = [f(1); df - dt/mu];

% dt = diff(x);
% s = zeros(length(x),1);
% for k = 2:length(dt)
%     s(k) = f(k) + f(k-1)*dt(k)/mu;
% end
% s(1) = 1;