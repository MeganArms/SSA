function F = assignframes(s,r,t1,tau,Nframes)

th = 0.6;

f = s:tau:s+r;
if length(f) > 1 
    F = ceil((f-t1)/tau);
    if abs(s+r-f(end)) >= th*tau
        F = [F,F(end)+1];
%     elseif abs(s-t1) < th*tau
%         F(1) = [];
    end    
elseif length(f) == 1 && abs(s+r-f) < th*tau
    F = [];
elseif length(f) == 1 && abs(s+r-f) >= th*tau
    F = ceil((f-t1)/tau);
end

F = F(F>0 & F<=Nframes); % cannot assign frames beyond the duration of the video
end