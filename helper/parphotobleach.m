function [B, steps] = parphotobleach(traces,et)
% Time is output in indices, must mulitply by lag time to get actual pb times
T = size(traces,1);
steps = cell(T,2);
B = zeros(T,1);
tic;
parfor q = 1:T
    % Get individual trace and add time data
	dat = [traces(q,:); et:et:length(traces(q,:))*et];
	% Detect photobleaching steps
	rv = stepdetect(dat,0,et); 
	% Indicate censoring as 1 if there are photobleaching steps
	if ~isempty(rv) || sum(rv{2}{2} < 0) > 1
		B(q) = 1;
	else
		B(q) = 0;
	end
	% Keep all results in one cell
    R = rv(1,:);
	steps(q,:) = R;
% 	% Display progress
% 	if q == round(0.1*T)
% 		t = toc;
% 		fprintf('10%% (%d traces) finished! Time cost: %f minutes\n',q,t/60);
% 	elseif q == round(0.5*T)
% 		t = toc;
% 		fprintf('50%% (%d traces) finished! Time cost: %f minutes\n',q,t/60);
% 	elseif q == round(0.9*T)
% 		t = toc;
% 		fprintf('90%% (%d traces) finished! Time cost: %f minutes\n',q,t/60);
% 	end
end