function [B, steps] = photobleach(traces,et)

T = length(traces);
steps = cell(T,2)
B = zeros(T,1);
for k = 1:T
	tic
	% Get individual trace and add time data
	dat = [traces(k,:); et:et:length(traces(k,:))*et];
	% Detect photobleaching steps
	rv = stepdetect(dat,k);
	% Indicate censoring as 1 if there are photobleaching steps
	if sum(rv{2}{2} < 0) > 0
		B(k) = 1
	else
		B(k) = 0;
	end
	% Keep all results in one cell
	steps(k,:) = rv(k,:);
	% Display progress
	if k == round(0.1*T)
		t = toc;
		fprintf('10%% (%d traces) finished! Time cost: %f minutes\n',k,t/60);
	elseif k == round(0.5*T)
		t = toc;
		fprintf('50%% (%d traces) finished! Time cost: %f minutes\n',k,t/60);
	elseif k == round(0.9*T)
		t = toc;
		fprintf('90%% (%d traces) finished! Time cost: %f minutes\n',k,t/60);
	end
end