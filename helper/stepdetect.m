%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Photobleaching Step Detection
%     Copyright (C) 2014
%     Venkatramanan Krishnamani, Rahul Chadda & ...
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   bleaching_step_detection(dat,fig_num)
% - modified from heaviside_step_detection_v1
%   to process bleaching trace (Jan 2012)
%
%   detect bleaching steps in intensity trace by fitting a sliding
%   Heaviside step funtion.
%
%   The output data matrix has the following format:
%   	steps(1,:) = step time.
%   	steps(2,:) = step size.
%   	steps(3,:) = baseline.
%
%   The input data contain both displacement and time points:
%    dat == displacement and time data.
%           dat(1,:) = displacement data
%           dat(2,:) = time information
%    fig_num == figure number; 0 for no display

function rv = stepdetect(dat,fig_num)

% Label trajectories that experience photobleaching using median filtering
% and convultion with the Heaviside function.
%
% INPUT
%	- TRACES: Output from individualTrace - intensity traces of connected
%	events.
% OUTPUT
%	- B: Vector of 1s and 0s indicating if photobleaching occured or not.
%	A 1 indicates photobleaching and a 0 indicates no photobleaching.

% Default scan window size
winsize = 9;
% Dwell rejection
f_reject = 2;
% Set bleaching step to be negative
forward = -1;
% Use 3 step size initial guesses to scan for different step sizes
gstep = forward*[0.25 0.5 0.75]*(max(dat(1,:)-min(dat(1,:))));
% Set initial window and values
inarg = [gstep mean(dat(1,1:5))]; % [step size guess, baseline estimate]
rawdat = dat(:,1:winsize); % data window for step detection
tindx = 1; % time point index for scanning data
winref = 1; % current start reference point of a window
cntindx = 1; % count number for scanned points in the current winsize window
tflag = 1; % flag for updating chi_sq_min
step_num = 0;
steps = [];
% Residual noise level for step detection
medfiltered = filter_median(dat(1,:),floor(winsize/2));
noise_level = std(dat(1,:) - medfiltered)^2;
% Coarse-grain filtering to locate ROIs for steps
th = median(movvar(dat(1,:),winsize,'Endpoints','shrink'));
padsize = ceil(winsize/2);
padd = padarray(dat(1,:),[0 padsize],'symmetric');
fpadd = nlfilter(padd,[1 winsize],@thcheck,th);
f_win = fpadd(padsize+1:end-padsize);
padsize = ceil(winsize/2)*10;
padd = padarray(dat(1,:),[0 padsize],'symmetric');
fpadd = nlfilter(padd,[1 10*winsize],@thcheck,th);
f_largewin = fpadd(padsize+1:end-padsize);
mask = max([f_win; f_largewin]);
timevec = dat(2,:);
tloc = timevec(logical(mask));
stepindices = 1:length(dat(1,:));
sloc = stepindices(logical(mask));

while tindx <= length(tloc)
	tstep = tloc(tindx);
    sindx = sloc(tindx); % frame number
	%--- define Heaviside fitting function ---
    intensities = rawdat(1,:); timevecwin = rawdat(2,:); 
    chi2_heaviside = @(x)sum((intensities - x(1)*(timevecwin > tstep) - x(2)).^2);
	%--- try different step guesses ---
	% gstep(1)
	inarg(1) = gstep(1);
	outarg = fminsearch(chi2_heaviside,inarg,optimset('MaxFunEvals',1e20,'MaxIter',1e10));
	chi_sq1 = sum((rawdat(1,:)-outarg(1)*(rawdat(2,:)>tstep)-outarg(2)).^2);
	% gstep(2)
	inarg(1) = gstep(2);
	outarg = fminsearch(chi2_heaviside,inarg,optimset('MaxFunEvals',1e20,'MaxIter',1e10));
	chi_sq2 = sum((rawdat(1,:)-outarg(1)*(rawdat(2,:)>tstep)-outarg(2)).^2);
	% gstep(3) 
	inarg(1) = gstep(3);
	outarg=fminsearch(chi2_heaviside,inarg, optimset('MaxFunEvals',1e20,'MaxIter',1e10));
    chi_sq3=sum((rawdat(1,:)-outarg(1)*(rawdat(2,:)>tstep)-outarg(2)).^2);
    chi_sq = min([chi_sq1 chi_sq2 chi_sq3]);
    dat_num = length(rawdat(1,:)); % number of points used in calculating the chi square.
    % This is used to normalize the global noise level.

    if tflag == 1
    	min_chi_sq = chi_sq; 
    end

    if chi_sq < min_chi_sq || tflag == 1 % no step detected
        inarg = outarg; % update min_chi_sq
        min_chi_sq = chi_sq;
        cntindx = cntindx + 1;
        tflag = 2; % update scan index
        if cntindx == winsize % update window flag for step scanning
            tflag = 1; % resize window when scan index reaches the end of scan window
            cntindx = 1;
            if sindx + winsize - 1 >= length(dat(2,:))
                rawdat = dat(:,winref:end);
            else
                rawdat = dat(:,winref:sindx + winsize - 1);
            end
        end
        % Step detection criteria are set here:
        % (1) step has to last more than f_reject frames
        % (2) chi_sq has to increase by more than noise_level above the min_chi_sq
        % Step detected
    elseif chi_sq > (min_chi_sq + noise_level * dat_num) && cntindx > f_reject
        step_num = step_num + 1;
        if step_num == 1 % record first stem
            dstep = outarg(1,1);
            baseline = mean(dat(1,1:sindx - 2)); % refine the first baseline
        else
            baseline = mean(dat(1,(sum((dat(2,:) < steps(1,step_num - 1))) + 2):sindx - 2));
            steps(2,step_num - 1) = baseline - steps(3,step_num - 1); % refine previous step size
        end
        inarg = [dstep, baseline + dstep]; % set initial guesses for the next round
        steps = [steps [mean(dat(2,sindx + [-1, 0])); dstep; baseline]]; % record a step detected
        
        tflag = 1; % reset windowflag for a new window
        winref = sindx;
        cntindx = 1; % reset scan index
        if sindx + winsize - 1 >= length(dat(2,:))
            rawdat = dat(:,sindx:length(dat(2,:)));
        else
            rawdat = dat(:,sindx:sindx+winsize-1);
        end
    else % no step is detected
        inarg = outarg; % no change in min_chi_sq
        cntindx = cntindx + 1;
        tflag = 2;
        if cntindx == winsize
            tflag = 1;
            cntindx = 1;
            if sindx + winsize - 1 >= length(dat(2,:))
                rawdat = dat(:,winref:length(dat(2,:)));
            else
                rawdat = dat(:,winref:sindx + winsize - 1);
            end
        end
    end
    tindx = tindx + 1; % advance to the next frame within the while loop
end

if step_num > 0
    % refine the last step size
    steps(2,step_num) = mean(dat(1,dat(2,:) > steps(1,step_num))) - steps(3,step_num);
    % filter out the steps that are increasing in intensity (photobleaching is strictly decreasing in intensity)
    filter_matrix = steps(2,:) < 0;
    filter_matrix = repmat(filter_matrix, [size(steps(:,1)),1]);
    steps = reshape(steps(logical(filter_matrix)),[size(steps(:,1)), sum(filter_matrix(1,:))]);
    step_num = length(steps(1,:));
end

if step_num > 0
    % remove extra deep steps
    for k = 1:length(steps(3,:)) - 1
        steps(2,k) = steps(3,k + 1) - steps(3,k);
    end
    
    % remove small steps that are smaller than 1x standard deviation of local noise level
    % get local noise level backing up noise_window steps
    x_range = steps(1,:);
    noise_window = winsize;
    filter_matrix = ones(size(steps(2,:)));
    for l = 1:length(x_range)
        if (x_range(l) - noise_window) <= 0
            start_point = ceil(x_range(l));
            end_point = floor(x_range(l) + noise_window);
        else
            start_point = ceil(x_range(l) - noise_window);
            end_point = floor(x_range(l));
        end
        local_noise = std(dat(1,start_point:end_point));
        if abs(steps(2,l)) <= abs(local_noise)
            filter_matrix(l) = 0;
        end
    end
    
    filter_matrix = repmat(filter_matrix, [size(steps(:,1)), 1]);
    steps = reshape(steps(logical(filter_matrix)),[size(steps(:,1)), sum(filter_matrix(1,:))]);
    
    % Update step size to account for deletion of small stpes
    % Remove extra deep steps
    for p = 1:length(steps(3,:)) - 1
        steps(2,p) = steps(3,p+1) - steps(3,p);
    end
    step_num = length(steps(1,:));
end

% Plot the detected steps and display step data
if step_num > 0 && fig_num  > 0
    line([dat(2,1) steps(1,1)],steps(3,1)*ones(1,2),'linewidth',2);
    for indx=1:step_num 
        line(steps(1,indx)*ones(1,2),[steps(3,indx) sum(steps(2:3,indx))],'linewidth',2);
        if indx==step_num 
            line([steps(1,indx) dat(2,end)],sum(steps(2:3,indx))*ones(1,2),'linewidth',2);
        else
            line([steps(1,indx) steps(1,indx+1)],steps(3,indx+1)*ones(1,2),'linewidth',2);
        end
    end
    % disp(steps);
end

%==== return detected steps ====
step_details = {};
if step_num > 0
    step_details = {steps(1,:),steps(2,:),steps(3,:)};
end
rv = {step_num,step_details};
% steps(1,:) = step time.
% steps(2,:) = step size.
% steps(3,:) = baseline.
%===============================


