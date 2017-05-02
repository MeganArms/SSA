function [lifetimes, censoring] = datacleanup(d0,dstart,obj,th,ff,et,exclude)

% DATACLEANUP creates the right data to input to the EDCF function.
% INPUT
%	D0 - experiment start time in the format 'yyyyMMdd HHmmss'
%	DSTART - observation start time in the format 'yyyyMMdd HHmmss'
%	OBJ - SSA object
%	TH - distance threshold (in pixels) for connecting two surface events
%	ET - time between the start of exposures, in seconds
%	FF - the minimum number of frames that an object must appear on to be
%	included in the lifetime data.
% 	EXCLUDE - exclude events of unknown lifetime ('yes') or include ('no')
% OUTPUT
% 	LIFETIMES - lifetime data to input into ecdf
% 	CENSORING - censoring vector to input into ecdf

t0 = datetime(d0,'InputFormat','yyyyMMdd HHmmss');
tstart = datetime(dstart,'InputFormat','yyyyMMdd HHmmss');
L = datenum(tstart - t0)*24*60*60;

spots = connectEvents(obj,th);
if nargout > 1
	[C,censoring] = getCoordinates(obj,spots,exclude,ff);
else
	C = getCoordinates(obj,spots,exclude,ff);
end
lifetimes = cellfun(@length,C)*et + L;
end