function spikes = GetSpikeTimes(units,varargin)

%GetSpikeTimes - Get spike timestamps.
%
%  USAGE
%
%    spikes = GetSpikeTimes(units,<options>)
%
%    units          optional list of units, i.e. [electrode group, cluster] pairs;
%                   special conventions:
%                     'all'          all units
%                     []             no units
%                     cluster = -1   all clusters except artefacts and MUA
%                     cluster = -2   all clusters except artefacts
%                     cluster = -3   all clusters
%                   (artefacts are assumed to be in cluster 0, and MUA in 1)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'output'      'time' returns only timestamps, 'full' returns timestamps,
%                   electrode groups and clusters, and 'numbered' returns
%                   timestamps and arbitrary unique identifiers corresponding
%                   to single units (see EXAMPLE below) (default = 'time')
%    =========================================================================
%
%  EXAMPLES
%
%    % timestamps for all spikes
%    s = GetSpikeTimes;
%    % or
%    s = GetSpikeTimes('all');
%
%    % timestamps for units [1 7] and [4 3]
%    s = GetSpikeTimes([1 7;4 3]);
%
%    % timestamps for all single units on electrode group 5, and unit [6 3]
%    s = GetSpikeTimes([5 -1;6 3]);
%
%    % timestamps for all units on electrode group 5, except artefacts
%    s = GetSpikeTimes([5 -2]);
%
%    % timestamps, electrode groups and clusters, for all spikes
%    s = GetSpikeTimes('output','full');
%    % or
%    s = GetSpikeTimes('all','output','full');
%
%    % timestamps and identifiers, for units [1 7], [4 3] and [2 5]
%    % unit [1 7] will be assigned number 1, unit [2 5] number 2, and
%    % unit [4 3] number 3, independent from the order in which they are listed
%    s = GetSpikeTimes([1 7;4 3;2 5],'output','numbered');
%
%  NOTE
%
%    An electrode group is an ensemble of closely spaced electrodes that record from
%    the same neurons, e.g. a single wire electrode, or a wire tetrode, or a multisite
%    silicon probe in octrode configuration, etc.
%
%  SEE
%
%    See also LoadSpikeTimes, PlotTicks.


% Copyright (C) 2004-2017 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Default values
output = 'time';
if nargin == 0, units = 'all'; end

% Optional parameter
if ischar(units) && ~strcmp(units,'all'),
	varargin = {units,varargin{:}};
	units = 'all';
else
	if ~strcmp(units,'all') && ~isempty(units) && (~isimatrix(units) || size(units,2) ~= 2),
		error('Incorrect list of units (type ''help <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a>'' for details).');
	end
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'output',
			output = lower(varargin{i+1});
			if ~isastring(output,'time','full','numbered'),
				error('Incorrect value for property ''output'' (type ''help <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a>'' for details).']);

  end
end

if isempty(units),
	spikes = [];
else
	spikes = DATA.spikes;
end
% Adjust output matrix size
if isempty(spikes),
	switch(output),
		case 'time',
			spikes = nan(0,1);
		case 'numbered',
			spikes = nan(0,2);
		case 'full',
			spikes = nan(0,3);
	end
	return
end

% Selected units only
if ~isastring(units,'all'),
	nUnits = size(units,1);
	selected = zeros(size(spikes(:,1)));
	for i = 1:nUnits,
		channel = units(i,1);
		cluster = units(i,2);
		switch cluster,
			case -1,
				selected = selected | (spikes(:,2) == channel & spikes(:,3) ~= 0 & spikes(:,3) ~= 1);
			case -2,
				selected = selected | (spikes(:,2) == channel & spikes(:,3) ~= 0);
			case -3,
				selected = selected | spikes(:,2) == channel;
			otherwise,
				selected = selected | (spikes(:,2) == channel & spikes(:,3) == cluster);
		end
	end
	spikes = spikes(selected,:);
end

if strcmp(output,'time'),
	spikes = spikes(:,1);
elseif strcmp(output,'numbered'),
	[units,~,i] = unique(spikes(:,2:end),'rows');
	nUnits = length(units);
	index = 1:nUnits;
	id = index(i)';
	spikes = [spikes(:,1) id];
end

spikes = sortrows(spikes,1);
