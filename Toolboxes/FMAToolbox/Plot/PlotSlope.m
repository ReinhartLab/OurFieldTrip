function p = PlotSlope(x,y,slope,len,varargin)

%PlotSlope - Plot line with given slope.
%
%  USAGE
%
%    p = PlotSlope(x,y,slope,length,<options>)
%
%    x,y            coordinates of central point
%    slope          line slope
%    length         segment length along the x axis
%    <options>      line and color specification (see <a href="matlab:help plot">plot</a>)
%

% Copyright (C) 2010-2016 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 4,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotSlope">PlotSlope</a>'' for details).');
end

X = x+[-len len]/2;
Y = y+[-len len]/2*slope;
plot(X,Y,varargin{:});
