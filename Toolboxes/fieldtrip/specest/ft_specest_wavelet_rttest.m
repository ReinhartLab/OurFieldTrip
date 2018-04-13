function [spectrum] = ft_specest_wavelet_rttest(cfg,data, varargin)

% FT_SPECEST_WAVELET performs time-frequency analysis on any time series trial data
% using the 'wavelet method' based on Morlet wavelets, doing convolution in the time
% domain by multiplication in the frequency domain.
%
% Use as
%   [spectrum,freqoi,timeoi] = ft_specest_wavelet(dat,time...)
% where
%   dat       = matrix of chan*sample
%   time      = vector, containing time in seconds for each sample
%   spectrum  = array of chan*freqoi*timeoi of fourier coefficients
%   freqoi    = vector of frequencies in spectrum
%   timeoi    = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include
%   pad       = number, total length of data after zero padding (in seconds)
%   padtype   = string, indicating type of padding to be used (see ft_preproc_padding, default = 'zero')
%   freqoi    = vector, containing frequencies of interest
%   timeoi    = vector, containing time points of interest (in seconds)
%   width     = number or vector, width of the wavelet, determines the temporal and spectral resolution
%   gwidth    = number, determines the length of the used wavelets in standard deviations of the implicit Gaussian kernel
%   verbose   = output progress to console (0 or 1, default 1)
%   polyorder = number, the order of the polynomial to fitted to and removed from the data prior to the fourier transform (default = 0 -> remove DC-component)
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMCONVOL, FT_SPECEST_TFR, FT_SPECEST_HILBERT, FT_SPECEST_MTMFFT

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

dat = cell2mat(permute(data.trial,[1 3 2]));
time = data.time{1};

ntrials = numel(data.trial);
trllength = zeros(1, ntrials);
for itrial = 1:ntrials
  trllength(itrial) = size(data.trial{itrial}, 2);
end
if strcmp(cfg.pad, 'maxperlen')
  padding = max(trllength);
  pad = padding/data.fsample;
elseif strcmp(cfg.pad, 'nextpow2')
  padding = 2^nextpow2(max(trllength));
  pad = padding/data.fsample;
else
  padding = cfg.pad*data.fsample;
  if padding<max(trllength)
    error('the specified padding is too short');
  end
end

freqoi       = ft_getopt(cfg, 'foi',       []);
fboi   = round(freqoi .* pad) + 1;
freqoi   = (fboi-1) ./ pad; % boi - 1 because 0 Hz is included in fourier output

width  = ft_getopt(cfg, 'width',  7);
gwidth = ft_getopt(cfg, 'gwidth', 3);
scc   = ft_getopt(cfg, 'scc', false);


% get the optional input arguments

timeoi    = ft_getopt(cfg, 'toi', 'all');

padtype   = ft_getopt(cfg, 'padtype',   'zero');
polyorder     = ft_getopt(cfg, 'polyremoval', 0);




% Set n's
[nchan,ndatsample,ntrials] = size(dat);

% This does not work on integer data
typ = class(dat);
if ~strcmp(typ, 'double') && ~strcmp(typ, 'single')
  dat = cast(dat, 'double');
end

% Remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% Determine fsample and set total time-length of data
fsample = 1./mean(diff(time));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
postpad    = round((pad - dattime) * fsample);
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data

% Set freqboi and freqoi
freqoiinput = freqoi;
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all') % if input was 'all'
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
% check for freqoi = 0 and remove it, there is no wavelet for freqoi = 0
if freqoi(1)==0
  freqoi(1)  = [];
  freqboi(1) = [];
end
nfreqboi = length(freqboi);
nfreqoi  = length(freqoi);

% throw a warning if input freqoi is different from output freqoi
if isnumeric(freqoiinput)
  % check whether padding is appropriate for the requested frequency resolution
  rayl = 1/endtime;
  if any(rem(freqoiinput,rayl)) % not always the case when they mismatch
    ft_warning('padding not sufficient for requested frequency resolution, for more information please see the FAQs on www.ru.nl/neuroimaging/fieldtrip');
  end
  if numel(freqoiinput) ~= numel(freqoi) % freqoi will not contain double frequency bins when requested
    ft_warning('output frequencies are different from input frequencies, multiples of the same bin were requested but not given');
  else
    if any(abs(freqoiinput-freqoi) >= eps*1e6)
      ft_warning('output frequencies are different from input frequencies');
    end
  end
end

% Set timeboi and timeoi
timeoiinput = timeoi;
offset = round(time(1)*fsample);
if isnumeric(timeoi) % if input is a vector
  timeoi   = unique(round(timeoi .* fsample) ./ fsample);
  timeboi  = round(timeoi .* fsample - offset) + 1;
  ntimeboi = length(timeboi);
elseif strcmp(timeoi,'all') % if input was 'all'
  timeboi  = 1:length(time);
  ntimeboi = length(timeboi);
  timeoi   = time;
end

% throw a warning if input timeoi is different from output timeoi
% if isnumeric(timeoiinput)
%   if numel(timeoiinput) ~= numel(timeoi) % timeoi will not contain double time-bins when requested
%     ft_warning('output time-bins are different from input time-bins, multiples of the same bin were requested but not given');
%   else
%     if any(abs(timeoiinput-timeoi) >= eps*1e6)
%       ft_warning('output time-bins are different from input time-bins');
%     end
%   end
% end


% Compute fft
if ~scc
    spectrum = complex(nan(nchan,nfreqoi,ntimeboi,ntrials),nan(nchan,nfreqoi,ntimeboi,ntrials));
    dat(:,end+1:end+postpad,:) = zeros([nchan postpad ntrials]);
    datspectrum = fft(dat, [], 2);
else
    spectrum = complex(nan(nchan,nfreqoi,ntimeboi),nan(nchan,nfreqoi,ntimeboi));
    datspectrum = fft(gpuArray(ft_preproc_padding(dat, padtype, 0, postpad)), [], 2);
end


% Creating wavelets
% expand width to array if constant width
if numel(width) == 1
  width = ones(1,nfreqoi) * width;
end

  dt = 1/fsample;
  sf = freqoi ./ width;
  st = 1 ./(2*pi*sf);
  A = 1 ./sqrt(st*sqrt(pi));
for ifreqoi = 1:nfreqoi
    
%   str = sprintf('frequency %d (%.2f Hz)', ifreqoi,freqoi(ifreqoi));
%     

  toi2 = -gwidth*st(ifreqoi):dt:gwidth*st(ifreqoi);
  tap = (A(ifreqoi) *exp(-toi2.^2/(2*st(ifreqoi)^2)))';
  acttapnumsmp = size(tap,1);
  ins = ceil(endnsample./2) - floor(acttapnumsmp./2);
  prezer = zeros(ins,1);
  pstzer = zeros(endnsample - ((ins-1) + acttapnumsmp)-1,1);
  
  % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
  ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./fsample) .* freqoi(ifreqoi));
  
  % create wavelet and fft it
 
%   [st, cws] = dbstack;
%   if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis') && verbose
%     % specest_convol has been called by ft_freqanalysis, meaning that ft_progress has been initialised
%     ft_progress(fbopt.i./fbopt.n, ['trial %d, ',str,'\n'], fbopt.i);
%   elseif verbose
%     fprintf([str, '\n']);
%   end
  
  % compute indices that will be used to extracted the requested fft output
  reqtimeboiind    = find((timeboi >=  (acttapnumsmp ./ 2)) & (timeboi < (ndatsample - (acttapnumsmp ./2))));
  reqtimeboi       = timeboi(reqtimeboiind);
  
  % compute datspectrum*wavelet, if there are reqtimeboi's that have data
  if ~isempty(reqtimeboi)
      if ~scc
          dum = [fftshift(ifft(datspectrum .* repmat(fft(complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer)),[],1)',[nchan 1]), [], 2),2)] .* sqrt(2 ./ fsample);
          spectrum(:,ifreqoi,reqtimeboiind,:) = dum(:,reqtimeboi,:);
      else
          dum = gather([fftshift(ifft(datspectrum .* repmat(fft(gpuArray(complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer))),[],1)',[nchan 1]), [], 2),2)] .* sqrt(2 ./ fsample));
          spectrum(:,ifreqoi,reqtimeboiind,:) = dum(:,reqtimeboim,:);
      end
  end
end

spectrum = permute(spectrum,[4 1 2 3]);
cumtapcnt = ones([1 nfreqoi]);

freq        = [];
freq.label  = data.label;
freq.dimord = dimord;
freq.freq   = foi;
