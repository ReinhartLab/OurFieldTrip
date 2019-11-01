% General-purpose data processing and statistical tests for FMAToolbox.
%
% Bins and histograms.
%
%   Accumulate                     - Accumulate repeated observations.
%   Bin                            - Assign each input value a bin number between 1 and N.
%
% General Statistics.
%
%   sem                            - Compute standard error of the mean (SEM).
%   semedian                       - Compute standard error of the median.
%   nansem                         - Compute standard error of the mean (SEM), ignoring NaNs.
%   mean2str                       - Convert mean and SEM to string.
%   npcdf                          - Non-parametric cumulative distribution function.
%   BartlettTest                   - Test if k groups of samples have equal variances (homogeneity of variances).
%   FisherTest                     - Test if two groups of samples have equal variances.
%   MultinomialConfidenceIntervals - Simultaneous multinomial confidence intervals.
%   CompareSlopes                  - Perform linear regression on two groups and compare slopes.
%
% Circular Statistics.
%
%   CircularANOVA                  - One or two-way ANOVA on circular data.
%   CircularConfidenceIntervals    - Compute circular mean and confidence intervals of circular data.
%   CircularDistribution           - Compute circular distribution and statistics.
%   CircularMean                   - Estimate the circular mean.
%   CircularRegression             - Linear-circular regression.
%   CircularVariance               - Estimate circular variance and standard deviation.
%   Concentration                  - Estimate the concentration parameter for circular data.
%   ConcentrationTest              - Test homogeneity of concentration parameters.
%   WatsonU2Test                   - Test if two samples (circular data) have different means / variances.
%
% Samples: selection.
%
%   Restrict                       - Keep only samples that fall in a given list of time intervals.
%   InIntervals                    - Test which values fall in a list of intervals.
%   ConsolidateIntervals           - Consolidate intervals.
%   SubtractIntervals              - Subtract intervals.
%   ExcludeIntervals               - Exclude intersecting intervals.
%   FindInInterval                 - Find samples that fall in a given interval.
%   CountInIntervals               - Count samples that fall in each of a list of intervals.
%   ToIntervals                    - Convert logical vector to a list of intervals.
%   IsExtremum                     - Identify local maxima or minima.
%   ZeroCrossings                  - Test zero crossings in a given time series.
%   IsFirstAfter                   - Identify first item after each of a list of timestamps.
%   IsLastBefore                   - Identify last item before each of a list of timestamps.
%   Threshold                      - Find periods above/below threshold.
%
% Samples: simple processing.
%
%   Detrend                        - Detrend signal.
%   Diff                           - Differentiate samples.
%   Filter                         - Filter samples.
%   Interpolate                    - Interpolate samples (positions, spikes, LFP, etc.) at given timestamps.
%
% General-purpose.
%
%   Array2Matrix                   - Transform N-dimensional array into a matrix.
%   Array2PagedMatrix              - Transform an N-dimensional array into a paged matrix.
%   CircularShift                  - Shift matrix rows or columns circularly.
%   Clip                           - Clip values.
%   CumSum                         - Cumulative sum of elements. Partial sums can also be computed.
%   DistanceTransform              - Compute distance transform of a binary matrix (distance to nearest 1 value).
%   ExtendArray                    - Extend existing array to a given (larger) size.
%   Insert                         - Insert lines in a matrix.
%   Match                          - Replace values in one list with closest values in a second list.
%   MatchPairs                     - Pair nearest values in two lists.
%   RunningAverage                 - Compute running linear or angular average.
%   SineWavePeaks                  - Find peaks (or troughs) in a sine wave.
%   Smooth                         - Smooth using a Gaussian kernel.
%   XCorr1                         - Compute cross-correlograms and their modes for pairs of curves.
%   ZeroToOne                      - Normalize values in [0,1].
%
