1/29/18: Edit to ft_checkdata line 874 to increase speed.

2/14/18 Deleted Matlab scripts folders in folder compat for incompatibility with PCT

2/16/18 
	commented out line 62-64 in univariate2bivariate for incompatibility with fourierspctrm for connectivity.

2/21/18 
	added " ...,'allowoverlap',true " for ft_definetrial line 248 for compatibilities with artifact rejection 

3/15/18 changed topoplot_common.m to " msk = mean(msk(sellab, xmin:xmax),2);" at line 604.
	   ft_plot_topo change 199 to become "maskimage = (maskimage .* maskimagetmp)
"

commented out the old 199 and line 241

added lines 253-255 "  h.AlphaData = maskimage;
 h.FaceColor = 'texturemap';
 h.FaceAlpha = 'texturemap';
"

added 253 "  if ~all(diff(maskimage(maskimage ~= 0)))
      maskimagetmp = find(maskimage(maskimage ~= 0));
      maskimage(maskimagetmp(1)) = 1;
  end
"
for issue when all nonzero alpha values were equal

saving of persistent maskimage was being incompatiable with changing alphadata;
changed line 159 to "if isequal(current_argin, previous_argin) || ~isempty(mask)
  % convert the mask into a binary image
  maskimage = zeros(gridscale);%false(gridscale);
"

3/30/18 

Added SCC capabilities to ft_freqanalysis
line 230 - 238 
"if isfield(cfg,'scc')
    if cfg.scc == 1
       scc = 1 ;
    else
        scc = 0;
    end
else
    scc = 0;
end
"

all ft_specest functions also edited to have gpuarray(fft)

4/5/18

Added read_brainvision_eeg.mat to have function to preprocess faster

4/6/18
add more gpuarrays in ft_specest wavelet for better GPU compatibilities.

6/6
fast-forwarded a few scripts from newer fieldtrip
topoplot_common with our previous change on 3/15/18
updated ft_warning
added new ft_notification
added defaultId
updated ft_defaults
updated ft_platform_supports

6/12 
updated ft_singleplotER/TFR

7/23/18
added .oldfoi in ft_freqanalysis to change pad-modified .freq back to regular range.
added ft_error 

8/3/2018
Modified SCD to be in uV/m^2 units instead of the V/m^2 units
