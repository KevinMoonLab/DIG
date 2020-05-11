function Y = multi_scale(PDX, ndim, series_distance)
% CMDS
Y = randmds(PDX, ndim);
% MMDS
% if ~(strcmp(series_distance,"IEG") || strcmp(series_distance,"IEG2"))
% opt = statset('display','iter');
% Y = mdscale(PDX,ndim,'options',opt,'start',Y,'Criterion','metricstress');
% end 
end

