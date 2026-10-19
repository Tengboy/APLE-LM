function items=aple_figure_catalog
% Version 0.2 figure IDs: journal order, with CDF directly after the SNR plot.
keys={'grid','snr','cdf','aperture','rf','range','coverage','iterations'};
suffix={'grid_precision_reliability','snr_performance','error_distribution', ...
    'array_grouping','rf_precision_reliability','distance_performance','spatial_coverage','iteration_behavior'};
families={'grid','snr','snr','aperture','rf','range','coverage','iterations'};
old={'M3','M1','M2','M4(a)','M5','M4(b)','M7','M6'};
journal={'Fig. 4','Figs. 5-6','CDF supplement','Fig. 7','Fig. 8','Fig. 9','Fig. 10','Fig. 11'};
items=struct('code',{},'key',{},'filename',{},'family',{},'previous_code',{},'journal',{});
for k=1:8
    items(k)=struct('code',sprintf('M%d',k),'key',keys{k},'filename',sprintf('M%d_%s',k,suffix{k}), ...
        'family',families{k},'previous_code',old{k},'journal',journal{k});
end
end
