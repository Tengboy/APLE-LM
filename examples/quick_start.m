% Run from the package root after extracting the complete folder.
setup_aple;
plot_results('reference','Figures',{'M5','M8'});
out=run_experiments('M2','smoke');
plot_results(out,'Figures',{'M2','M3'});
