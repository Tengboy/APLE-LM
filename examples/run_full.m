% This is a long Monte Carlo computation. Repeat to resume after interruption.
setup_aple;
out=run_experiments('all','full','BatchSize',25);
plot_results(out);
% Optional: run_experiments('all','full','UseParallel',true,'Workers',2);
