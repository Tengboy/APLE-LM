function bundle=aple_system(d,case_id)
% One active dictionary in RAM; no system/cache MAT files are written.
thread_guard=aple_threads; %#ok<NASGU>
c=d.cases(case_id);
b=struct('counter_ini_max',20,'counter_coarse_max1',20,'counter_coarse_max2',200, ...
    'counter_fine_max',20,'threshold_ini',1e-3,'threshold_coarse',1e-4,'threshold_fine',1e-5, ...
    'iter_max',200,'freq',6e9,'damping_fact',.7,'lambda',.05,'array_dir_x',[0;1;0], ...
    'array_dir_y',[0;0;1],'array_n',[1;0;0],'element_gap',.025, ...
    'source_num',c.users,'rmin',c.range(1),'rmax',c.range(2),'theta_range',120);
sys=build_array_system(b,c.M,c.Ns,c.RF,c.W_seed,c.grid);
bundle=struct('accurate',sys);
if any(strcmp(d.methods,'APLE-LM-ACM'))
    acm=sys;
    [~,acm.steer_vec_noise_inv,acm.steer_vec_noise_inv2]= ...
        prepare_grid_dictionary(sys.pos_grid,sys.W,sys.pos_element,sys.lambda,sys.array_n,true);
    bundle.acm=acm;
end
if any(ismember(d.methods,{'ES-GA','P-SIGW'}))
    baseline=fig5_full_baseline(sys); bundle.esga=baseline; bundle.psigw=baseline;
    bundle.esga.counter_ini_max=200; bundle.psigw.counter_ini_max=0;
end
end
