function [positions,info]=fig5_estimate(mont,sys,use_acm)
% The reviewed main_240803_old loop, with the same ACM counterpart schedule.
if use_acm
    initialize=@pos_ini_appro_est; coarse=@pos_coarse_appro_est;
    joint=@pos_coarse_appro_est_joint; gain=@rho_appro_est; fine=@pos_fine_appro_est;
else
    initialize=@pos_ini_est; coarse=@pos_coarse_est;
    joint=@pos_coarse_est_joint; gain=@rho_est; fine=@pos_fine_est;
end
total=tic; timer=tic;
[p,q]=initialize(mont,sys); initialization_seconds=toc(timer);
it=struct('pos_coarse',p,'pos_para_coarse',q,'pos_fine',p,'pos_para_fine',q, ...
    'pos_damping',p,'cal_Hess',0,'turbo_iter_counter',0);
previous=p; coarse_seconds=0; fine_seconds=0;
for iteration=1:sys.iter_max
    it.turbo_iter_counter=iteration;
    if iteration==1
        timer=tic;
        [it.pos_coarse,it.pos_para_coarse]=coarse(it,mont,sys);
        [it.pos_coarse,it.pos_para_coarse,it.Hess_para_coarse]=joint(it,mont,sys);
        it.pos_fine=it.pos_coarse; it.pos_para_fine=it.pos_para_coarse;
        it.pos_damping=it.pos_coarse; coarse_seconds=toc(timer);
    end
    timer=tic;
    [it.rho_est,it.rho_invcov]=gain(it,mont,sys);
    [it.pos_fine,it.pos_para_fine]=fine(it,mont,sys);
    damping=max(sys.damping_fact^iteration,.01);
    it.pos_damping=(1-damping)*it.pos_damping+damping*it.pos_fine;
    change=norm(it.pos_fine-previous,'fro')/norm(previous,'fro');
    fine_seconds=fine_seconds+toc(timer);
    if change<1e-4 || iteration>=sys.iter_max, break; end
    previous=it.pos_fine;
end
it.cal_Hess=1;
[~,~,it.Hess_para_fine]=fine(it,mont,sys);
[~,coarse_status]=chol(it.Hess_para_coarse); [~,fine_status]=chol(it.Hess_para_fine);
positions=pos_comb_est(it);
assert(isreal(positions) && all(isfinite(positions),'all'),'APLE:InvalidFusion','Non-finite or complex fusion position.');
info=struct('seconds',toc(total),'initialization_seconds',initialization_seconds, ...
    'coarse_seconds',coarse_seconds,'fine_seconds',fine_seconds,'iterations',iteration, ...
    'stopped_by_tolerance',change<1e-4,'relative_change',change, ...
    'coarse_positive_definite',coarse_status==0,'fine_positive_definite',fine_status==0);
end

