function [positions,trace]=fig11_trace(mont,sys,steps,take_snapshots)
% Figure 11 executes every requested update. Snapshot state never feeds back.
% take_snapshots=false is a regression control, used only by tests.
if nargin<3, steps=20; end
if nargin<4, take_snapshots=true; end
assert(steps>=1 && steps==floor(steps)); K=sys.source_num; P=3*K;
positions=nan(3,K,steps);
trace=struct('executed',false(1,steps),'executed_iterations',0,'first_tolerance_iteration',NaN, ...
    'fine_positions',nan(3,K,steps),'damped_positions',nan(3,K,steps), ...
    'fine_parameters',nan(3,K,steps),'coarse_parameters',nan(3,K), ...
    'coarse_hessian',nan(P,P),'fine_hessians',nan(P,P,steps), ...
    'relative_change',nan(1,steps),'coarse_positive_definite',false, ...
    'fine_positive_definite',false(1,steps),'error',{repmat({''},1,steps)}, ...
    'fatal_error','','initialization_seconds',NaN,'optimization_seconds',zeros(1,steps), ...
    'snapshot_seconds',zeros(1,steps));
timer=tic;
try
    [p,q]=pos_ini_est(mont,sys);
catch err
    trace.fatal_error=getReport(err,'extended','hyperlinks','off'); trace.error(:)={trace.fatal_error}; return;
end
trace.initialization_seconds=toc(timer);
it=struct('pos_coarse',p,'pos_para_coarse',q,'pos_fine',p,'pos_para_fine',q, ...
    'pos_damping',p,'cal_Hess',0,'turbo_iter_counter',0); previous=p;
for iteration=1:steps
    timer=tic; it.turbo_iter_counter=iteration;
    try
        if iteration==1
            [it.pos_coarse,it.pos_para_coarse]=pos_coarse_est(it,mont,sys);
            [it.pos_coarse,it.pos_para_coarse,it.Hess_para_coarse]=pos_coarse_est_joint(it,mont,sys);
            it.pos_fine=it.pos_coarse; it.pos_para_fine=it.pos_para_coarse; it.pos_damping=it.pos_coarse;
            trace.coarse_parameters=it.pos_para_coarse; trace.coarse_hessian=it.Hess_para_coarse;
            [~,status]=chol(it.Hess_para_coarse); trace.coarse_positive_definite=status==0;
        end
        [it.rho_est,it.rho_invcov]=rho_est(it,mont,sys);
        [it.pos_fine,it.pos_para_fine]=pos_fine_est(it,mont,sys);
        damping=max(sys.damping_fact^iteration,.01);
        it.pos_damping=(1-damping)*it.pos_damping+damping*it.pos_fine;
        change=norm(it.pos_fine-previous,'fro')/norm(previous,'fro');
        assert(isreal(it.pos_fine) && all(isfinite(it.pos_fine),'all') && isfinite(change), ...
            'F11:InvalidState','An optimizer update produced an invalid state.');
        trace.executed(iteration)=true; trace.fine_positions(:,:,iteration)=it.pos_fine;
        trace.damped_positions(:,:,iteration)=it.pos_damping; trace.fine_parameters(:,:,iteration)=it.pos_para_fine;
        trace.relative_change(iteration)=change;
    catch err
        trace.fatal_error=getReport(err,'extended','hyperlinks','off');
        trace.error(iteration:end)={trace.fatal_error}; trace.optimization_seconds(iteration)=toc(timer); break;
    end
    trace.optimization_seconds(iteration)=toc(timer);
    if take_snapshots || iteration==steps
        timer=tic;
        try
            snapshot=it; snapshot.cal_Hess=1;
            [~,~,snapshot.Hess_para_fine]=pos_fine_est(snapshot,mont,sys);
            trace.fine_hessians(:,:,iteration)=snapshot.Hess_para_fine;
            [~,status]=chol(snapshot.Hess_para_fine); trace.fine_positive_definite(iteration)=status==0;
            combined=pos_comb_est(snapshot);
            assert(isreal(combined) && all(isfinite(combined),'all'),'F11:InvalidFusion','Invalid fusion snapshot.');
            positions(:,:,iteration)=combined;
        catch err
            trace.error{iteration}=getReport(err,'extended','hyperlinks','off');
        end
        trace.snapshot_seconds(iteration)=toc(timer);
    end
    % Deliberately continue even after the normal tolerance has been reached.
    previous=it.pos_fine;
end
trace.executed_iterations=nnz(trace.executed);
first=find(trace.relative_change<1e-4,1); if ~isempty(first), trace.first_tolerance_iteration=first; end
end
