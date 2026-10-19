function [positions,q,info]=fig5_psigw(q,mont,sys)
% Historical repository P-SIGW refinement, with explicit diagnostics.
% Retain sign directions, one angle/two range updates, range-only stopping,
% and the historical counter > limit check. This is not a new SIGW algorithm.
objective=@(x)MCRB_polor_obj_joint_bias_new(x,sys,mont);
alpha_angle=.01; alpha_range=.1;
trace=nan(sys.counter_coarse_max2+3,1); trace(1)=objective(q(:));
angle_step=NaN; relative=Inf;
for iteration=1:sys.counter_coarse_max2+2
    gradient=grad_calculate_joint1_bias(q(:),mont.y_compress_noise,sys.W, ...
        mont.inv_C_n,sys.pos_element,sys.lambda,sys.array_n);
    direction=sign(gradient);
    delta=reshape(direction,3,[]);
    if mod(iteration,3)==1
        alpha_angle=armijo_part1(q(:),-direction,objective,1e-6,2*alpha_angle);
        delta(1:2,:)=alpha_angle.*delta(1:2,:); delta(3,:)=0;
        angle_step=norm(delta,'fro');
    else
        alpha_range=armijo_part2(q(:),-direction,objective,1e-6,2*alpha_range);
        delta(1:2,:)=0; delta(3,:)=alpha_range.*delta(3,:);
    end
    q=reshape(q(:)-delta(:),3,[]);
    trace(iteration+1)=objective(q(:));
    normal=sqrt(1-q(1,:).^2-q(2,:).^2);
    positions=[q(3,:).*normal;q(3,:).*q(1,:);q(3,:).*q(2,:)];
    if mod(iteration,3)~=1
        relative=norm(delta,'fro')/norm(positions,'fro');
        if relative<sys.threshold_coarse || iteration>sys.counter_coarse_max2
            break;
        end
    end
end
assert(isreal(positions) && all(isfinite(positions),'all'),'Invalid P-SIGW position.');
info=struct('iterations',iteration,'stopped_by_tolerance',relative<sys.threshold_coarse, ...
    'iteration_cap',iteration>sys.counter_coarse_max2,'relative_range_step',relative, ...
    'last_angle_step',angle_step,'objective_trace',trace(1:iteration+1), ...
    'stopping_rule','historical range-only relative step; not full gradient stationarity');
end
