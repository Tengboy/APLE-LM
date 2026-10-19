function [pos_coarse,pos_para_coarse,H_o] = pos_coarse_appro_est_joint(iterPara,montPara,sysPara)
% ACM counterpart: phase-only channel; same verified update schedule and state fixes.
% Joint coarse update used once by main_240803_old.
% Keep the original sign-gradient schedule: angle, range, range.
pos_coarse = iterPara.pos_coarse;
pos_para_coarse = iterPara.pos_para_coarse;
q = pos_para_coarse(:);
gradient = @(x)grad_calculate_joint1(x,montPara.y_compress_noise, ...
    sysPara.W,montPara.inv_C_n,sysPara.pos_element,sysPara.lambda);
curvature = @(x)hess_calculate_joint1(x,montPara.y_compress_noise, ...
    sysPara.W,montPara.inv_C_n,sysPara.pos_element,sysPara.lambda);
if isfield(iterPara,'cal_Hess') && iterPara.cal_Hess == 1
    H_o = curvature(q);
    return;
end
objective = @(x)MCRB_polor_obj_joint_approx(x,sysPara,montPara);
alpha_angle = 0.01;
alpha_range = 0.1;
range_step = 1;
maximum = sysPara.counter_coarse_max2;
% As in the original, stopping is checked on range updates. The counter
% starts at 1 here (the preceding per-user warm start starts at 2).
for counter = 1:maximum+2
    direction = reshape(sign(gradient(q)),3,[]);
    parameters = reshape(q,3,[]);
    angle_update = mod(counter,3)==1;
    if angle_update
        alpha_angle = armijo_part1(q,-direction(:),objective,1e-6,2*alpha_angle);
        parameters(1:2,:) = parameters(1:2,:)-alpha_angle*direction(1:2,:);
    else
        alpha_range = armijo_part2(q,-direction(:),objective,1e-6,2*alpha_range);
        range_step = alpha_range*direction(3,:);
        parameters(3,:) = parameters(3,:)-range_step;
    end
    q = parameters(:);
    position = [parameters(3,:).*sqrt(1-sum(parameters(1:2,:).^2,1)); ...
        parameters(3,:).*parameters(1:2,:)];
    stopping = ~angle_update && (norm(range_step)/norm(position,'fro') ...
        <sysPara.threshold_coarse || counter>maximum);

    % Preserve the two original Newton opportunities and their step limits.
    % A Newton direction need not descend when the Hessian is indefinite.
    if stopping || (max(abs(range_step))<1e-2 && counter>15)
        H_step = curvature(q);
        delta = H_step\gradient(q);
        if stopping
            limit = (sysPara.rmax-sysPara.rmin)/2;
        else
            limit = 0.5;
        end
        delta(abs(delta)>limit) = 0;
        candidate = reshape(q-delta,3,[]);
        if isreal(candidate) && all(isfinite(candidate),'all') && ...
                all(sum(candidate(1:2,:).^2,1)<1) && all(candidate(3,:)>0) && ...
                objective(candidate(:))<objective(q)
            q = candidate(:);
        end
    end
    if stopping, break; end
end
pos_para_coarse = reshape(q,3,[]);
pos_coarse = [pos_para_coarse(3,:).*sqrt(1-sum(pos_para_coarse(1:2,:).^2,1)); ...
    pos_para_coarse(3,:).*pos_para_coarse(1:2,:)];
% The original returned the Hessian from BEFORE its last Newton step.
% Fusion requires curvature evaluated at the exact returned coarse mean.
H_o = curvature(q);
end
