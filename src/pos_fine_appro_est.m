function [pos_fine,pos_para_fine,H_o2] = pos_fine_appro_est(iterPara,montPara,sysPara) %#ok<INUSD>
% ACM counterpart: phase-only channel; same verified update schedule and state fixes.
% Joint fine-stage update, with the original sign-gradient/AO schedule.
pos_fine = iterPara.pos_fine;
pos_para_fine = iterPara.pos_para_fine;
rho = iterPara.rho_est;
Q = iterPara.rho_invcov;
q = pos_para_fine(:);
gradient = @(x)grad_calculate_joint2_new(x,rho,Q, ...
    sysPara.pos_element,sysPara.lambda);
curvature = @(x)hess_calculate_joint2_new(x,rho,Q, ...
    sysPara.pos_element,sysPara.lambda);
if iterPara.cal_Hess == 1
    H_o2 = curvature(q);
    return; % Curvature and the fusion mean now refer to exactly the same point.
end
objective = @(x)MCRB_polor_obj_joint2_approx(x,iterPara,sysPara);
alpha_angle = 0.01;
alpha_range = 0.1;
range_step = 1;
for counter = 2:sysPara.counter_fine_max+2
    direction = reshape(sign(gradient(q)),3,[]);
    parameters = reshape(q,3,[]);
    if mod(counter,2)==0
        alpha_angle = armijo_part1(q,-direction(:),objective,1e-6,2*alpha_angle);
        parameters(1:2,:) = parameters(1:2,:)-alpha_angle*direction(1:2,:);
    else
        alpha_range = armijo_part2(q,-direction(:),objective,1e-6,2*alpha_range);
        range_step = alpha_range*direction(3,:);
        parameters(3,:) = parameters(3,:)-range_step;
    end
    q = parameters(:);
    pos_fine = [parameters(3,:).*sqrt(1-sum(parameters(1:2,:).^2,1)); ...
        parameters(3,:).*parameters(1:2,:)];
    if mod(counter,2)~=0 && (norm(range_step)/norm(pos_fine,'fro')<sysPara.threshold_fine ...
            || counter>sysPara.counter_fine_max)
        break;
    end

    % Optional original Newton refinement. Accept only physical, decreasing steps.
    if max(abs(range_step))<1e-2 && counter>20
        H_o2 = curvature(q);
        delta = H_o2\gradient(q);
        delta(abs(delta)>(sysPara.rmax-sysPara.rmin)/2) = 0;
        candidate = reshape(q-delta,3,[]);
        if all(isfinite(candidate),'all') && isreal(candidate) && ...
                all(sum(candidate(1:2,:).^2,1)<1) && all(candidate(3,:)>0) && ...
                objective(candidate(:)) < objective(q)
            q = candidate(:);
            if norm(delta)<sysPara.threshold_fine, break; end
        end
    end
end
pos_para_fine = reshape(q,3,[]);
pos_fine = [pos_para_fine(3,:).*sqrt(1-sum(pos_para_fine(1:2,:).^2,1)); ...
    pos_para_fine(3,:).*pos_para_fine(1:2,:)];
% Curvature is requested explicitly by the final cal_Hess call, so a Hessian
% calculated before a Newton step is never advertised as the final curvature.
H_o2 = zeros(numel(q),numel(q));
end
