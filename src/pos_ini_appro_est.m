function [pos_ini,pos_para_ini] = pos_ini_appro_est(montPara,sysPara)
% ACM counterpart: phase-only channel; same verified update schedule and state fixes.
% Sequential grid initialization followed by angle/range refinement.
% No true user positions are needed; only the known number of users is used.
K = sysPara.source_num;
pos_ini = zeros(3,K);
pos_para_ini = zeros(3,K);
y = montPara.y_compress_noise;
Q = montPara.inv_C_n;
for user = 1:K
    B = subarray_steering(pos_ini(:,1:user-1),sysPara,true);
    y_res = y;
    if ~isempty(B)
        weighted_B = B'*Q;
        gain = (weighted_B*B)\(weighted_B*y);
        y_res = y-B*gain;
    end

    % Evaluate all stored grid projection scores in one page operation.
    projections = pagemtimes(sysPara.steer_vec_noise_inv2,y_res);
    scores = abs(sum(conj(projections).* ...
        pagemtimes(sysPara.steer_vec_noise_inv,projections),1));
    [~,grid_index] = max(scores(:));
    position = sysPara.pos_grid(:,grid_index);
    q = [position(2:3)/norm(position);norm(position)];
    previous_position = position;
    local_mont = montPara;
    local_mont.y_compress_noise = y_res;
    objective = @(x)MCRB_polor_obj_joint_approx(x,sysPara,local_mont);
    alpha_angle = 0.01;
    alpha_range = 0.1;

    % Preserve the original alternating schedule and stopping criterion.
    for counter = 2:sysPara.counter_ini_max+1
        gradient = grad_calculate_joint1(q,y_res,sysPara.W,Q, ...
            sysPara.pos_element,sysPara.lambda);
        direction = sign(gradient);
        if mod(counter,2)==0
            alpha_angle = armijo_part1(q,-direction,objective,1e-6,2*alpha_angle);
            q(1:2) = q(1:2)-alpha_angle*direction(1:2);
        else
            alpha_range = armijo_part2(q,-direction,objective,1e-6,2*alpha_range);
            q(3) = min(max(q(3)-alpha_range*direction(3),sysPara.rmin),sysPara.rmax);
            position = q(3)*[sqrt(1-q(1)^2-q(2)^2);q(1:2)];
            if norm(position-previous_position)/norm(previous_position)<sysPara.threshold_ini
                break;
            end
            previous_position = position;
        end
    end
    % Also synchronize on an angular last iteration (odd iteration limit).
    pos_para_ini(:,user) = q;
    pos_ini(:,user) = q(3)*[sqrt(1-q(1)^2-q(2)^2);q(1:2)];
end
end
