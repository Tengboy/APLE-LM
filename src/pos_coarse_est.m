function [pos_coarse,pos_para_coarse,Hessian,alpha_saver] = pos_coarse_est(iterPara,montPara,sysPara)
% Per-user coarse update used by main_240803 (not a joint paper reimplementation).
K = sysPara.source_num;
S = sysPara.array_num_x*sysPara.array_num_y;
pos_coarse = iterPara.pos_coarse;
pos_para_coarse = iterPara.pos_para_coarse;
if isfield(sysPara,'counter_coarse_max')
    maximum = sysPara.counter_coarse_max;
elseif iterPara.turbo_iter_counter == 1
    maximum = sysPara.counter_coarse_max1;
else
    maximum = sysPara.counter_coarse_max2;
end
Hessian = zeros(3*K);
alpha_saver = zeros(maximum,2,K);
Q = montPara.inv_C_n;
y = montPara.y_compress_noise;

% Damped positions are fixed throughout this call: build their responses once.
B_all = subarray_steering(iterPara.pos_damping,sysPara);
for user = 1:K
    other_columns = true(1,S*K);
    other_columns((user-1)*S+(1:S)) = false;
    B = B_all(:,other_columns);
    weighted_B = B'*Q;
    gain = (weighted_B*B)\(weighted_B*y);
    y_res = y-B*gain;
    q = pos_para_coarse(:,user);

    if iterPara.cal_Hess == 1
        % Main discards the returned means in this call. Do not optimize again.
        indices = (user-1)*3+(1:3);
        Hessian(indices,indices) = hess_calculate_joint1_bias(q,y_res, ...
            sysPara.W,Q,sysPara.pos_element,sysPara.lambda,sysPara.array_n);
        continue;
    end

    local_mont = montPara;
    local_mont.y_compress_noise = y_res;
    objective = @(x)MCRB_polor_obj_joint_bias_new(x,sysPara,local_mont);
    alpha_angle = 0.01;
    alpha_range = 0.1;
    angle_count = 0;
    range_count = 0;
    % Keep the legacy schedule: one angle update per three counter values;
    % the stopping test is applied on range updates, including the final one.
    for counter = 2:maximum+2
        gradient = grad_calculate_joint1_bias(q,y_res,sysPara.W,Q, ...
            sysPara.pos_element,sysPara.lambda,sysPara.array_n);
        direction = sign(gradient);
        if mod(counter,3)==0
            alpha_angle = armijo_part1(q,-direction,objective,1e-6,2*alpha_angle);
            q(1:2) = q(1:2)-alpha_angle*direction(1:2);
            angle_count = angle_count+1;
            alpha_saver(angle_count,1,user) = alpha_angle;
        else
            alpha_range = armijo_part2(q,-direction,objective,1e-6,2*alpha_range);
            range_step = alpha_range*direction(3);
            q(3) = q(3)-range_step;
            range_count = range_count+1;
            alpha_saver(range_count,2,user) = alpha_range;
        end
        position = q(3)*[sqrt(1-q(1)^2-q(2)^2);q(1:2)];
        if mod(counter,3)~=0 && ...
                (abs(range_step)/norm(position)<sysPara.threshold_coarse || counter>maximum)
            break;
        end
    end
    pos_para_coarse(:,user) = q;
    pos_coarse(:,user) = position;
end
end
