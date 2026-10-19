function [FIM_partial_mat,Transition_mat] = cal_FIM(rho_pseTrue,sysPara)

para.sub_element_num_x = sysPara.sub_element_num_x;
para.sub_element_num_y = sysPara.sub_element_num_y;
para.array_num_x = sysPara.array_num_x;
para.array_num_y = sysPara.array_num_y;
para.source_num = sysPara.source_num;
para.array_n = sysPara.array_n;
para.pos_element = sysPara.pos_element;
para.W_vec = sysPara.W_vec;
para.lambda = sysPara.lambda;
variable_set = [sysPara.source_pos(:);abs(rho_pseTrue);angle(rho_pseTrue)];

receive_fun = @(x)gene_y(x,para);
FIM_partial_mat = grad_numerical_vec(receive_fun,variable_set,1e-7);

%% Transition matrix to calculate CRLB of H 
% main_240803 only requests the received-signal Jacobian.
if nargout > 1
    H_fun = @(x)gene_H(x,para);
    Transition_mat = grad_numerical_vec(H_fun,variable_set,1e-7);
end
end
