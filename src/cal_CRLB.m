function [pos_bound,inv_C_n] = cal_CRLB(y_vec,SNR,FIM_partial_mat,sysPara)
% Conditional deterministic CRLB; no finite Bayesian prior is included.
% Keep the historical API: pos_bound is SUMMED over all users. Divide by K
% when comparing with the per-user position MSE (see main_240803 summary).
signal_pow = sum(abs(y_vec).^2,'all')/numel(y_vec);
noise_pow = signal_pow/10^(SNR/10);
if isfield(sysPara,'noise_precision')
    % Valid only for the fixed W used to prepare the dictionary.
    inv_C_n = sysPara.noise_precision/noise_pow;
else
    inv_C_n = sysPara.W_vec_2\eye(size(sysPara.W_vec_2))/noise_pow;
end
FIM_mat = 2*real(FIM_partial_mat'*inv_C_n*FIM_partial_mat);
% Only position diagonal entries are needed; solve for those columns.
position_count = 3*sysPara.source_num;
selection = eye(size(FIM_mat,1),position_count);
position_covariance_columns = FIM_mat\selection;
pos_bound = sum(diag(position_covariance_columns(1:position_count,:)));
end
