function [pos_bound,H_bound,inv_C_n,H_bound_per_user] = cal_CRLB_n(y_vec,SNR,FIM_partial_mat,Transition_mat,sysPara)
% Conditional CRLB for position and the compressed equivalent channel W*h.
% First three outputs preserve the original API and TOTAL-variance convention.
% Optional fourth output is a 1-by-K vector for mean per-user NMSE bounds.
noise_pow = mean(abs(y_vec(:)).^2)/10^(SNR/10);
if isfield(sysPara,'noise_precision')
    inv_C_n = sysPara.noise_precision/noise_pow;
else
    inv_C_n = sysPara.W_vec_2\eye(size(sysPara.W_vec_2))/noise_pow;
end
FIM_mat = 2*real(FIM_partial_mat'*inv_C_n*FIM_partial_mat);
covariance = FIM_mat\eye(size(FIM_mat));
position_count = 3*sysPara.source_num;
pos_bound = sum(diag(covariance(1:position_count,1:position_count)));
% Only diag(T*C*T') is needed; avoid a (2*RF*K)-by-(2*RF*K) product.
channel_variance = real(sum((Transition_mat*covariance).*conj(Transition_mat),2));
H_bound = sum(channel_variance);
if nargout > 3
    % gene_H stacks real(H(:)), then imag(H(:)); each user occupies RF rows.
    channel_variance = reshape(channel_variance,sysPara.RF_chain,sysPara.source_num,2);
    H_bound_per_user = reshape(sum(channel_variance,[1,3]),1,sysPara.source_num);
end
end
