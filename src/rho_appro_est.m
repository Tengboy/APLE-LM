function [rho_est,rho_invcov] = rho_appro_est(iterPara,montPara,sysPara)
% ACM counterpart: phase-only channel; same verified update schedule and state fixes.
% Weighted least squares for subarray gains, paper Eq. (39).
% Keep user-major column order, consistent with the fine-stage block diagonal.
B = subarray_steering(iterPara.pos_damping,sysPara,true);
weighted_B = B'*montPara.inv_C_n;
rho_invcov = weighted_B*B;
rho_est = rho_invcov\(weighted_B*montPara.y_compress_noise);
end
