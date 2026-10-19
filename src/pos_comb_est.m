function [pos_comb,pos_para_comb]= pos_comb_est(iterPara)
Hess_para_coarse = iterPara.Hess_para_coarse;
Hess_para_fine = iterPara.Hess_para_fine;
pos_para_coarse = iterPara.pos_para_coarse;
pos_para_fine = iterPara.pos_para_fine;
% Deliberate Gaussian fusion in [chi_x,chi_y,r]; see report Q2 for its approximation.
inv_cov_comb_iter1 = Hess_para_coarse + Hess_para_fine;
pos_para_comb = inv_cov_comb_iter1\(Hess_para_coarse*pos_para_coarse(:) + Hess_para_fine*pos_para_fine(:));
pos_para_comb = reshape(pos_para_comb,3,[]);
temp_iter1 = sqrt(1 - pos_para_comb(1,:).^2 - pos_para_comb(2,:).^2);
pos_comb = [pos_para_comb(3,:).*temp_iter1;...
    pos_para_comb(3,:).*pos_para_comb(1,:);...
    pos_para_comb(3,:).*pos_para_comb(2,:);];

end
