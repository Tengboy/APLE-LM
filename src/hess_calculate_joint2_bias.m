function [hess,rho_est] = hess_calculate_joint2_bias(para,y,inv_C,pos_element,lambda,array_n)
% Curvature of the fine objective with the same cos(theta)^(3/2) radiation
% model as its objective and gradient. The original Hessian used exponent 1/2.
gradient = @(q)grad_calculate_joint2_bias(q,y,inv_C,pos_element,lambda,array_n);
hess = hessian_from_gradient(gradient,para);
if nargout > 1
    [~,rho_est] = gradient(para);
end
end
