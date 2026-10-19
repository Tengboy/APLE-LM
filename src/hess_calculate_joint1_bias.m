function [hess,rho_est] = hess_calculate_joint1_bias(para,y,W,inv_C,pos_element,lambda,array_n)
% Curvature of the profiled coarse objective in [chi_x,chi_y,r] coordinates.
% Differentiate the verified analytical gradient. This also includes the
% second-order chain rule omitted by the original Cartesian Hessian pullback.
gradient = @(q)grad_calculate_joint1_bias(q,y,W,inv_C,pos_element,lambda,array_n);
hess = hessian_from_gradient(gradient,para);
if nargout > 1
    [~,rho_est] = gradient(para);
end
end
