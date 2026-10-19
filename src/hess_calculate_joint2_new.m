function [hess,rho_est] = hess_calculate_joint2_new(para,y,inv_C,pos_element,lambda)
% Full parameter curvature of the phase-only fine profile objective.
gradient=@(q)grad_calculate_joint2_new(q,y,inv_C,pos_element,lambda);
hess=hessian_from_gradient(gradient,para);
if nargout>1, [~,rho_est]=gradient(para); end
end
