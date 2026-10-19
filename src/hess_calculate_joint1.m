function [hess,rho_est] = hess_calculate_joint1(para,y,W,inv_C,pos_element,lambda)
% Full parameter curvature of the phase-only coarse profile objective.
gradient=@(q)grad_calculate_joint1(q,y,W,inv_C,pos_element,lambda);
hess=hessian_from_gradient(gradient,para);
if nargout>1, [~,rho_est]=gradient(para); end
end
