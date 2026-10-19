function H = hessian_from_gradient(gradient,parameters)
% Central differences of an analytical gradient, with relative steps.
% Used for final message curvature (and the existing optional Newton step),
% not on every gradient iteration. No eigenvalue clipping or damping is added.
parameters = parameters(:);
q = reshape(parameters,3,[]);
if ~isreal(q) || any(~isfinite(q),'all') || ...
        any(sum(q(1:2,:).^2,1)>=1) || any(q(3,:)<=0)
    error('APLE:InvalidParameters','Hessian requires real, physical position parameters.');
end
n = numel(parameters);
H = zeros(n,n);
for k = 1:n
    step = eps^(1/3)*max(1,abs(parameters(k)));
    user = ceil(k/3);
    component = mod(k-1,3)+1;
    if component == 3
        step = min(step,parameters(k)/4);
    else
        % Both perturbations must stay strictly inside the direction-cosine disk.
        other = 3-component;
        margin = sqrt(1-q(other,user)^2)-abs(parameters(k));
        step = min(step,margin/4);
    end
    if step <= eps(max(1,abs(parameters(k))))
        error('APLE:HessianBoundary','Position is too close to the physical boundary for central differences.');
    end
    delta = zeros(n,1);
    delta(k) = step;
    H(:,k) = (gradient(parameters+delta)-gradient(parameters-delta))/(2*step);
end
H = real((H+H.')/2);
end
