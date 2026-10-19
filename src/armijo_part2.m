function [alpha_armijo,counter] = armijo_part2(gama,d,f0,thres,alpha_ini)
% Backtracking for range; retain the original strict-decrease rule.
alpha_armijo = 0;
counter = 0;
direction = reshape(d,3,[]);
direction(1:2,:) = 0;
direction = direction(:);
if ~any(direction), return; end
f_ori = f0(gama);
alpha = alpha_ini;
while alpha > thres
    candidate = gama(:) + alpha*direction;
    q = reshape(candidate,3,[]);
    counter = counter + 1;
    if all(sum(q(1:2,:).^2,1)<1) && all(q(3,:)>0)
        value = f0(candidate);
        if isreal(value) && isfinite(value) && value < f_ori
            alpha_armijo = alpha;
            return;
        end
    end
    alpha = 0.5*alpha;
end
end
