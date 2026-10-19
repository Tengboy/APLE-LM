function [alpha_armijo,counter] = armijo_part1(gama,d,f0,thres,alpha_ini)
% Backtracking for direction cosines; retain the original strict-decrease rule.
% Despite its historical name, this is not the full Armijo inequality.
alpha_armijo = 0;
counter = 0;
direction = reshape(d,3,[]);
direction(3,:) = 0;
direction = direction(:);
if ~any(direction), return; end
if alpha_ini == 0, alpha_ini = 0.001; end
f_ori = f0(gama);
alpha = alpha_ini;
while alpha > thres
    candidate = gama(:) + alpha*direction;
    q = reshape(candidate,3,[]);
    counter = counter + 1;
    % Never evaluate the channel at complex coordinates or nonpositive range.
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
