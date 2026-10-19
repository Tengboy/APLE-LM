function [index_perm,position_mse] = match_user_positions(estimated,truth)
% One-to-one minimum-squared-distance assignment, independent of toolboxes.
% index_perm(k) identifies the estimated user assigned to true user k.
% Hungarian shortest augmenting path algorithm: O(K^3) time, O(K^2) memory.
K = size(truth,2);
if ~isequal(size(estimated),size(truth)) || size(truth,1)~=3 || K==0 || ...
        ~isreal(estimated) || ~isreal(truth) || ...
        any(~isfinite(estimated),'all') || any(~isfinite(truth),'all')
    error('APLE:InvalidPositions','Matching requires two finite real 3-by-K position arrays.');
end
cost = zeros(K);
for user = 1:K
    cost(user,:) = sum((estimated-truth(:,user)).^2,1);
end
row_potential = zeros(K,1);
column_potential = zeros(K+1,1);
row_for_column = zeros(K+1,1); % Index 1 is the temporary unmatched column.
predecessor = zeros(K+1,1);
for new_row = 1:K
    row_for_column(1) = new_row;
    column = 1;
    distance = inf(K+1,1);
    visited = false(K+1,1);
    while true
        visited(column) = true;
        row = row_for_column(column);
        delta = inf;
        next_column = 0;
        for candidate = 2:K+1
            if ~visited(candidate)
                reduced = cost(row,candidate-1)-row_potential(row)-column_potential(candidate);
                if reduced < distance(candidate)
                    distance(candidate) = reduced;
                    predecessor(candidate) = column;
                end
                if distance(candidate) < delta
                    delta = distance(candidate);
                    next_column = candidate;
                end
            end
        end
        for candidate = 1:K+1
            if visited(candidate)
                row_potential(row_for_column(candidate)) = row_potential(row_for_column(candidate))+delta;
                column_potential(candidate) = column_potential(candidate)-delta;
            else
                distance(candidate) = distance(candidate)-delta;
            end
        end
        column = next_column;
        if row_for_column(column)==0, break; end
    end
    while column~=1
        previous = predecessor(column);
        row_for_column(column) = row_for_column(previous);
        column = previous;
    end
end
index_perm = zeros(1,K);
for column = 2:K+1
    index_perm(row_for_column(column)) = column-1;
end
position_mse = sum(cost(sub2ind([K,K],1:K,index_perm)))/K;
end
