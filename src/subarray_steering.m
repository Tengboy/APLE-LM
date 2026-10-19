function B = subarray_steering(positions,sysPara,phase_only)
% B in paper Eq. (18): one column per subarray, then per user.
% Code coordinates are [array normal; first array axis; second array axis].
% Each subarray response is normalized by its own reference antenna.
% Optional phase-only response implements the ACM model in paper Eq. (8).
if nargin<3, phase_only=false; end
R = sysPara.RF_chain;
nx = sysPara.sub_element_num_x;
ny = sysPara.sub_element_num_y;
S = sysPara.array_num_x*sysPara.array_num_y;
N = nx*ny;
K = size(positions,2);
reference = floor((nx+1)/2) + nx*(floor((ny+1)/2)-1);
elements = reshape(sysPara.pos_element,3,N*S);
W = reshape(sysPara.W,R,N,S);
B = complex(zeros(R,S*K));
for user = 1:K
    offset = positions(:,user)-elements;
    distance = vecnorm(offset,2,1);
    if phase_only
        response = exp(-1j*2*pi*distance/sysPara.lambda);
    else
        radiation = sqrt((sysPara.array_n.'*offset)./distance).^3;
        response = radiation./distance.*exp(-1j*2*pi*distance/sysPara.lambda);
    end
    response = reshape(response,N,S);
    response = response./response(reference,:);
    B(:,(user-1)*S+(1:S)) = reshape( ...
        pagemtimes(W,reshape(response,N,1,S)),R,S);
end
end
