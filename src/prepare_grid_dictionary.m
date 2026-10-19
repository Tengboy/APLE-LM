function [steer_vec,gram_inverse,weighted_steering,noise_precision] = prepare_grid_dictionary(pos_grid,W,pos_element,lambda,array_n,phase_only)
% Cache the same grid projections as the original main program.
% W is fixed, so factor W*W' once rather than twice at every grid point.
if nargin<6, phase_only=false; end
[R,N,ax,ay] = size(W);
S = ax*ay;
G = size(pos_grid,2);
W_vec = reshape(W,R,[]);
noise_chol = chol(W_vec*W_vec','lower');
W_pages = reshape(W,R,N,S);
steer_vec = complex(zeros(R,S,G));
gram_inverse = complex(zeros(S,S,G));
weighted_steering = complex(zeros(S,R,G));
for grid = 1:G
    distance = squeeze(vecnorm(pos_grid(:,grid)-pos_element,2,1));
    if phase_only
        response = exp(-1j*2*pi*distance/lambda);
    else
        radiation = sqrt(squeeze(pagemtimes(array_n.',pos_grid(:,grid)-pos_element))./distance).^3;
        response = radiation./distance.*exp(-1j*2*pi*distance/lambda);
    end
    B = reshape(pagemtimes(W_pages,reshape(response,N,1,S)),R,S);
    BtQ = (noise_chol'\(noise_chol\B))';
    steer_vec(:,:,grid) = B;
    weighted_steering(:,:,grid) = BtQ;
    gram_inverse(:,:,grid) = (BtQ*B)\eye(S);
end
if nargout > 3
    noise_precision = noise_chol'\(noise_chol\eye(R));
end
end
