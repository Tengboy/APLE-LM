function [compressed,physical]=fig6_channel_basis(positions,sys,approximate)
% One steering column per user, normalized at the same reference as gene_H.
elements=reshape(sys.pos_element,3,[]);
shape=[sys.sub_element_num_x,sys.sub_element_num_y,sys.array_num_x,sys.array_num_y];
centre=floor((shape+1)/2);
reference=sub2ind(shape,centre(1),centre(2),centre(3),centre(4));
K=size(positions,2); physical=complex(zeros(size(elements,2),K));
assert(size(positions,1)==3 && isreal(positions) && all(isfinite(positions),'all'));
for user=1:K
    displacement=positions(:,user)-elements;
    distance=sqrt(sum(displacement.^2,1));
    assert(all(distance>0),'F6:InvalidDistance','A user coincides with an antenna.');
    phase=exp(-1j*2*pi*distance/sys.lambda)/exp(-1j*2*pi*distance(reference)/sys.lambda);
    if approximate
        column=phase; % ACM omits distance attenuation and element gain variation.
    else
        cosine=(sys.array_n.'*displacement)./distance;
        assert(all(cosine>0),'F6:InvalidDirection','Accurate model requires the front half-space.');
        element_gain=sqrt(cosine).^3;
        column=(distance(reference)./distance).*(element_gain/element_gain(reference)).*phase;
    end
    physical(:,user)=column.';
end
compressed=sys.W_vec*physical;
end
