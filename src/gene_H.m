function H = gene_H(input,para)

sub_element_num_x = para.sub_element_num_x;
sub_element_num_y = para.sub_element_num_y;
array_num_x = para.array_num_x;
array_num_y = para.array_num_y;
source_num = para.source_num;
array_n = para.array_n;
pos_element = para.pos_element;
W_vec = para.W_vec;
lambda = para.lambda;

source_pos_vec = input(1:3*source_num);
rho_pseTrue_amp = input(3*source_num+1:4*source_num);
rho_pseTrue_ang = input(4*source_num+1:5*source_num);


source_pos = reshape(source_pos_vec,3,[]);
rho_pseTrue = rho_pseTrue_amp .* exp(1j*rho_pseTrue_ang);

dis_element = zeros(source_num,sub_element_num_x,sub_element_num_y,array_num_x,array_num_y);
rho_element_bias = zeros(size(dis_element));
for sor_idx = 1:source_num
    dis_element(sor_idx,:,:,:,:) = vecnorm(source_pos(:,sor_idx) - pos_element,2,1);
    rho_element_bias(sor_idx,:,:,:,:) = (sqrt(squeeze(pagemtimes(array_n.',source_pos(:,sor_idx) - pos_element))./squeeze(dis_element(sor_idx,:,:,:,:)))).^3;
end
A_part = permute(reshape(dis_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
    floor((array_num_x+1)/2),floor((array_num_y+1)/2))./dis_element, source_num,[]),[2,1]);
B_part = permute(reshape(rho_element_bias./rho_element_bias(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
    floor((array_num_x+1)/2),floor((array_num_y+1)/2)),source_num,[]),[2,1]);
C_part = permute(reshape(exp(-1j*2*pi*dis_element/lambda) ./ exp(-1j*2*pi*dis_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
    floor((array_num_x+1)/2),floor((array_num_y+1)/2))/lambda),source_num,[]),[2,1]);
Array_mat = A_part.*B_part.*C_part;
% Array_mat = C_part;
H_all = W_vec*Array_mat * diag(rho_pseTrue);
% y_compress = y_vec;
H = [real(H_all(:));imag(H_all(:))];
end