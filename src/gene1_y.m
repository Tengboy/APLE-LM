function [y_vec,H_true,rho_pseTrue] = gene1_y(sysPara)

source_num = sysPara.source_num;
sub_element_num_x = sysPara.sub_element_num_x;
sub_element_num_y = sysPara.sub_element_num_y;
array_num_x = sysPara.array_num_x;
array_num_y = sysPara.array_num_y;
source_pos = sysPara.source_pos;
pos_element = sysPara.pos_element;
array_n = sysPara.array_n;
lambda = sysPara.lambda;
W_vec = sysPara.W_vec;


dis_element = zeros(source_num,sub_element_num_x,sub_element_num_y,array_num_x,array_num_y);
rho_element_bias = zeros(size(dis_element));
for sor_idx = 1:source_num
    dis_element(sor_idx,:,:,:,:) = vecnorm(source_pos(:,sor_idx) - pos_element,2,1);
    rho_element_bias(sor_idx,:,:,:,:) = (sqrt(squeeze(pagemtimes(array_n.',source_pos(:,sor_idx) - pos_element))./squeeze(dis_element(sor_idx,:,:,:,:)))).^3;
end

% 2 real signal
rho_element = 1./dis_element;
% modified rho
A_part = permute(reshape(dis_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
    floor((array_num_x+1)/2),floor((array_num_y+1)/2))./dis_element, source_num,[]),[2,1]);
B_part = permute(reshape(rho_element_bias./rho_element_bias(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
    floor((array_num_x+1)/2),floor((array_num_y+1)/2)),source_num,[]),[2,1]);
C_part = permute(reshape(exp(-1j*2*pi*dis_element/lambda) ./ exp(-1j*2*pi*dis_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
    floor((array_num_x+1)/2),floor((array_num_y+1)/2))/lambda),source_num,[]),[2,1]);
Array_mat = A_part.*B_part.*C_part;
rho_pseTrue = rho_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
    floor((array_num_x+1)/2),floor((array_num_y+1)/2)).*rho_element_bias(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
    floor((array_num_x+1)/2),floor((array_num_y+1)/2)).* exp(-1j*2*pi*dis_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
    floor((array_num_x+1)/2),floor((array_num_y+1)/2))/lambda);
% power control
x = exp(1j*2*pi*rand(source_num,1))./abs(rho_pseTrue);
y_vec = Array_mat * (rho_pseTrue .* x);
% Scale each user column without constructing a dense diagonal matrix.
H_true = (W_vec * Array_mat) .* (rho_pseTrue.*x).';

rho_pseTrue = rho_pseTrue.* x;
end
