% L 子阵列个数
% W 波束赋形矩阵
% y 接收信号
% M 射频链个数
function out = MCRB_polor_obj_joint2_approx(pos_para,iterPara,para)
inv_noise_cov = iterPara.rho_invcov;
y_res = iterPara.rho_est;
pos_element = para.pos_element;
lambda = para.lambda;
%%
sor_num = length(pos_para)/3;
pos_para = reshape(pos_para,3,sor_num);
zth = sqrt(1- pos_para(1,:).^2 - pos_para(2,:).^2);
pos = [pos_para(3,:).*zth;pos_para(3,:).*pos_para(1,:);pos_para(3,:).*pos_para(2,:)];
[~,sub_element_num_x,sub_element_num_y,array_num_x,array_num_y] = size(pos_element);

dis_element_grid = zeros(sor_num,array_num_x,array_num_y);
for sor_idx = 1:sor_num
    dis_element_grid(sor_idx,:,:) = squeeze(vecnorm(pos(:,sor_idx) - pos_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),:,:),2,1));
end
Y_temp = zeros(sor_num,array_num_x,array_num_y);
steer_mat = [];

for sor_idx = 1:sor_num
    Y_temp(sor_idx,:,:) = squeeze(exp(-1j*2*pi*dis_element_grid(sor_idx,:,:)/lambda));
    steer_mat = blkdiag(steer_mat,reshape(Y_temp(sor_idx,:,:),[],1));
end

rho_temp = (steer_mat'*(inv_noise_cov)*steer_mat)\steer_mat'*(inv_noise_cov)*y_res;
temp = y_res-steer_mat*rho_temp;
%             temp = y_compress-Mat_B/(Mat_B'*Mat_B)*Mat_B'*y_compress;
out = abs(temp' *(inv_noise_cov)* temp);
                
end