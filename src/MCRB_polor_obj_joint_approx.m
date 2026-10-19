% L 子阵列个数
% W 波束赋形矩阵
% y 接收信号
% M 射频链个数
function out = MCRB_polor_obj_joint_approx(pos_para,para,montPara)
inv_noise_cov = montPara.inv_C_n;
y_res = montPara.y_compress_noise;
W = para.W;
pos_element = para.pos_element;
lambda = para.lambda;
%%
sor_num = length(pos_para)/3;
pos_para = reshape(pos_para,3,sor_num);
zth = sqrt(1- pos_para(1,:).^2 - pos_para(2,:).^2);
pos = [pos_para(3,:).*zth;pos_para(3,:).*pos_para(1,:);pos_para(3,:).*pos_para(2,:)];
[RF_chain,~,array_num_x,array_num_y] = size(W);
[~,sub_element_num_x,sub_element_num_y,~,~] = size(pos_element);
dis_element_grid = zeros(sor_num,sub_element_num_x,sub_element_num_y,array_num_x,array_num_y);
for sor_idx = 1:sor_num
    dis_element_grid(sor_idx,:,:,:,:) = vecnorm(pos(:,sor_idx) - pos_element,2,1);
end
Y_temp = zeros(sor_num,sub_element_num_x,sub_element_num_y,array_num_x,array_num_y);
for sor_idx = 1:sor_num
    Y_temp(sor_idx,:,:,:,:) = squeeze(exp(-1j*2*pi*dis_element_grid(sor_idx,:,:,:,:)/lambda));
end

steer_mat = zeros(RF_chain,array_num_x,array_num_y,sor_num);

for sor_idx = 1:sor_num
    for x_arr_idx = 1:array_num_x
        for y_arr_idx = 1:array_num_y
            W_temp = W(:,:,x_arr_idx,y_arr_idx);
            % W_temp = reshape(W_temp,RF_chain,[]);
            steer_mat(:,x_arr_idx,y_arr_idx,sor_idx) = W_temp*reshape(Y_temp(sor_idx,:,:,x_arr_idx,y_arr_idx)...
                ,[],1);
        end
    end
end
steer_vec = reshape(steer_mat,RF_chain,[]);
rho_temp = (steer_vec'*(inv_noise_cov)*steer_vec)\steer_vec'*(inv_noise_cov)*y_res;
temp = y_res-steer_vec*rho_temp;
%             temp = y_compress-Mat_B/(Mat_B'*Mat_B)*Mat_B'*y_compress;
out = abs(temp' *(inv_noise_cov)* temp);
                
end