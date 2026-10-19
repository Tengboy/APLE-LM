function [grad,rho_est,rho_invcov]= grad_calculate_joint1_bias(para,y,W,inv_C,pos_element,lambda,array_n)
sor_num = length(para)/3;
para = reshape(para,3,sor_num);
zth = sqrt(1- para(1,:).^2 - para(2,:).^2);
pos = [para(3,:).*zth;para(3,:).*para(1,:);para(3,:).*para(2,:)];
[RF_chain,~,array_num_x,array_num_y] = size(W);
[~,sub_element_num_x,sub_element_num_y,~,~] = size(pos_element);
dis_element_grid = zeros(sor_num,sub_element_num_x,sub_element_num_y,array_num_x,array_num_y);
rho_element_bias_grid = zeros(size(dis_element_grid));
for sor_idx = 1:sor_num
    dis_element_grid(sor_idx,:,:,:,:) = vecnorm(pos(:,sor_idx) - pos_element,2,1);
    rho_element_bias_grid(sor_idx,:,:,:,:) = (sqrt(squeeze(pagemtimes(array_n.',pos(:,sor_idx) - pos_element))./squeeze(dis_element_grid(sor_idx,:,:,:,:)))).^3;
end
Y_temp = zeros(sor_num,sub_element_num_x,sub_element_num_y,array_num_x,array_num_y);
for sor_idx = 1:sor_num
    Y_temp(sor_idx,:,:,:,:) = 1./dis_element_grid(sor_idx,:,:,:,:).*rho_element_bias_grid(sor_idx,:,:,:,:)...
        .*(exp(-1j*2*pi*dis_element_grid(sor_idx,:,:,:,:)/lambda));
%     Y_temp(sor_idx,:,:,:,:) = squeeze(exp(-1j*2*pi*dis_element_grid(sor_idx,:,:,:,:)/lambda));
end

steer_mat = zeros(RF_chain,sor_num,array_num_x,array_num_y);

for sor_idx = 1:sor_num
    for y_arr_idx = 1:array_num_y
        for x_arr_idx = 1:array_num_x
            W_temp = (W(:,:,x_arr_idx,y_arr_idx));
            % W_temp = reshape(W_temp,RF_chain,[]);
            steer_mat(:,sor_idx,x_arr_idx,y_arr_idx) = W_temp*reshape(Y_temp(sor_idx,:,:,x_arr_idx,y_arr_idx)...
                /Y_temp(sor_idx,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),x_arr_idx,y_arr_idx),[],1);
        end
    end
end
steer_vec = reshape(steer_mat,RF_chain,[]);

temp_vec1 = (steer_vec'*inv_C*steer_vec)\steer_vec'*inv_C*y;    % q_1

temp_mat1 = zeros(3,3,sor_num);
for sor_idx = 1:sor_num
    temp_mat1(:,:,sor_idx) = [-para(3,sor_idx)*para(1,sor_idx)/sqrt(1- para(1,sor_idx)^2 - para(2,sor_idx)^2),...
        -para(3,sor_idx)*para(2,sor_idx)/sqrt(1- para(1,sor_idx)^2 - para(2,sor_idx)^2),...
        sqrt(1- para(1,sor_idx)^2 - para(2,sor_idx)^2);...
        para(3,sor_idx),0,para(1,sor_idx);...
        0,para(3,sor_idx),para(2,sor_idx)];
end

%% grad2
temp_mat2 = zeros(RF_chain,length(temp_vec1),3*sor_num);
counter = 0;
for y_arr_idx = 1:array_num_y
    for x_arr_idx = 1:array_num_x
        A_part = permute(reshape(dis_element_grid(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
            x_arr_idx,y_arr_idx)./dis_element_grid(:,:,:,x_arr_idx,y_arr_idx), sor_num,[]),[2,1]);
        B_part = permute(reshape(rho_element_bias_grid(:,:,:,x_arr_idx,y_arr_idx)./rho_element_bias_grid(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
            x_arr_idx,y_arr_idx),sor_num,[]),[2,1]);
        C_part = permute(reshape(exp(-1j*2*pi*dis_element_grid(:,:,:,x_arr_idx,y_arr_idx)/lambda) ./ exp(-1j*2*pi*dis_element_grid(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
            x_arr_idx,y_arr_idx)/lambda),sor_num,[]),[2,1]);
        for sor_idx = 1:sor_num
            counter = counter + 1;
            W_temp = (W(:,:,x_arr_idx,y_arr_idx));
            reshape_dis_element_grid = reshape(dis_element_grid(sor_idx,:,:,x_arr_idx,y_arr_idx),1,[]);
            % W_temp = reshape(W_temp,RF_chain,[]);
%             temp_mat3 = reshape(pos(:,sor_idx) - pos_element(:,:,:,x_arr_idx,y_arr_idx),3,[]);
%             temp_mat3 = temp_mat3./vecnorm(temp_mat3);
            % 由于导向矢量的重构，此处temp_mat2需要重新写
            e_dir_mat = reshape(pos(:,sor_idx) - pos_element(:,:,:,x_arr_idx,y_arr_idx),3,[]);
            e_dir_mat = e_dir_mat./reshape_dis_element_grid;
            e_dir_ref = pos(:,sor_idx) - squeeze(pos_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
                x_arr_idx,y_arr_idx));
            e_dir_ref = e_dir_ref/norm(e_dir_ref);
            % 对位置求导 可分为三部分
            A_part_p = (1./reshape_dis_element_grid.').^2 .* (reshape_dis_element_grid.'*e_dir_ref.' - ...
                dis_element_grid(sor_idx,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
                x_arr_idx,y_arr_idx)*e_dir_mat.');
            B_part_p = permute(...
                ((array_n - squeeze(pagemtimes(pagemtimes(reshape(e_dir_mat,3,1,[]),reshape(e_dir_mat,1,3,[])),array_n)))./...
                    reshape_dis_element_grid * (e_dir_ref.'* array_n)- ....
                (array_n - e_dir_ref*e_dir_ref.'*array_n)./dis_element_grid(sor_idx,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
                x_arr_idx,y_arr_idx)*(e_dir_mat.'*array_n).')/(e_dir_ref.'* array_n)^2/2 *3.* (B_part(:,sor_idx).^(1/3)).'...
                ,[2,1]);
            C_part_p = -1j*2*pi/lambda*C_part(:,sor_idx).*(e_dir_mat - e_dir_ref).';
            temp_mat2(:,counter,:) = W_temp * (B_part(:,sor_idx) .* C_part(:,sor_idx).* A_part_p +...
                A_part(:,sor_idx).* C_part(:,sor_idx) .* B_part_p  + A_part(:,sor_idx) .* B_part(:,sor_idx) .* C_part_p)*...
                [zeros(3,(sor_idx-1)*3),temp_mat1(:,:,sor_idx),zeros(3,(sor_num-sor_idx)*3)];
        end
    end
end
grad = real(squeeze(pagemtimes(-2*y'*inv_C+2*temp_vec1'*steer_vec'*inv_C,pagemtimes(temp_mat2,temp_vec1))));
rho_est = temp_vec1;
rho_invcov = (steer_vec'*inv_C*steer_vec);
end