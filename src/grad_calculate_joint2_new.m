function [grad,rho_est,rho_invcov]= grad_calculate_joint2_new(para,y,inv_C,pos_element,lambda)
sor_num = length(para)/3;
para = reshape(para,3,sor_num);
zth = sqrt(1- para(1,:).^2 - para(2,:).^2);
pos = [para(3,:).*zth;para(3,:).*para(1,:);para(3,:).*para(2,:)];
[~,sub_element_num_x,sub_element_num_y,array_num_x,array_num_y] = size(pos_element);

dis_element_grid = zeros(sor_num,array_num_x,array_num_y);
for sor_idx = 1:sor_num
    dis_element_grid(sor_idx,:,:) = vecnorm(pos(:,sor_idx) - pos_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),:,:),2,1);
end
Y_temp = zeros(sor_num,array_num_x,array_num_y);
for sor_idx = 1:sor_num
    Y_temp(sor_idx,:,:) = (exp(-1j*2*pi*dis_element_grid(sor_idx,:,:)/lambda));
    %     Y_temp(sor_idx,:,:,:,:) = squeeze(exp(-1j*2*pi*dis_element_grid(sor_idx,:,:,:,:)/lambda));
end

steer_mat = [];

for sor_idx = 1:sor_num
    steer_mat = blkdiag(steer_mat,reshape(Y_temp(sor_idx,:,:)./Y_temp(sor_idx,floor((array_num_x+1)/2),floor((array_num_y+1)/2)),[],1));
end

temp_vec1 = (steer_mat'*inv_C*steer_mat)\steer_mat'*inv_C*y;

temp_mat1 = zeros(3,3,sor_num);
for sor_idx = 1:sor_num
    temp_mat1(:,:,sor_idx) = [-para(3,sor_idx)*para(1,sor_idx)/sqrt(1- para(1,sor_idx)^2 - para(2,sor_idx)^2),...
        -para(3,sor_idx)*para(2,sor_idx)/sqrt(1- para(1,sor_idx)^2 - para(2,sor_idx)^2),...
        sqrt(1- para(1,sor_idx)^2 - para(2,sor_idx)^2);...
        para(3,sor_idx),0,para(1,sor_idx);...
        0,para(3,sor_idx),para(2,sor_idx)];
end

%% grad2
temp_mat2 = zeros(length(y),length(temp_vec1),3*sor_num);
counter = 0;
C_part = permute(reshape(exp(-1j*2*pi*dis_element_grid/lambda) ./ exp(-1j*2*pi*dis_element_grid(:,floor((array_num_x+1)/2),floor((array_num_y+1)/2))/lambda),sor_num,[]),[2,1]);
for sor_idx = 1:sor_num
    counter = counter + 1;
    %             temp_mat3 = reshape(pos(:,sor_idx) - pos_element(:,:,:,x_arr_idx,y_arr_idx),3,[]);
    %             temp_mat3 = temp_mat3./vecnorm(temp_mat3);
    % 由于导向矢量的重构，此处temp_mat2需要重新写
    e_dir_mat = reshape(pos(:,sor_idx) - pos_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),:,:),3,[]);
    e_dir_mat = e_dir_mat./reshape(dis_element_grid(sor_idx,:,:),1,[]);
    e_dir_ref = pos(:,sor_idx) - squeeze(pos_element(:,floor((sub_element_num_x+1)/2),floor((sub_element_num_y+1)/2),...
        floor((array_num_x+1)/2),floor((array_num_y+1)/2)));
    e_dir_ref = e_dir_ref/norm(e_dir_ref);
    % 对位置求导 可分为三部分

    C_part_p = -1j*2*pi/lambda*diag(C_part(:,sor_idx))*(e_dir_mat - e_dir_ref).';
    temp_mat2((sor_idx-1)*array_num_x*array_num_y+1:sor_idx*array_num_x*array_num_y,sor_idx,:) = (C_part_p)*...
        [zeros(3,(sor_idx-1)*3),temp_mat1(:,:,sor_idx),zeros(3,(sor_num-sor_idx)*3)];
end
grad = real(squeeze(pagemtimes(-2*y'*inv_C+2*temp_vec1'*steer_mat'*inv_C,pagemtimes(temp_mat2,temp_vec1))));
rho_est = temp_vec1;
rho_invcov = (steer_mat'*inv_C*steer_mat);
end