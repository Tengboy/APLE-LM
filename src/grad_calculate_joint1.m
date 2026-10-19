function [grad,rho_est,rho_invcov]= grad_calculate_joint1(para,y,W,inv_C,pos_element,lambda)
sor_num = length(para)/3;
para = reshape(para,3,sor_num);
zth = sqrt(1- para(1,:).^2 - para(2,:).^2);
pos = [para(3,:).*zth;para(3,:).*para(1,:);para(3,:).*para(2,:)];
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

temp_vec1 = (steer_vec'*inv_C*steer_vec)\steer_vec'*inv_C*y;

temp_mat1 = zeros(3,3,sor_num);
for sor_idx = 1:sor_num
    temp_mat1(:,:,sor_idx) = [-para(3,sor_idx)*para(1,sor_idx)/sqrt(1- para(1,sor_idx)^2 - para(2,sor_idx)^2),...
        -para(3,sor_idx)*para(2,sor_idx)/sqrt(1- para(1,sor_idx)^2 - para(2,sor_idx)^2),...
        sqrt(1- para(1,sor_idx)^2 - para(2,sor_idx)^2);...
        para(3,sor_idx),0,para(1,sor_idx);...
        0,para(3,sor_idx),para(2,sor_idx)];
end

% temp_mat2 = zeros(RF_chain,3*sor_num,length(temp_vec1));
% counter = 0;
% for sor_idx = 1:sor_num
%     for y_arr_idx = 1:array_num_x
%         for x_arr_idx = 1:array_num_y
%             counter = counter + 1;
%             W_temp = W(:,:,:,x_arr_idx,y_arr_idx);
%             W_temp = reshape(W_temp,RF_chain,[]);
%             temp_mat3 = reshape(pos(:,sor_idx) - pos_element(:,:,:,x_arr_idx,y_arr_idx),3,[]);
%             temp_mat3 = temp_mat3./vecnorm(temp_mat3);
%             temp_mat2(:,:,counter) = temp_vec1(counter) * (-1j*2*pi/lambda)* W_temp*...
%                 (reshape(Y_temp(sor_idx,:,:,x_arr_idx,y_arr_idx),[],1).* (temp_mat3)')...
%                 *[zeros(3,(sor_idx-1)*3),temp_mat1(:,:,sor_idx),zeros(3,(sor_num-sor_idx)*3)];
%         end
%     end
% end
% grad = -2*real(y'*inv_C*sum(temp_mat2,3)- temp_vec1'*steer_vec'*inv_C*sum(temp_mat2,3));
% grad1 = grad.';
% for idx = 1:length(temp_vec1)
%     sor_idx = ceil(idx/array_num_x/array_num_y);
%     pos(:,sor_idx) - pos_element(:,:,x_arr_idx,y_arr_idx)
%     temp_mat2(:,:,idx) = temp_vec1(idx) * (-1j*2*pi/lambda)* steer_vec(:,idx)* (1)'...
%         *[zeros(3,(sor_idx-1)*3),temp_mat1(:,:,sor_idx),zeros(3,(sor_num-sor_idx)*3)];
% end
%% grad2
temp_mat2 = zeros(RF_chain,length(temp_vec1),3*sor_num);
counter = 0;
for sor_idx = 1:sor_num
    for y_arr_idx = 1:array_num_y
        for x_arr_idx = 1:array_num_x
            counter = counter + 1;
            W_temp = W(:,:,x_arr_idx,y_arr_idx);
            % W_temp = reshape(W_temp,RF_chain,[]);
            temp_mat3 = reshape(pos(:,sor_idx) - pos_element(:,:,:,x_arr_idx,y_arr_idx),3,[]);
            temp_mat3 = temp_mat3./vecnorm(temp_mat3);
            temp_mat2(:,counter,:) = (-1j*2*pi/lambda)* W_temp*...
                (reshape(Y_temp(sor_idx,:,:,x_arr_idx,y_arr_idx)...
                ,[],1).* (temp_mat3)')...
                *[zeros(3,(sor_idx-1)*3),temp_mat1(:,:,sor_idx),zeros(3,(sor_num-sor_idx)*3)];
        end
    end
end
grad = real(squeeze(pagemtimes(-2*y'*inv_C+2*temp_vec1'*steer_vec'*inv_C,pagemtimes(temp_mat2,temp_vec1))));
rho_est = temp_vec1;
rho_invcov = (steer_vec'*inv_C*steer_vec);
end