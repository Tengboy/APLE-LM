% theta_range: 扇区范围，单位为角度
% rmin,rmax: 扇区半径范围

function [source_pos,source_pos_para] = gene_pos(sysPara)

source_num = sysPara.source_num;
sub_element_num_x = sysPara.sub_element_num_x;
sub_element_num_y = sysPara.sub_element_num_y;
rmin = sysPara.rmin;
rmax = sysPara.rmax;
theta_range = sysPara.theta_range;


theta_range_rad = theta_range/2/180*pi;

% Direction-cosine spacing follows each subarray axis. The default 15-by-15
% case keeps the original positions and random-number sequence unchanged.
gap_x = 2/sub_element_num_x;
gap_y = 2/sub_element_num_y;

source_pos = zeros(3,source_num);
xAoA_set = [];
yAoA_set = [];
for sor_idx = 1:source_num
    flag = 1;
    while flag
        th1 = theta_range_rad*rand();
        th2 = 2*pi*rand();
        r = rmin + (rmax-rmin)*rand();

        pos_x = r*sin(th1)*cos(th2);
        pos_y = r*sin(th1)*sin(th2);
        pos_z = r*cos(th1);
        pos_temp = [pos_z;pos_x;pos_y];

        xAoA_temp = pos_temp(2,:)./vecnorm(pos_temp);
        yAoA_temp = pos_temp(3,:)./vecnorm(pos_temp);

        if sor_idx == 1
            xAoA_set = [xAoA_set, xAoA_temp];
            yAoA_set = [yAoA_set, yAoA_temp];
            source_pos(:,sor_idx) = pos_temp;
            flag = 0;
        elseif min(abs(xAoA_set-xAoA_temp))>gap_x ...
                || min(abs(yAoA_set-yAoA_temp))>gap_y
        % elseif min(abs(xAoA_set-xAoA_temp))>2/sub_element_num_x ...
        %         && min(abs(yAoA_set-yAoA_temp))>2/sub_element_num_y
            xAoA_set = [xAoA_set, xAoA_temp];
            yAoA_set = [yAoA_set, yAoA_temp];
            source_pos(:,sor_idx) = pos_temp;
            flag = 0;
        end
    end
end
%     source_pos = [10,5]';
source_pos_para = [source_pos(2,:)./vecnorm(source_pos);...
    source_pos(3,:)./vecnorm(source_pos);...
    vecnorm(source_pos)];

end
