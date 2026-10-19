function baseline=fig5_full_baseline(accurate)
% View the SAME physical array and W as one full array for ES-GA/P-SIGW.
% Original column order is [ix,iy,subarray_x,subarray_y]. Canonical full-array
% order is [ix,subarray_x,iy,subarray_y]; permute coordinates AND W together.
nx=accurate.sub_element_num_x; ny=accurate.sub_element_num_y;
mx=accurate.array_num_x; my=accurate.array_num_y;
indices=reshape(1:accurate.element_num,nx,ny,mx,my);
order=reshape(permute(indices,[1,3,2,4]),1,[]);
baseline=accurate;
baseline.array_num_x=1; baseline.array_num_y=1; baseline.array_num=1;
baseline.sub_element_num_x=nx*mx; baseline.sub_element_num_y=ny*my;
baseline.array_gap_x=nx*mx*accurate.element_gap;
baseline.array_gap_y=ny*my*accurate.element_gap;
elements=reshape(accurate.pos_element,3,[]);
baseline.pos_element=reshape(elements(:,order),3,nx*mx,ny*my);
baseline.pos_subarray_c=zeros(3,1);
baseline.W=accurate.W_vec(:,order); baseline.W_vec=baseline.W;
baseline.element_permutation=order;
% Keep the physical observation covariance; column permutation changes no noise.
baseline.W_vec_2=accurate.W_vec_2;
[~,baseline.steer_vec_noise_inv,baseline.steer_vec_noise_inv2]= ...
    prepare_grid_dictionary(baseline.pos_grid,baseline.W,baseline.pos_element, ...
    baseline.lambda,baseline.array_n);
end
