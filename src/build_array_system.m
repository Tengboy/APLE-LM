function sys=build_array_system(base,M,N,RF,W_seed,angular_grid_size)
% Square physical array, shallow parameter struct, unchanged numerical kernels.
% M subarrays per axis, N elements per subarray axis. W follows stored element order.
assert(all([M,N,RF,angular_grid_size]>=1) && all(mod([M,N,RF,angular_grid_size],1)==0));
sys=base; sys.array_num_x=M; sys.array_num_y=M; sys.array_num=M*M;
sys.sub_element_num_x=N; sys.sub_element_num_y=N; sys.RF_chain=RF;
sys.array_gap_x=N*sys.element_gap; sys.array_gap_y=N*sys.element_gap;
sys.element_num_x=M*N; sys.element_num_y=M*N; sys.element_num=(M*N)^2;
origins=zeros(3,M,M);
for ax=1:M
    for ay=1:M
        origins(:,ax,ay)=sys.array_gap_x*sys.array_dir_x*(ax-1)+sys.array_gap_y*sys.array_dir_y*(ay-1);
    end
end
centres=origins+sys.element_gap*sys.array_dir_x*(N-1)/2+sys.element_gap*sys.array_dir_y*(N-1)/2;
array_centre=sum(centres,[2,3])/sys.array_num;
elements=zeros(3,N,N,M,M);
for ax=1:M
    for ay=1:M
        for ix=1:N
            for iy=1:N
                elements(:,ix,iy,ax,ay)=origins(:,ax,ay)+sys.element_gap*sys.array_dir_x*(ix-1) ...
                    +sys.element_gap*sys.array_dir_y*(iy-1);
            end
        end
    end
end
sys.pos_element=elements-array_centre; sys.pos_subarray_c=centres-array_centre;
% Even arrays retain the same floor-index reference in the existing signal models.
theta_beg=(90-sys.theta_range/2)/180*pi; theta_end=(90+sys.theta_range/2)/180*pi;
xs=cos(theta_beg):-2/angular_grid_size:cos(theta_end); ys=xs;
radii=linspace(sys.rmin,sys.rmax,2); grid=zeros(3,numel(xs)*numel(ys)*2); count=0;
for x=xs
    for y=ys
        if x*x+y*y>(1-cos(sys.theta_range/2/180*pi)^2), continue; end
        for r=radii
            count=count+1; grid(:,count)=r*[sqrt(1-x*x-y*y);x;y];
        end
    end
end
sys.pos_grid=grid(:,1:count);
stream=RandStream('mt19937ar','Seed',W_seed,'NormalTransform','Ziggurat');
sys.W=exp(1j*2*pi*rand(stream,RF,N*N,M,M));
sys.W_vec=reshape(sys.W,RF,[]); sys.W_vec_2=sys.W_vec*sys.W_vec';
[~,sys.steer_vec_noise_inv,sys.steer_vec_noise_inv2,sys.noise_precision]= ...
    prepare_grid_dictionary(sys.pos_grid,sys.W,sys.pos_element,sys.lambda,sys.array_n);
end
