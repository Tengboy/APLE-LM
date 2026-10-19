function [g]=grad_numerical_vec(f,x0,dx)
    %计算数量场f在x0处的梯度g（数值计算，差分距离dx）仅适用于实数
    y = f(x0);
    output_len = length(y);
    n=length(x0);
    g=zeros(output_len,n);
    for i=1:n
        x1=x0;
        x1(i)=x1(i)+dx;
        g(:,i)=(f(x1)-y)/dx; % Reuse the unchanged base evaluation.
    end
end
