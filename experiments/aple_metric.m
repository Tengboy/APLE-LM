function s=aple_metric(mse,threshold)
% Unknown/nonfinite outcomes count as failures, never silently as successes.
v=mse(:); n=numel(v); assert(n>0 && threshold>0);
valid=isfinite(v) & v>=0; keep=valid & v<=threshold^2;
k=sum(~keep); z=1.959963984540054; p=k/n; den=1+z*z/n;
mid=(p+z*z/(2*n))/den; half=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den;
s=struct('attempted',n,'valid',sum(valid),'included',sum(keep),'exceeded',k, ...
    'rmse',NaN,'conditional_rmse',NaN,'outage',p,'wilson_low',max(0,mid-half),'wilson_high',min(1,mid+half));
if all(valid), s.rmse=sqrt(mean(v)); end
if any(keep), s.conditional_rmse=sqrt(mean(v(keep))); end
end
