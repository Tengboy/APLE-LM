function r=aple_trial(seed,d,cid,bundle)
% Compute in memory, then return only per-trial scalars and a short status.
thread_guard=aple_threads; %#ok<NASGU>
sys=bundle.accurate; c=d.cases(cid); np=numel(d.points); nm=numel(d.methods);
r=struct('mse',nan(nm,np),'bound',nan(1,np),'nmse',nan(nm,np), ...
    'channel_bound',nan(1,np),'iterations',zeros(nm,np,'uint16'),'error','');
if strcmp(d.family,'coverage')
    data=fig10_data(seed,sys,15,c.theta,c.radius);
elseif strcmp(d.family,'iterations')
    % Consume the same six noise vectors as F9, then select the fourth (14 dB).
    data=fig5_full_data(seed,sys,5:3:20);
    mont=struct('y_compress_noise',data.observations(:,4),'inv_C_n',sys.noise_precision/data.noise_power(4));
    [p,tr]=fig11_trace(mont,sys,20,true);
    for j=1:20
        if all(isfinite(p(:,:,j)),'all'), [~,r.mse(j)]=match_user_positions(p(:,:,j),data.truth); end
    end
    r.iterations(1,:)=uint16(1:20); r.bound(:)=data.bound(4);
    bad=find(~cellfun(@isempty,tr.error),1); if ~isempty(bad), r.error=short_error(tr.error{bad}); end
    return;
else
    data=fig5_full_data(seed,sys,d.points);
end
r.bound=data.bound;
channel=strcmp(d.family,'snr');
if channel
    truth_sys=sys; truth_sys.source_pos=data.truth;
    [Btrue,physical]=fig6_channel_basis(data.truth,truth_sys,false);
    htrue=Btrue.*data.true_gain(:).'; signal=physical*data.true_gain(:);
    power=sum(abs(htrue).^2,1); [J,T]=cal_FIM(data.true_gain(:),truth_sys);
end
for j=1:np
    mont=struct('y_compress_noise',data.observations(:,j),'inv_C_n',sys.noise_precision/data.noise_power(j));
    if channel
        [~,~,Q,variance]=cal_CRLB_n(signal,d.points(j),J,T,truth_sys);
        r.channel_bound(j)=mean(variance./power);
    end
    for m=1:nm
        try
            method=d.methods{m}; info=struct('iterations',0);
            switch method
                case 'APLE-LM', [p,info]=fig5_estimate(mont,sys,false);
                case 'APLE-LM-ACM', [p,info]=fig5_estimate(mont,bundle.acm,true);
                case 'ES-GA', [p,~]=pos_ini_est(mont,bundle.esga);
                case 'P-SIGW'
                    [~,q]=pos_ini_est(mont,bundle.psigw); [p,~,info]=fig5_psigw(q,mont,bundle.psigw);
            end
            [permutation,r.mse(m,j)]=match_user_positions(p,data.truth);
            r.iterations(m,j)=uint16(info.iterations);
            if channel
                B=fig6_channel_basis(p,sys,m==2); weighted=B'*Q;
                gain=(weighted*B)\(weighted*data.observations(:,j)); estimate=B.*gain.';
                r.nmse(m,j)=mean(sum(abs(estimate(:,permutation)-htrue).^2,1)./power);
                assert(isfinite(r.nmse(m,j)),'APLE:InvalidChannel','Nonfinite channel error.');
            end
        catch err
            r.error=short_error([err.identifier ': ' err.message]);
        end
    end
end
end
function s=short_error(s)
s=s(1:min(numel(s),300));
end
