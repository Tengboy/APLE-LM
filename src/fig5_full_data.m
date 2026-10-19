function trial=fig5_full_data(seed,accurate,SNR)
% Generate all observations before invoking any estimator. Same data as fig5_trial.
previous=RandStream.getGlobalStream;
cleanup=onCleanup(@()RandStream.setGlobalStream(previous)); %#ok<NASGU>
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed,'NormalTransform','Ziggurat'));
[p,q]=gene_pos(accurate);
accurate.source_pos=p; accurate.source_pos_para=q;
[y,~,gain]=gene1_y(accurate);
J=cal_FIM(gain,accurate);
trial=struct('seed',seed,'SNR',SNR(:).','truth',p,'truth_parameters',q, ...
    'true_gain',gain,'true_signal',y,'observations',complex(zeros(accurate.RF_chain,numel(SNR))), ...
    'bound',nan(1,numel(SNR)),'noise_power',nan(1,numel(SNR)));
signal_power=sum(abs(y).^2,'all')/numel(y);
for s=1:numel(SNR)
    trial.observations(:,s)=add_noise(y,SNR(s),accurate);
    bound=cal_CRLB(y,SNR(s),J,accurate);
    trial.bound(s)=bound/accurate.source_num;
    trial.noise_power(s)=signal_power/10^(SNR(s)/10);
end
assert(all(isfinite(trial.bound) & trial.bound>0),'Invalid position bound.');
end
