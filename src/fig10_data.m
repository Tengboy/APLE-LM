function record=fig10_data(seed,sys,SNR,theta,radius)
% Fixed single-user location: only gain phase and noise change across seeds.
assert(sys.source_num==1 && theta>=0 && theta<=60 && radius>=sys.rmin && radius<=sys.rmax);
previous=RandStream.getGlobalStream;
cleanup=onCleanup(@()RandStream.setGlobalStream(previous)); %#ok<NASGU>
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed,'NormalTransform','Ziggurat'));
truth=radius*[cosd(theta);sind(theta);0]; q=[sind(theta);0;radius];
sys.source_pos=truth; sys.source_pos_para=q;
[signal,~,gain]=gene1_y(sys); J=cal_FIM(gain,sys);
observation=add_noise(signal,SNR,sys); bound=cal_CRLB(signal,SNR,J,sys);
power=sum(abs(signal).^2,'all')/numel(signal)/10^(SNR/10);
assert(isfinite(bound) && bound>0 && isfinite(power) && power>0,'F10:InvalidInput','Invalid conditional bound/noise power.');
record=struct('seed',seed,'SNR',SNR,'theta_deg',theta,'radius_m',radius,'azimuth_deg',0, ...
    'truth',truth,'truth_parameters',q,'true_gain',gain,'true_signal',signal, ...
    'observations',observation,'bound',bound,'noise_power',power);
end
