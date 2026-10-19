function y_compress_noise = add_noise(y_vec,SNR,sysPara)
W_vec = sysPara.W_vec;

signal_pow = sum(abs(y_vec).^2,'all')/(numel(y_vec));
noise_pow = signal_pow/(10^(SNR/10));
noise = sqrt(noise_pow/2) * (randn(size(y_vec)) + 1j*randn(size(y_vec)));
y_vec_noise = y_vec + noise;
y_compress_noise = W_vec*y_vec_noise;

end