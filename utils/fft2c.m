function res = fft2c(x)
% Calculate fast Fourier transform of first two dimensions
%
%     Args:
%         x           : Time domain/image space data
%     Return:
%         res         : Frequency domain/k space data
%
% (c) Zheyuan_Yi 2018

%% fft
% scaling factor
S = size(x);
fctr = S(1) * S(2);

x = reshape(x, S(1), S(2), prod(S(3:end)));
res = zeros(size(x));

for n = 1:size(x, 3)
    res(:, :, n) = 1 / sqrt(fctr) * fftshift(fft2(ifftshift(x(:, :, n))));
end

res = reshape(res, S);

end
