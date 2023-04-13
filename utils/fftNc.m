function res = fftNc(x, Dim, noShiftDim, pcSize)
% N-dimensional fast Fourier transform
%
%     Args:
%         x             : time domain/image space data
%         *Dim          : fft along the dimension (default: all)
%         *noShiftDim   : dimension without fftshift (default: none)
%         *pcSize       : crop/pad array with zeros in frequency domain (default: none)
%     Return:
%         res           : frequency domain/k space data
% 
% Note:
% 1. if no need 'fftshift' along specific dimension, use 'noShiftDim' 
%         [1]   --> first dimension 
%         [1,3] --> first and third dimension 
% 
% 2. data has been scaled by its size 
% 
% 3. length of 'pcSize' should match dimension of 'x' 
% 
% 4. crop/Pad array would automatically re-scale data to keep same energy level
% 
% Require:
%    cropp, zpadd
% 
% (c) Zheyuan Yi 2021

%% Input check
if ~exist('noShiftDim', 'var')
    noShiftDim = 0;
end

if ~exist('Dim', 'var')
    Dim = 1:ndims(x);
end

if ~exist('pcSize', 'var')
    pcSize = size(x);
end

xSize = size(x); cSize = pcSize; 
cSize(cSize>xSize) = xSize(cSize>xSize); 

%% FFT
res = x;
for n = 1:ndims(x)
    if ismember(n, Dim)
        if ismember(n, noShiftDim)
            res = fft(ifftshift(res, n), [], n);
        else
            res = fftshift(fft(ifftshift(res, n), [], n), n);
        end
    end
end

res = cropp(res, cSize);
res = zpadd(res, pcSize);

%% Re-scale
fctr_fft = prod(size(x, Dim))^(1/length(Dim));
fctr_pc = (prod(pcSize)/numel(x))^(1-1/length(Dim));

res = res * (1/fctr_fft) * fctr_pc;

end
