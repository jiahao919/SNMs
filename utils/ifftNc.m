function res = ifftNc(x, Dim, noShiftDim, pcSize)
% N-dimensional inverse fast Fourier transform
%
% Args:
%         x             : frequency domain / k space data
%         *Dim          : ifft along the dimension (default: all)
%         *noShiftDim   : dimension without ifftshift (default: none)
%         *pcSize       : crop/pad array with zeros in frequency domain (default: none)
% Return:
%         res           : time domain / image space data
% 
% Note:
% 1. if no need 'ifftshift' along specific dimension, use 'noShiftDim' 
%         [1]   --> first dimension 
%         [1,3] --> first and third dimension 
% 
% 2. data has been scaled by its size 
% 
% 3. length of 'PadSize' should match dimension of 'x' 
% 
% 4. crop/pad array would automatically re-scale data to keep same energy level
% 
% Require:
%    cropp, zpadd
% 
% (c) Zheyuan_Yi 2021

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

%% IFFT
res = cropp(x, cSize);
res = zpadd(res, pcSize);

for n = 1:ndims(x)
    if ismember(n, Dim)
        if ismember(n, noShiftDim)
            res = ifft(ifftshift(res, n), [], n);             
        else
            res = fftshift(ifft(ifftshift(res, n), [], n), n);
        end   
    end
end

%% Re-scale
fctr_ifft = prod(size(res, Dim))^(1/length(Dim));
fctr_pc = (prod(pcSize)/numel(x))^(1-1/length(Dim));

res = res * fctr_ifft * fctr_pc;

end