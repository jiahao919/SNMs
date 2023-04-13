close all
clear
addpath('.\utils\')
addpath('.\ESPIRiT utils\')

%% load kData
fprintf('Reading and preprocessing data... \n')

load('.\data\brain.mat')
img_fully_sampled = img_6ch;
ksp_fully_sampled = fft2c(img_fully_sampled);
[nx, ny, nz, nc, nt] = size(img_fully_sampled);


ksp_fully_sampled_normalized = ksp_fully_sampled./prctile(vect(img_fully_sampled),99.995); % normalization

%% parameter selection
tol = 1e-4;
ncalib = 24;
R = 2;
lambda = 2e-3;
eigThresh_k = (8e-4).^(1/2); % 4e-4
ksize2 = [6 6]; 

%% Undersampling
mask=zeros([nx ny]);
mask(:,[ny/2:-R:1 ny/2:R:end])=1; % CC359 dataset
mask_nocalib = mask;
calib = cropp(ksp_fully_sampled,[nx,ncalib,nz,nc]);
calib_region = ones([nx,ncalib,nz,nc]);
calib_region = zpadd(calib_region,[nx,ny,nz,nc]);
mask_nocalib = mask_nocalib + calib_region;
mask_nocalib(mask_nocalib>=1) = 1;
mask = mask_nocalib;

tmp_slice_u = ifft2c(ksp_fully_sampled_normalized.*mask);
tmp_slice_f = ifft2c(ksp_fully_sampled_normalized);
                 
mask0 = mask; % original mask
mask = 1 - mask; % the mask indicating the unacquired points

%% 2d recon
sl_2d = 1; % slice to be reconstructed

mask0_2d = squeeze(mask0(:,:,sl_2d,:)); % 2D mask for the selected slice
mask_2d = squeeze(mask(:,:,sl_2d,:)); % 2D (unacquired point) mask for the selected slice
tmp_slice_2d = squeeze(tmp_slice_u(:,:,sl_2d,:)); % 2D alising image corresponding to the undersampled k-space
tmp_slice_f_2d = squeeze(tmp_slice_f(:,:,sl_2d,:)); % 2D reference image corresponding to the fully sampled k-space

ksp_calib = cropp(fft2c(tmp_slice_2d),[nx,ncalib,nc]);
[V_null_2d, S_2d, Mac_2d] = est_null_2d(ksp_calib, ksize2); % Calculate V matrix 
[kernel_matrix,S] = dat2Kernel(ksp_calib,ksize2); % Calculate singular values
idx = find(S >= S(1)*eigThresh_k, 1, 'last' ); % Find the cut-off according to the singular value threshold, i.e. sigma
Null_2d = maps_v2(V_null_2d(:,idx:end), size(tmp_slice_2d,1), size(tmp_slice_2d,2),size(tmp_slice_2d,3), ksize2); % v2 of generating nulling maps by constructing positive-definite matrix

% figure;immontage(abs(Null_2d),[0 1]) % show the positive-definite matrix of sptial nulling maps
    
tmp_recon = fft2c(tmp_slice_2d); % Acquired k-space part for data consistency
max_iter = 400;
S1.type='()';
S1.subs{:} = find(vect(mask_2d)); % Idx for the unaccquired (i.e., to be reconstructed) points
N = length(find(vect(mask_2d)));

tmp_zero_matrix = zeros(size(tmp_slice_2d));
K2I = @(x) ifft2c(reshape(x,size(tmp_slice_2d))); % k-space to image domain transform
I2K = @(x) fft2c(reshape(x,size(tmp_slice_2d))); % Image domain to k-space transform
rep = @(x) reshape(x,[size(tmp_slice_2d,1:3) 1 size(tmp_slice_2d,4)]); 
A = @(x) subsref(I2K(sum(Null_2d.*rep(K2I(subsasgn(tmp_zero_matrix,S1,x))),3)),S1); % Maps*Im 
b = - subsref(I2K(sum(Null_2d.*rep(K2I(tmp_recon)),3)),S1); % -Maps*Ia

tic
tmp_res = pcg(A, b, tol, max_iter, speye(N,N), speye(N,N), subsref(tmp_recon,S1)); % Maps*Ia = -Maps*Im
tmp_res = reshape(subsasgn(tmp_recon, S1, tmp_res),size(tmp_recon)); % Data consistency: assign the reconstructed points to the missing k-space part (same as assigning the acquired points to the reconstructed k-space)
tmp_res_im = ifft2c(tmp_res);


%% l1 regularization
imSize = size(tmp_res_im);
if length(imSize)>2
    imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2)))),imSize(3)];
else
    imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2))))];
end

XOP = Wavelet('Daubechies_TI',4,6);

if lambda > 0 
    tmp = zpadd(tmp_res_im,imSize_dyd);
    tmp = XOP*tmp;
    tmp = SoftThresh(tmp,lambda);
    tmp_res_im = XOP'*tmp;
    tmp_res_im = reshape(tmp_res_im,imSize_dyd);
    tmp_res_im = cropp(tmp_res_im,imSize);
end                                                       
tmp_res = fft2c(tmp_res_im);
toc

%% Calculate errors
recon_slice_2d = ifft2c(tmp_res);
err_map = tmp_slice_f_2d - recon_slice_2d;

%% Save
path = strcat('.\results\');
mkdir(path)
figure('Visible','off');immontage(sos(recon_slice_2d),[0 1])
exportgraphics(gcf,strcat(path,'recon.tiff'),'Resolution',300)
figure('Visible','off');immontage(sos(tmp_slice_f_2d),[0 1])
exportgraphics(gcf,strcat(path,'ref.tiff'),'Resolution',300)
figure('Visible','off');immontage(sos(err_map),[0 .2])
exportgraphics(gcf,strcat(path,'error x 5.tiff'),'Resolution',300)
figure('Visible','off');immontage(Null_2d,[0 1])
exportgraphics(gcf,strcat(path,'positive-definite matrix of sptial nulling maps.tiff'),'Resolution',300)

%% calculate NRMSE and SSIM
load('.\mask_MSE\brain_mask_mse.mat') % load mask within the object region for caculating metrics                            
quantitative_measurements = CalcPerf(sos(tmp_slice_f_2d).*mask_mse(:,:,sl_2d),sos(recon_slice_2d).*mask_mse(:,:,sl_2d));quantitative_measurements

%% functions
function [V_null, S, Mac] = est_null_2d(res, kernel_size)
% Estimate null-space basis of Hankel matrix

Mac = search_ACS_2d(res, kernel_size);

% compute null-space via SVD     
[V_null, S, ~] = eigs(Mac'*Mac, size(Mac,2)); 
% [~, ~, V_null] = svd(Mac, 'econ'); 

end

function Mac = search_ACS_2d(data, kernel_size)
% Mac: calibration Matrix constructed from ACS

Mac = im2row(data,kernel_size);
Mac = catSplit(Mac,2,3);

end

function [Null, M, rank] = est_rank_3d(V_null, Mac, rank, imsize, ksize, ref, tmp_rank, mask)

switch rank
    case 'Auto'
        S = svds(Mac,size(Mac,2));
        S = -gradient(S); 
        rank = find(S<0.01, 1);   
    case 'Auto_LMaFit'
        rank = lmafit_estRank(Mac);
    case 'Auto_GCV'
        rank = GCV_estRank(Mac,9);
    case 'NRMSE'
        if tmp_rank == 0 
            tmp_rank = 300:700;
        else
            if tmp_rank - 3 < 1
                tmp_rank = 1:tmp_rank+3;
            else
                tmp_rank = tmp_rank-3:tmp_rank+3;
            end
        end
        rank = NRMSE_estRank_eig_withmask_3d(V_null, ref, tmp_rank, imsize, ksize, mask);
%         rank = NRMSE_estRank_eig(V_null, ref, tmp_rank, imsize, ksize, 3);
end

[Null, M] = spatial_support_3d(V_null(:,rank:end), imsize(1), imsize(2), imsize(3), ksize);

end

function [Null, M] = spatial_support_3d(V_null, Nx, Ny, Nz, Nc, R)

numflt = size(V_null,2);       % number of filters
fltlen = size(V_null,1)/Nc;    % filter length


% back to k-space sets of kernels as in ESPIRiT
sps = reshape(V_null,[fltlen, Nc, numflt]);
sps = reshape(sps,[R(1), R(2), R(3), Nc, numflt]);


% calculating image-space spatial support via Eigen-decom 
[Null, M] = kernelEig_v1_3d(sps,[Nx Ny Nz]);
M = 1-M;

end

function Null_op = Null_Hermitian(Null)
    Null_op = conj(Null).*permute(Null,[1 2 3 5 4]);
    Null_op = -Null_op;
    Null_op(:,:,:,1:size(Null,4)+1:end) = 1 + Null_op(:,:,:,1:size(Null,4)+1:end);
end

function Null = maps_v2(V_null, Nx, Ny, Nc, R)

numflt = size(V_null,2);       % number of filters
fltlen = size(V_null,1)/Nc;    % filter length


% back to k-space sets of kernels 
sps = reshape(V_null,[fltlen, Nc, numflt]);
sps = reshape(sps,[R(1), R(2), Nc, numflt]); % ker_x*ker_y*num_ch*num_ker e.g. 6x6x8x245

% calculating image-space maps operator 
Null = fft2(conj(sps),2*R(1)-1,2*R(2)-1)/sqrt(prod(R)); % image-domain conj sym
Null = permute(Null,[1 2 5 3 4]);% Nx Ny 1 Nc Nk

Null_c = fft2(flip(flip(sps,1),2),2*R(1)-1,2*R(2)-1)/sqrt(prod(R));%image-domain sym
Null_c = permute(Null_c,[1 2 3 5 4]);% Nx Ny Nc 1 Nk

Null = sum(Null.*Null_c,5); %image domain
Null = ifft2(Null); %k-space
Null = zpadd(Null,[Nx, Ny, Nc, Nc])*sqrt(prod([Nx Ny]));
Null = ifft2c(Null);

end

function x = SoftThresh(data,thre)
    cut = thre/sqrt(5e-2);
	res = abs(data) - cut;
	res = 1/2*(res + abs(res));
	x   = exp(1j*angle(data)).*res;
end
