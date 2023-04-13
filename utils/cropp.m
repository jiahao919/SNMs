function res = cropp(x,cx)
% crop N-dimensional array arround its center
%
%       Args:
%               x         [ndarray]  : array to crop
%               cx        [num]      : array size after crop
%       Return:
%               res       [ndarray]  : cropped array
% 
% Note:
% 1. Array will be squeezed after crop
% 
% 2. Consider N/2 + 1 as center for array with even size, 
% 
% (c) Zheyuan Yi 2018

%% Input check
if nargin < 2 
	error('please specify the target size')
end

if length(cx) ~= ndims(x)
	warning('crop size do not match the dimension of array')
end

%%  specify the size of cropping
sx = size(x); N =length(sx);M = length(cx);
if  M < N
    cx = [cx, ones(1,N-M)];
    if sx(1) == 1
        cx = fliplr(cx); M = 2;
    end
    if M == 1
        M = 2;
    end
end
    
if sum(sx<cx)>0
    error('size after crop must be smaller than before')
end
    
%%  cropping
idx = '';
for n = 1:M
    if mod(sx(n),2)          
        idx = sprintf('%s,%d:%d',idx,ceil(sx(n)/2)+1+floor(-cx(n)/2),ceil(sx(n)/2)+floor(cx(n)/2));             
    else
        idx = sprintf('%s,%d:%d',idx,floor(sx(n)/2)+1+ceil(-cx(n)/2),floor(sx(n)/2)+ceil(cx(n)/2));
    end
end
res = sprintf('res = squeeze(x(%s));',idx(2:end));
eval(res);

end