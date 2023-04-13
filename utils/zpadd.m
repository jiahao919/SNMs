function res = zpadd(x, zx)
% Padding zeros arround its center
%
%     Args:
%         x(ND)   : Center data
%         zx(1,M) : Size of data after zeros padding
%     Return:
%         res(Nd) : Padded data
%     Warning:
%         1. Dimension size not specified or equals to origin will not be padded
%         2. zpad([1 2 3],[6]) = [0 0 1 2 3 0]
%
% (c) Zheyuan Yi 2018

%% Input check
if nargin ~= 2
    error('must have a target size')
end

%%  specify the size of padding
sx = size(x); N = length(sx); M = length(zx);

if M < N

    if sx(1) == 1
        zx = [zx, sx(1)];
        zx = fliplr(zx); M = 2;
    else
        zx = [zx, sx(M + 1:end)];
    end

    if M == 1
        M = 2;
    end

end

% if sum(sx>zx)>0
%     error('size after padding must not smaller than before')
% end

%%  zeros padding
res = zeros(zx);
idx = '';

for n = 1:M

    if mod(siz(x, n), 2)
        idx = sprintf('%s,%d:%d', idx, floor(zx(n) / 2) + 1 + ceil(-siz(x, n) / 2), floor(zx(n) / 2) + ceil(siz(x, n) / 2));
    else
        idx = sprintf('%s,%d:%d', idx, ceil(zx(n) / 2) + 1 + floor(-siz(x, n) / 2), ceil(zx(n) / 2) + floor(siz(x, n) / 2));
    end

end

cmd = sprintf('res(%s) = x;', idx(2:end));
eval(cmd);

end
