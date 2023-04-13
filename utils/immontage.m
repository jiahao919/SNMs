function immontage(img, range, shape, rot, fl)
% Display series of images as a montage.
%
%     Args:
%         x          : Input images in series
%         shape(1,2) : The shape of the montage
%         range(1,2) : Window level (similarly to imshow)
%         rot        : rotation of images (1: clockwise 90 degree) (default = 0)
%         fl         : flip images (1: left2right 2:up2down) (default = 0)
%     Inputs(easy mode):
%         Variable name
%     Return:
%         Imshow the montage image
%     Note:
%         Dimension after second will all treat as series
%         Range will auto-determined to [min max] if not specified
%         Shape will be auto-padded to square if not specified
%
% (c) Zheyuan Yi 2018

%% Input easy mode
if nargin == 0
    str = input('VarName: ', 's');
    img = evalin('caller', str);
end

img = img(:, :, :);


% display range
if ~exist('range', 'var') || isempty(range)

    if min(abs(img(:))) ~= max(abs(img(:)))
        range = [min(abs(img(:))), max(abs(img(:)))];
    else
        range = [0 1];
    end

end

if ~exist('rot', 'var') || isempty(rot)
    rot = 0;
end

if ~exist('fl', 'var') || isempty(fl)
    fl = 0;
end

img = rot90(img, rot);

switch fl
    case 0
    case 1
        img = fliplr(img);
    case 2
        img = flipud(img);
end

[sx, sy, nc] = size(img);

% specify shape of display
if ~exist('shape', 'var') || isempty(shape)

    if ceil(sqrt(nc))^2 ~= nc
        nc = ceil(sqrt(nc))^2;
        img(end, end, nc) = 0;
    end

    img = reshape(img, sx, sy * nc);
    img = permute(img, [2, 3, 1]);
    img = reshape(img, sy * sqrt(nc), sqrt(nc), sx);
    img = permute(img, [3, 2, 1]);
    img = reshape(img, sx * sqrt(nc), sy * sqrt(nc));

else
    img = reshape(img, sx, sy * nc);
    img = permute(img, [2, 3, 1]);
    img = reshape(img, sy * shape(2), shape(1), sx);
    img = permute(img, [3, 2, 1]);
    img = reshape(img, sx * shape(1), sy * shape(2));
end

imshow(abs(img), 'DisplayRange', range);
% imshow(img, 'DisplayRange', range);
% imcontrast;

end
