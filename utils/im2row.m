function res = im2row(im,winSize)
% Reshape image to row vectors from window
%
%     Args:
%         im(2D/3D)    : Images
%         winsize(1,2) : Reshape window size
%     Return:
%         res(2D/3D)   : row vectors matrix
%     Warning: 
%         1. 3D images is reshaped slice-wise. 
%         2. Window slide from top to bottom, then left to right.
%         3. Data point read from top to bottom, then left to right.
%         4. Example:
%            [10 11 12 13                       [10 14 11 15  row1
%             14 15 16 17       ---->            14 18 15 19  row2
%             18 19 20 21   winsize = [2 2]      18 22 19 23  row3
%             22 23 24 25]                       ...........]   
% (c) Zheyuan_Yi 2018

%% Input check
% if prod(winSize)>prod(size(im)-winSize+1)
%     error('Row number is less than column after transformation. Comment this error if no problem');
% end

%% Transform to row
[sx,sy,sz,scont] = size(im);

res = zeros((sx-winSize(1)+1),(sy-winSize(2)+1),prod(winSize),sz,scont);
count=0;
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;  
        res(:,:,count,:,:) = im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:,:);       
    end
end
res = reshape(res,[(sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz,scont]);

end

