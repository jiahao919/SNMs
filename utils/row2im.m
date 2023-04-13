function res = row2im(mtx,imSize,winSize)
% Reshape row vector matrix to image from window
%
%     Args:
%         mtx(2D)      : Row vector matrix
%         imsize(1,2)  : Original image size
%         winsize(1,2) : Reshape window size
%     Return:
%         res(2D)   : Row vectors matrix
% (c) Zheyuan_Yi 2018

%% Reshape
sx = imSize(1);sy = imSize(2);
sz = imSize(3);scont = imSize(4);
res = zeros(sx,sy,sz,scont);

W1 = (1:winSize(1))'*(1:winSize(2));
W2 = repmat(W1(:,end),[1 sy-winSize(2)*2]);
W3 = [W1 W2 fliplr(W1)];
W4 = repmat(W3(end,:),[sx-winSize(1)*2 1]);
W5 = [W3;W4;flipud(W3)];
W = repmat(W5,[1 1 sz scont]);

mtx = reshape(mtx,[sx-winSize(1)+1,sy-winSize(2)+1,prod(winSize),sz,scont]);

count=0;
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;
        res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:,:) = res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:,:) + squeeze(mtx(:,:,count,:,:));
    end
end

res = res./W;

end
