function y = catSplit(x,Dim_cat,Dim_split,sizSplit)
% Concatenate to eliminate one dimension or split to form one dimension along specific dimension.
%
%     Args:
%         x(ND)          : Input array
%         Dim_cat(num)   : Cat along dimension
%         Dim_split(num) : Split along dimension
%         sizSplit(num)  : Size after split of Dim_Split, 1 means concatenate  
%     Return:
%         y(MD)          : Reshaped array
%     Warning: 
%         Only use to concatenate and split along one dimension at once.
%         Along multiple dimension, please call in multiple times.
%     Example:
%         1. catenate dimension 3 to dimension 1
%         size(A) = [256,256,32,4]
%         B = catSplit(A,1,3)
%         size(B) = [256*32,256,4]
%         2. catenate dimension 3 to dimension 1
%         size(A) = [256,256,32,4]
%         B = catSplit(A,1,3,4)
%         size(B) = [256*32/4,256,4,4]
% (c) Zheyuan_Yi 2018

%% Input Check
% if nargin<3
%     warning('please specify the dimension');
% end
if nargin < 4
    sizSplit = 1;  
end
% if size(x,Dim_split)/siz ~= fix(size(x,Dim_split)/siz)
%     warning('siz input is not a divisor of dim_split')
% end

%% Split and Concatnate
siz = size(x);
if sizSplit>1
    % Split
    index = 1:ndims(x)+1;
    index = [index(1:Dim_cat) Dim_split+1 index(Dim_cat:end)];
    index(Dim_cat+2*(Dim_cat~=1)) = [];
    index(Dim_split+1+(Dim_cat<Dim_split)) = [];
    
    siz_new = siz(Dim_split)/sizSplit;
    siz = [siz(1:Dim_split) siz_new siz(Dim_split:end)];
    siz(Dim_split+2) = [];
    siz(Dim_split) = sizSplit;       

    y = reshape(x,siz);
    y = permute(y,index);
else
    % Concatenate
    if ~(Dim_split>ndims(x))
    index = 1:ndims(x);
    index = [index(1:Dim_cat) Dim_split index(Dim_cat:end)];
    index(Dim_cat+2) = [];
    index(Dim_split+(Dim_cat<Dim_split)) = []; 

    siz(Dim_cat) = siz(Dim_cat)*siz(Dim_split);
    siz(Dim_split) = [];
    
    y = permute(x,index);
    y = reshape(y,siz);
    else 
    y = x;
    end
end

y = squeeze(y);

end