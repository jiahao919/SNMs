function res = vect(x)
% Vectorize N dimensional input in one column
% This is the wrapper of matlab vectorization
%
%     Args:
%         x(Nd)     : Input data
%     Return:
%         res(1D)   : Vectorized output
% (c) Zheyuan_Yi 2018


%%  vectorization
    res = x(:);
end