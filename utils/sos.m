function res = sos(x, dim, pnorm)
    % Computes the square root of sum-of-squares along dim
    %
    %     Args:
    %         x       : multi-channel data
    %         dim     : Along specified dimension (default: last dimension)
    %         pnorm   : Power root and sum of power (default = 2)
    %     Return:
    %         res     : Square root of sum-of-squares data
    %
    % (c) Zheyuan_Yi 2018

    %% Input check
    if nargin < 3
        pnorm = 2;

        if nargin < 2
            dim = ndims(x);
        end

    end

    %% SOS
    res = squeeze((sum(abs(x.^pnorm), dim)).^(1 / pnorm));

end
