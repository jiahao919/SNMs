function y_size = siz(x, Dim)
    % Calculate size of specific dim
    %
    %     Args:
    %         x        : variable
    %         Dim(1,M) : Specified dimension (default: all dimensions)
    %     Inputs(easy mode):
    %         Variable name
    %     Return:
    %         y_size   : Size of variable
    %
    % (c) Zheyuan Yi 2018

    %% Input easy mode
    if nargin == 0
        str = input('Varname: ', 's');
        x = evalin('caller', str);
    end

    %% Get size
    if nargin <= 1
        y_size = size(x);
    end

    if nargin > 1

        if strcmp(Dim, 'end')
            Dim = length(y_size);
        end

        y_size = ones(1, length(Dim));

        for i = 1:length(Dim)
            y_size(i) = size(x, Dim(i));
        end

    end

end
