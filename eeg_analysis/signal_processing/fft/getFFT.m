function [xfft, f] = getFFT(x, srate, varargin)
% [xfft, f] = getFFT(x, srate, 'nx', 'OutputType')
% data (num_of_channels, length of time)

% default_vars
nx = size(x, 2);
isr = true;
% read vars
if ~isempty(varargin)
    % real / complex
    i = read_args(varargin, 'OutputType');
    if i
        output_type = varargin{i+1};
        if strcmp(output_type, 'real')
            isr = true;
        elseif strcmp(output_type, 'complex')
            isr = false;
        else
            error('');
        end
    end
    i = read_args(varargin, 'nx');
    if i
        nx = varargin{i+1};
    end
end
xfft = fft(x, nx, 2) / nx;
if isr
    xfft = abs(xfft);
end
    if mod(nx, 2) == 0 % even #
        f = srate*(0:nx/2)/nx;
        xfft = xfft(:, 1:nx/2+1);
        xfft(:, 2:end-1) = 2 * xfft(:, 2:end-1);
    else
        f = (0:(nx+1)/2)/(nx+1)*srate;
        xfft = xfft(:, 1:(nx+1)/2);
        xfft(:, 2:end) = 2 * xfft(:, 2:end);
    end
end
