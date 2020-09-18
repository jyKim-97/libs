function [xffts, t_new, f_new] = getSTFFT(x, t, srate, varargin)
% [xffts, t_new, f_new] = getSTFFT(x, t, srate, 'nx', 'wbin', 'mbin', 'maxf', 'nf', 'OutputType', 'Display')
% x (num_of_ch, length_of_times, num_of_trials)

% default vars
wbin = srate / 2;
mbin = 10;
maxf = nan;
nf = nan;
show = false;
nx = wbin;
output_type = 'real';
% read params
if ~isempty(varargin)
    % wbin
    i = read_args(varargin, 'wbin');
    if i
        wbin = varargin{i+1};
    end
    % mbin
    i = read_args(varargin, 'mbin');
    if i
        mbin = varargin{i+1};
    end
    % output type
    i = read_args(varargin, 'OutputType');
    if i
        output_type = varargin{i+1};
    end
    % nx
    i = read_args(varargin, 'nx');
    if i
        nx = varargin{i+1};
    end
    % nx
    i = read_args(varargin, 'nf');
    if i
        nf = varargin{i+1};
    end
    % maxf
    i = read_args(varargin, 'maxf');
    if i
        maxf = varargin{i+1};
    end
    % grphic option
    i = read_args(varargin, 'Display');
    if i
        show = varargin{i+1};
    end
end
% make window width
window = hann(wbin)'; 
% window = ones(1, wbin+1);
% run FFT
% ind = wbin/2:mbin:length(t)-wbin/2;
ind = 1:mbin:length(t);
t_new = t(ind);
sz = size(x);
if size(sz) < 3
    sz(3) = 1;
end
f_new = srate*(0:nx/2)/nx;
xffts = zeros(sz(1), length(f_new), sz(3), length(t_new));
for i = 1:length(ind)
    if ind(i)-wbin/2+1 < 0
        x2 = zeros(size(x, 1), wbin);
        len = ind(i)+wbin/2;
        x2(:, end-len+1:end) = x(:, 1:len);
    elseif ind(i)+wbin/2 > length(t)
        x2 = zeros(size(x, 1), wbin);
        len = length(t) - (ind(i)-wbin/2+1) + 1;
        x2(:, 1:len) = x(:, end-len+1:end);
    else
        x2 = x(:, ind(i)-wbin/2+1:ind(i)+wbin/2, :);
    end
    x2 = bsxfun(@times, x2, window);
    [xffts(:, :, :, i), ~] = getFFT(x2, srate, 'nx', nx, 'OutputType', output_type);
end
xffts = permute(xffts, [4, 2, 1, 3]); % (t, f, ch, trial)
% rescailing
if ~isnan(maxf)
    if isnan(nf)
        nf = 1000;
    end
    sz = size(xffts);
    sz(2) = nf;
    xffts_q = zeros(sz);
    f_new_q = linspace(f_new(1), maxf, nf);
    [X, Y] = meshgrid(t_new, f_new);
    [Xq, Yq] = meshgrid(t_new, f_new_q);
    for ch = 1:size(xffts, 3)
        for trial = 1:size(xffts, 4)
            tmp = interp2(X, Y, xffts(:, :, ch, trial)', Xq, Yq, "spline");
            xffts_q(:, :, ch, trial) = tmp';
        end
    end
    xffts = xffts_q;
    f_new = f_new_q;
end

if show
    figure;
    contourf(t_new, f_new, xffts(:, :, 1, 1)', 100, 'EdgeColor', 'none');
    colormap jet
    c = colorbar;
    c.Label.String = 'Amplitude';
    c.Label.FontSize = 12;
    xlabel('time', 'FontSize', 12)
    ylabel('frequency', 'FontSize', 12)
end
end