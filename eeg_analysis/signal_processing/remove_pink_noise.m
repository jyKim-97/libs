function varargout = remove_pink_noise(x, f, varargin)
% [xr, pink_noise(optional)] = remove_pink_noise(x, f, id(optional))
% x: fft result (len_t, len_f)
% f: frequency (recommend to use f > 1 Hz)
% id: baseline index, default: all index
% updated by JungYoung
switch nargin
    case 2
        id = true(1, size(x, 1));
    case 3
        id = varargin{1};
    otherwise
        error("too many arguments\n")
end
base = mean(x(id, :), 1);
%% need f > 1
id_f = f > 1;
f1 = f(id_f);
base = base(id_f);
p = polyfit(log(f1), log(base), 1);
pink_noise = exp(p(1)*log(f)+p(2));
%% subtract fitted line
xr = bsxfun(@minus, x, pink_noise);
xr(:, ~id_f) = nan;
%% return variables
varargout{1} = xr;
if nargout == 2
    varargout{2} = pink_noise;
end
end