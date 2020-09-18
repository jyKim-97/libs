function ax_hd_eeg = get_hdeeg_axis(gap)
% if isexist(varargin, 'FigureSize')
%     figsize = varargin{isexist(varargin, 'FigureSize')};
%     fig = figure('units', 'normalized', 'position', [0.01, 0.01, figsize(1), figsize(2)]);
% else
%     fig = figure('position', get(0, 'ScreenSize'));
% end
if numel(gap) == 1
    gap = [gap, gap];
end
%% prepare figure
sz = [9, 6];
location = [3, 4, 9, 10, 8, 11, 15, 16, 14, 17, 21, 22, 19, 24, ...
    27, 28, 26, 29, 25, 30, 33, 34, 32, 35, 31, 36, 39, 40, ...
    38, 41, 37, 42, 45, 46, 44, 47, 51, 52];
% margin (bottom, left)
mh = 0.02;
mw = 0.02;
% axes height and width
axh = (1 - 2 * mh - (sz(1) - 1) * gap(1)) / sz(1);
axw = (1 - 2 * mw - (sz(2) - 1) * gap(2)) / sz(2);
% initialize
ph = 1 - axh - mw;
ax_hd_eeg = cell(1, 38);
% calculate
loc = 0;
ch = 0;
for i = 1:38 % num of channels
    quot = floor(location(i)/6);
    remd = mod(location(i), 6);
    if remd == 0
        remd = 6;
        quot = quot - 1;
    end
    pw = mw+(axw+gap(2))*(remd-1);
    ph = 1-axh-mw-(axh+gap(1))*(quot);
    
    ax_hd_eeg{i} = axes('Units', 'normalized',...
        'Position', [pw, ph, axw, axh]);
end

end

