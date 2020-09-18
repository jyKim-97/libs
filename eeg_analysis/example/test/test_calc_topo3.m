%% load eeg data
load('test_data.mat');
topo = load('topo3_data.mat');
srate = 2000;
t = -1:1/srate:2;

%% test - get_hdeeg_axis.m
xffts = zeros(ceil(6001/(srate*0.05)), 128, 38);
for ch = 1:38
    [tmp, t_new, f_new] = getSTFFT(data(1, :), t, srate,...
        'nbin', 1000, 'maxf', 125, 'mbin', srate * 0.05, 'nf', 128);
    xffts(:, :, ch) = tmp;
end
%%
figure('Units', 'Normalized', 'pos', [0.05, 0.05, 0.9, 0.9]);
ax = get_hdeeg_axis(0.05);
for ch = 1:38
    axes(ax{ch})
    contourf(t_new, f_new, xffts(:, :, ch)', 20, 'EdgeColor', 'none')
    ylim([0, 60]);
    colormap jet
end

%% test - calc_topo3.m
% n = floor(size(data, 2)/10);
n = 20;
facecolors = zeros(length(topo.ind_f), n);
tic;
for i = 1:n
    facecolors(:, i) = calc_topo3(data(:, 20*i), topo, bad_ch, 250, 250);
end
fprintf('Done, execute time = %.3f\n', toc);

%%
plot_topo3(facecolors, topo, 'badchannels', bad_ch);

%%
[xx, yy, face] = get_griddata(data(:, 1), topo, 1:38);
figure; hold on
plot(topo.boundary(:, 1), topo.boundary(:, 2), 'k')
surf(xx, yy, face, 'EdgeColor', 'none', 'FaceAlpha', 0.9)
% axis equal
view(3);
colormap jet
set(gca, 'DataAspectRatio', [1, 1, 40])