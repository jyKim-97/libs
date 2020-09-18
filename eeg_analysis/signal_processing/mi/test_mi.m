%% make test_signals
srate = 2000;
t = 0:1/srate:10;
y1 = cos(2*pi*6*t);
% y2 = exp(y1-1) .* cos(2*pi*50*t)/10;

id = y1 < -0.8;
y2 = id .* cos(2*pi*50*t)/5;
y = y1 + y2 + normrnd(0, 0.5, [1, length(t)])/10;

%% bandpass signal
f_low = 1:1:15; len_low = length(f_low)-1;
f_high = 20:5:150; len_high = length(f_high)-1;

passed_data = struct('low', [], 'high', []);

for i = 1:len_low
    passed_data.low(i, :) = bandpass_signal(y, f_low(i:i+1), srate, 10);
end
for i = 1:len_high
    passed_data.high(i, :) = bandpass_signal(y, f_high(i:i+1), srate, 10);
end

%% 
phs_signal = zeros(size(passed_data.low));
for i = 1:len_low
    phs_signal(i, :) = angle(hilbert(passed_data.low(i, :)));
end
amp_signal = zeros(size(passed_data.high));
for i = 1:len_high
    amp_signal(i, :) = abs(hilbert(passed_data.high(i, :)));
end

%%
bins = bin_phase(phs_signal, 40);
mi = getMI(bins, amp_signal, 40);
%%
mid_pts = @(x, i) (x(i+1) + x(i)) / 2;
f_low_mids = arrayfun(@(i) mid_pts(f_low, i), 1:len_low);
f_high_mids = arrayfun(@(i) mid_pts(f_high, i), 1:len_high);

figure("Units", "Normalized", "pos", [0.5, 0.1, 0.4, 0.8]);
subplot(2, 1, 1)
plot(t, y, "k")

subplot(2, 1, 2)
% imagesc(f_low_mids, f_high_mids, mi);
% ax.YDir = "normal";
contourf(f_low_mids, f_high_mids, mi', 100, "EdgeColor", "none");
colormap jet
colorbar