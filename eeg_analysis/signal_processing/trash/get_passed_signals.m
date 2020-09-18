function [xfp, yfa] = get_passed_signals(x, srate, low_freqs, high_freqs, fo_l, fo_h)
% hsignals = get_passed_signals(data, low_freqs, high_freqs, srate, fo)
% Last Updated: 20.03.19
% ###### input
% data, 1 dim vector
% low_freqs, low freqs edges ([1, 2, 3] -> [1, 2], [2, 3])
% high_freqs, high freqs edges 
% srate, signal rate
% fo, filter order
% ##### output
% xfp
% yfa

len_low = length(low_freqs)-1;
len_high = length(high_freqs)-1;
%
if size(x, 1) > 1
    y = x(2, :);
    x = x(1, :);
end
%
xfp = zeros(len_low, length(x));
yfa = zeros(len_high, length(x));
for nl = 1:len_low
    pass_signal = bandpass_signal(x, [low_freqs(nl), low_freqs(nl+1)], srate, fo_l);
    xfp(nl, :) = angle(hilbert(pass_signal));
end
for nh = 1:len_high
    pass_signal = bandpass_signal(y, [high_freqs(nh), high_freqs(nh+1)], srate, fo_h);
    yfa(nh, :) = abs(hilbert(pass_signal));
end
end