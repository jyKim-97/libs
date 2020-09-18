function passed_signals = bandpass_signal_all(x, f_edges, srate, fo)
% passed_signals = bandpass_signal_all(x, f_range, srate, fo)
% x: signal, 1Xn
% f_edges: frequency edges
len_f = length(f_edges)-1;
passed_signals = zeros([len_f, length(x)]);
for i = 1:len_f
    passed_signals(i) = bandpass_signal(x, f_edges(i:i+1), srate, fo);
end
end