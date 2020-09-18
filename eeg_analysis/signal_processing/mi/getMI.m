function varargout = getMI(low_bins, amp_signals, n_bins)
% mi, probs = getMI(low_bins, amp_signals, n_bins);
% use bin_phase.m first
% mi: len_low x len_high
% probs: nbins x len_low x len_high

%% get MI
len_low = size(low_bins, 1);
len_high = size(amp_signals, 1);

probs = zeros(len_high, len_low, n_bins);
for i = 1:len_low
    for j = 1:n_bins
        id = low_bins(i, :) == j;
        probs(:, i, j) = nanmean(amp_signals(:, id), 2);
    end
end

probs = probs ./ nansum(probs, 3);

%% get entropy
H = -nansum(probs .* log(probs), 3)';
H = squeeze(H);
%% get modulation index (MI)
Hmax = log(n_bins);
mi = (Hmax - H) ./ Hmax;
%% return
varargout{1} = mi;
if nargout == 2
    varargout{2} = permute(probs, [3, 2, 1]);
end
end
