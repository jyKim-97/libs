function ps = get_amp_distrib(xfp, yfa, nbin)
% ps = get_amp_distrib(xfp, yfa, nbin)
% ##### input
% xfp, low freqs phase data
% yfa, high freqs amp data
% nbin, number of bins (distribution)
% ##### output
% ps, probability distribution, (nbin, nl, nh)

len_phase = size(xfp, 1);
len_amp = size(yfa, 1);
ps = zeros(nbin, len_phase, len_amp);
ind_low = discretize(xfp, linspace(-pi, pi, nbin+1));
for nl = 1:len_phase
    for nh = 1:len_amp
        for nb = 1:nbin
            ps(nb, nl, nh) = nanmean(yfa(nh, ind_low(nl, :) == nb));
        end
        ps(:, nl, nh) = ps(:, nl, nh) / sum(ps(:, nl, nh));
    end
end
end
