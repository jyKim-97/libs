function bins = bin_phase(phs_signals, n_bins)
edges = linspace(-pi, pi, n_bins+1);
bins = discretize(phs_signals, edges);
end