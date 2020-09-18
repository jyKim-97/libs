function [xx, yy, dat2d] = get_griddata(data, topo, alive_ch)
% check outlier
sgm = std(data);
avg = mean(data);
ind1 = data < avg-2*sgm;
ind2 = data > avg+2*sgm;
data(ind1) = avg-2*sgm;
data(ind2) = avg+2*sgm;
% x -> column, y -> row
bval = prctile(data, 10);
xx = linspace(min(topo.boundary(:, 1)), max(topo.boundary(:, 1)), 100);
yy = linspace(min(topo.boundary(:, 2)), max(topo.boundary(:, 2)), 100);
[xx, yy] = meshgrid(xx, yy);
dat2d = griddata([topo.loc(alive_ch, 1); topo.boundary(:, 1)], [topo.loc(alive_ch, 2); topo.boundary(:, 2)],...
    [data; bval*ones(size(topo.boundary, 1), 1)], ...
    xx, yy, 'cubic');
end