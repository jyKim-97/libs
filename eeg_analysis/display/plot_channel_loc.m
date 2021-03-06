topo = load('topo3_data.mat');

label = {'FP1','FP2','AF3','AF4','AF7','AF8',...
    'F1','F2','F5','F6','FC1','FC2','FC5','FC6',...
    'C1','C2','C3','C4','C5','C6',...
    'CP1','CP2','CP3','CP4','CP5','CP6',...
    'P1','P2','P3','P4','P5','P6',...
    'PO3','PO4','PO7','PO8','O1','O2'};
figure;
plot(topo.boundary(:, 1), topo.boundary(:, 2), 'k.', 'MarkerSIze', 12);
axis equal;
ylim([-7.5, 6]);
for i = 1:38
    text(topo.loc(i, 1), topo.loc(i, 2), [num2str(i), '.', label{i}],...
        'HorizontalAlign', 'center', 'FontWeight', 'bold');
end
set(gca, 'XTick', []);
set(gca, 'YTick', []);