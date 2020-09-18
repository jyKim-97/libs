function txt = live_plot_fn(tmp, event_obj, handles, vars)
pos = get(event_obj, 'Position');
% get(gca, 'Childern')
y_loc = pos(2); 
handles.p.YData = vars.y(y_loc, :);
txt = {sprintf('time: %.3f s', pos(1)),...
    sprintf('freq: %d Hz', pos(2)),...
    sprintf('val: %.3f', pos(3))};
end