t = linspace(0, 3, 1000);
f_sets = 1:20;
y = zeros(length(f_sets), length(t));
for i = 1:length(f_sets)
    y(i, :) = sin(2*pi*f_sets(i)*t);
end

%%
fig = figure('Units', 'Normalized', 'pos', [0.2, 0.3, 0.6, 0.3]);
ax1 = subplot(121);
contourf(t, f_sets, y, 20, 'EdgeColor', 'none'); colormap jet
datacursormode on
dcm_obj = datacursormode(fig);

ax2 = subplot(122);
p = plot(t, nan(1, length(t)), 'k');

handles.ax2 = ax2;
handles.p = p;
vars.y = y;
update_fcn = @(tmp, event_obj) live_plot_fn(tmp, event_obj, handles, vars);
set(dcm_obj, 'UpdateFcn', update_fcn); % if you use function script, add @ at the first. ex) @my_func

