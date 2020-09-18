function ind = assign_args(args, target)
% ind = read_args(args, target)
% #####input
% args, arguments (cell)
% target, target vars name (char)
% #####output
% ind, index of the target vars name (int or [])
% ############
% Last Updated: 2020-03-03
% JeongYeong
ind = cellfun(@(x) check_args(x, target), args);
if any(ind)
    ind = find(ind == 1, 1);
else
    ind = [];
end
end

function equiv = check_args(x, target)
% equiv = read_args(x, target)
if isstring(x) || ischar(x)
    equiv = strcmp(x, target);
else
    equiv = false;
end
end