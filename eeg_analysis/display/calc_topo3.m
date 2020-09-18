function facecolor = calc_topo3(chdata, topo, varargin)
% data : EEG channel(1~38)에서 얻은 data 값
% bad_ch : data의 bad channel
% fname : 2D topography map을 그릴때 사용한 channel(loc) & boundary coordination이 들어있는 mat 파일 이름
% res_2D : 2D topography map resolution (default : 250)
% res_3D : 3D topography map resolution (default : 250)
res_2D = 250;
res_3D = 250;
bad_ch = [];
switch nargin
    case 3
        bad_ch = varargin{1};
    case 4
        res_2D = varargin{2};
    case 5
        res_3D = varargin{3};
end

b_value = prctile(chdata, 10);
alivech = setdiff(1:38, bad_ch);
% 2d interpolation
[~, ~, Z] = interp_topo(topo, chdata, alivech, res_2D, b_value);
%% make object
Xlin = linspace(min(topo.V_top(:,1)), max(topo.V_top(:,1)), res_3D);
Zlin = linspace(min(topo.V_top(:,3)), max(topo.V_top(:,3)), res_3D);
[xx, zz] = meshgrid(Xlin, Zlin);
% (0, 0)방향에서 축 (x, z) -> 2D 평면에서 (x, y)에 대응
Zlin_2D = linspace(min(topo.boundary(:,2)), max(topo.boundary(:,2)), res_2D);
% 보통 맨 앞 2ch은 ground로   쓰이므로 topography mapping 영역에서 제외
z_max = max(topo.loc(3:end, 2));
z_max_mesh = find(Zlin_2D > z_max, 1, 'first');
z_min = min(topo.loc(3:end, 2));
z_min_mesh = find(Zlin_2D < z_min, 1, 'last');
%% calculate
Z_new = nan(res_3D, res_3D); % res_3D
for i = 2:res_3D-1
    if Zlin(i) < z_min % posterior
        temp_mesh = z_min_mesh*(Zlin(i)-Zlin(1))/(z_min-Zlin(1));
        Z_mesh1 = fix(temp_mesh)+1;
        Z_mesh2 = fix(temp_mesh)+2;
        w_z1 = fix(temp_mesh)+1-temp_mesh;
        w_z2 = temp_mesh-fix(temp_mesh);
    elseif Zlin(i) > z_max % anterior
        temp_mesh = (res_3D-z_max_mesh)*(Zlin(i)-z_max)/(Zlin(res_3D)-z_max);
        Z_mesh1 = fix(temp_mesh)+z_max_mesh;
        Z_mesh2 = fix(temp_mesh)+z_max_mesh+1;
        w_z1 = fix(temp_mesh)+1-temp_mesh;
        w_z2 = temp_mesh-fix(temp_mesh);
    else % inside
        temp_mesh = (z_max_mesh-z_min_mesh)*(Zlin(i)-z_min)/(z_max-z_min);
        Z_mesh1 = fix(temp_mesh)+z_min_mesh; % 1 x 1
        Z_mesh2 = fix(temp_mesh)+z_min_mesh+1; % 1 x 1
        w_z1 = fix(temp_mesh)+1-temp_mesh;
        w_z2 = temp_mesh-fix(temp_mesh);
    end
    ind = (topo.V_top(:,3) > Zlin(i-1)) & (topo.V_top(:,3) < Zlin(i+1));
    x1 = min(topo.V_top(ind, 1));
    x2 = max(topo.V_top(ind, 1));
    ind1 = find(Xlin >= x1, 1, 'first');
    ind2 = find(Xlin <= x2, 1, 'last');
    Z_new(i, :) = w_z1 * interp_x(Z, Z_mesh1, ind1, ind2, res_3D)...
        + w_z2 * interp_x(Z, Z_mesh2, ind1, ind2, res_3D);
end
nst_grid = zeros(length(topo.ind_f), 2);
x_v = xx(~isnan(Z_new));
z_v = zz(~isnan(Z_new));
facecolor = nan(length(topo.ind_f), 1);
a = topo.g.vertices(topo.g.faces(topo.ind_f, 1), :);
b = topo.g.vertices(topo.g.faces(topo.ind_f, 2), :);
c = topo.g.vertices(topo.g.faces(topo.ind_f, 3), :);
% topo.ind_f가 무엇인지 아는 것이 먼저
% 왜 faces에 column이 3개나 존재하는거지?
top_f_g = (a+b+c) / 3;
for i = 1:length(topo.ind_f)
    [~, tmp_ind] = min((x_v - top_f_g(i, 1)).^2 + (z_v - top_f_g(i, 3)).^2);
    nst_grid(i, 1) = find(zz(:, 1) == z_v(tmp_ind));
    nst_grid(i, 2) = find(xx(1, :) == x_v(tmp_ind));
    facecolor(i, 1) = Z_new(nst_grid(i, 1), nst_grid(i, 2));
end
end

function [xx, yy, zz] = interp_topo(topo, data, ch, res, b_value)
% code for plot_topo3
num_ch=38;

eeg_ch = 1:num_ch; %% 실제 recording channel을 앞 2ch을 제외한 3~40번 channel로 정함
ind_b = length(eeg_ch)+1:length(topo.boundary)+length(eeg_ch);

% boundary와 electrode 위치 결합
X = topo.loc(eeg_ch,1); 
X(ind_b) = topo.boundary(:,1); 
Y = topo.loc(eeg_ch,2);
Y(ind_b) = topo.boundary(:,2);
% bounary를 감싸는 mesh 생성
Xlin = linspace(min(X), max(X), res);
Ylin = linspace(min(Y), max(Y), res);
[xx, yy] = meshgrid(Xlin, Ylin);
% grid data -> b_value(boundary value)을 boundary값으로 지정
zz = griddata([X(ch)', X(ind_b)'], [Y(ch)', Y(ind_b)'], [data(ch)', b_value*ones(1,length(ind_b))], xx, yy, 'cubic');
% boundary밖의 점들을 nan처리
in = inpolygon(xx(:), yy(:), topo.boundary(:, 1), topo.boundary(:, 2));
zz(~in) = nan;
end

function Z_new = interp_x(Z, z_mesh, ind1, ind2, res)
% code for plot_topo3
Z_new = ones(1, res);
ind3 = find(~isnan(Z(z_mesh, :)));
if ~isempty(ind3)
    ind = ind2 - ind1 + 1;
    for j = 1:ind
        v = length(ind3)/double((ind2-ind1+1)) * double(j);
        if v < 1
            Z_new(j+ind1-1) = Z(z_mesh,ind3(1));
        elseif  j == (ind2-ind1+1)
            Z_new(j+ind1-1) = Z(z_mesh,ind3(end));
        else
            Z_new(j+ind1-1) = Z(z_mesh,fix(v)+ind3(1)-1)*(fix(v)+1-v) + Z(z_mesh,fix(v)+ind3(1))*(v-fix(v));
        end
    end
end
end