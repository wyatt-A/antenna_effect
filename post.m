
clear; clc;

filename = 'wire_full/E_tip_pos.h5';
e_t = readmatrix("wire_full/et");

%% Find the time-domain group

info = h5info(filename, '/FieldData');

disp({info.Groups.Name});

% Find TD group
idx = find(endsWith(string({info.Groups.Name}), '/TD'), 1);

if isempty(idx)
    error('Cannot find time-domain group.');
end

td_path = info.Groups(idx).Name;
td_info = h5info(filename, td_path);

n_frames = numel(td_info.Datasets);

fprintf('Found %d time-domain frames\n', n_frames);

%% Load time-domain E-field

% First frame
ds0 = td_info.Datasets(1);
E0 = h5read(filename, [td_path '/' ds0.Name]);

% Allocate [Nx, Ny, Nz, 3, Nt]
E_td = zeros([size(E0), n_frames], 'single');

t = zeros(n_frames,1);
step = zeros(n_frames,1);

for k = 1:n_frames

    ds = td_info.Datasets(k);
    path = [td_path '/' ds.Name];

    E_td(:,:,:,:,k) = h5read(filename, path);

    % Dataset name is FDTD iteration
    step(k) = str2double(ds.Name);

    % Read timestamp with full precision
    t(k) = double(h5readatt(filename, path, 'time'));

end

%% Sort by time

[t, order] = sort(t);

E_td = E_td(:,:,:,:,order);
step = step(order);

%% Load spatial mesh

x = h5read(filename, '/Mesh/x');
y = h5read(filename, '/Mesh/y');
z = h5read(filename, '/Mesh/z');

%% Report

fprintf('Field dimensions: %s\n', mat2str(size(E_td)));
fprintf('Time range: %.6f to %.6f ns\n', ...
    t(1)*1e9, t(end)*1e9);

fprintf('Mean sampling interval: %.4f ps\n', ...
    mean(diff(t))*1e12);

fprintf('Sampling frequency: %.4f GHz\n', ...
    1/mean(diff(t))/1e9);

%%

%sweep = linspace(0.1e9,5e9,100);
sweep = 3.7e9;
e_total = [];
for f = sweep

    % Fourier weights
    w = exp(-1i*2*pi*f*t);

    % Fourier transform along time dimension
    E_mat = reshape(E_td, [], numel(t));

    dt = mean(diff(t));

    E_fd = reshape( ...
        dt * (double(E_mat) * w), ...
        size(E_td,1), ...
        size(E_td,2), ...
        size(E_td,3), ...
        3);

    t_exc = e_t(:,1);
    et = e_t(:,2);

    dt = mean(diff(t_exc));

    Et_fd = dt * sum( ...
        et .* exp(-1i*2*pi*f*t_exc));

    % E_fd contains complex field components
    E_CW = E_fd / Et_fd;

    % Magnitude of vector field
    Emag_CW = sqrt(sum(abs(E_CW).^2, 4));

    e_total(end+1) = sum(Emag_CW,"all");

    niftiwrite(Emag_CW,"ECW");
    
    
    niftiwrite(abs(E_CW(:,:,:,1)),"ECW_x");
    niftiwrite(abs(E_CW(:,:,:,2)),"ECW_y");
    niftiwrite(abs(E_CW(:,:,:,3)),"ECW_z");

end

hold on
plot(sweep,e_total)


