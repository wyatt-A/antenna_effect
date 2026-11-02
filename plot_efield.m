
clear; clc; close all;

%% Configuration
Sim_Path = 'thin_wire_air';
dump_name = 'E_tip_pos.h5';  % Update to match your file
filename = fullfile(Sim_Path, dump_name);

L = 100;                 % Wire length (mm)
tip_x = L/2;             % Right tip at +50 mm
window = 5;              % Display region +/- 5 mm

%% Load 3D field
[field, mesh] = ReadHDF5Dump(filename);

% First frequency-domain dataset
E = field.FD.values{1};

% Mesh coordinates (assuming meters -> mm)
x = mesh.lines{1}(:)*1e3;
y = mesh.lines{2}(:)*1e3;
z = mesh.lines{3}(:)*1e3;

% Complex vector magnitude
Emag = sqrt( ...
    abs(E(:,:,:,1)).^2 + ...
    abs(E(:,:,:,2)).^2 + ...
    abs(E(:,:,:,3)).^2);

fprintf('Field dimensions: %s\n', mat2str(size(E)));

%% Find slices through the tip

[~, ix] = min(abs(x-tip_x));
[~, iy] = min(abs(y));
[~, iz] = min(abs(z));

% Orthogonal slices
E_xy = squeeze(Emag(:,:,iz));
E_xz = squeeze(Emag(:,iy,:));
E_yz = squeeze(Emag(ix,:,:));

% Shared color scale in dB
Eref = max(Emag(:));
to_db = @(a) 20*log10(max(a/Eref,1e-8));

%% Plot
figure('Color','w','Position',[100 100 1300 420]);

subplot(1,3,1);
imagesc(x,y,to_db(E_xy).');
axis xy equal;
xlim(tip_x + [-window window]);
ylim([-window window]);
clim([-60 0]);
colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
title('x-y plane (z = 0)');

subplot(1,3,2);
imagesc(x,z,to_db(E_xz).');
axis xy equal;
xlim(tip_x + [-window window]);
ylim([-window window]);
clim([-60 0]);
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title('x-z plane (y = 0)');

subplot(1,3,3);
imagesc(y,z,to_db(E_yz).');
axis xy equal;
xlim([-window window]);
ylim([-window window]);
clim([-60 0]);
colorbar;
xlabel('y (mm)');
ylabel('z (mm)');
title(sprintf('y-z plane (x = %.2f mm)',x(ix)));

colormap(turbo);
sgtitle('Electric-field magnitude near wire tip');

%%

%% E-field quiver visualization
% E: complex array [Nx, Ny, Nz, 3]
% x, y, z: spatial coordinates in mm

% Settings
phase_deg = 0;          % RF phase in degrees
arrow_spacing = 2;      % Plot every 2nd vector
tip_x = 50;             % Wire tip position (mm)
window = 5;             % +/- 5 mm around tip

%% Extract x-z slice at y = 0

[~, iy] = min(abs(y));

Ex = squeeze(E(:,iy,:,1));
Ez = squeeze(E(:,iy,:,3));

% Full 3D phasor magnitude on this plane
Emag = squeeze(sqrt(sum(abs(E(:,iy,:,:)).^2, 4)));

%% Reconstruct instantaneous E-field

phi = deg2rad(phase_deg);

Ex_t = real(Ex * exp(1i*phi));
Ez_t = real(Ez * exp(1i*phi));

%% Create coordinate grids

[X,Z] = ndgrid(x,z);

%% Plot magnitude background

figure('Color','w');

imagesc(x,z,20*log10(max(Emag,eps)).');
axis xy equal;
hold on;

colormap(turbo);
colorbar;

xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('E-field at RF phase %.0f deg', phase_deg));

%% Downsample arrows

ix = 1:arrow_spacing:numel(x);
iz = 1:arrow_spacing:numel(z);

Xq = X(ix,iz);
Zq = Z(ix,iz);

Uq = Ex_t(ix,iz);
Vq = Ez_t(ix,iz);

%% Normalize arrows to emphasize direction
% Set normalize_arrows = false to retain relative magnitudes.

normalize_arrows = true;

if normalize_arrows
    A = hypot(Uq,Vq);
    Uq = Uq ./ max(A,eps);
    Vq = Vq ./ max(A,eps);
end

quiver(Xq,Zq,Uq,Vq,0.5, ...
    'k','LineWidth',1);

%% Zoom around wire tip

xlim(tip_x + [-window window]);
ylim([-window window]);

plot(tip_x,0,'wo', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',7);

%%
% Sample at z approximately 0.5 mm
[~, iz] = min(abs(z - 0.5));

Ex_line = squeeze(E(:,iy,iz,1));
Ez_line = squeeze(E(:,iy,iz,3));

figure;

plot(x, abs(Ex_line), 'LineWidth', 1.5);
hold on;
plot(x, abs(Ez_line), 'LineWidth', 1.5);

xline(50, '--k', 'Wire tip');

xlabel('x (mm)');
ylabel('Field magnitude');
legend('|E_x|','|E_z|');
grid on;