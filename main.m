
clear; clc;
physical_constants;

%% User parameters
unit = 1e-3;              % mm
L = 100;                  % wire length, mm
wire_radius = 0.25;       % wire radius, mm

f0 = 1.5e9;
fc = 1.0e9;
f_eval = linspace(0.5e9, 2.5e9, 401);

% Mesh parameters (mm)
dx_center = 2.0;
dx_tip = 0.25;
tip_region = 5.0;
dyz_center = 0.5;
dx_outer = 4.0;

% Probe parameters
probe_spacing = 5.0;
probe_halfwidth = 2.0;

c0=3e8;
% Simulation geometry
lambda_min = c0/(f0+fc)/unit;
domain_half = max(L/2 + 70, lambda_min);
tfsf_half = domain_half - 30;

Sim_Path = 'air';
Sim_CSX = 'air.xml';
% Sim_Path = 'wire_full';
% Sim_CSX = 'wire_full.xml';

if ~exist(Sim_Path, 'dir')
    mkdir(Sim_Path);
end

%% Initialize FDTD
FDTD = InitFDTD( ...
    'NrTS', 10000, ...
    'EndCriteria', 1e-6);

FDTD = SetGaussExcite(FDTD, f0, fc);

FDTD = SetBoundaryCond(FDTD, ...
    {'PML_8','PML_8','PML_8', ...
     'PML_8','PML_8','PML_8'});

CSX = InitCSX();

%% Construct nonuniform mesh

% Uniform central region
x_center = (-L/2 + tip_region):dx_center: ...
           ( L/2 - tip_region);

% Refined regions near each tip
x_tip_neg = (-L/2-tip_region):dx_tip: ...
            (-L/2+tip_region);

x_tip_pos = (L/2-tip_region):dx_tip: ...
            (L/2+tip_region);

% Important positions
x_fixed = [-domain_half, -tfsf_half, ...
           -L/2, 0, L/2, ...
           tfsf_half, domain_half];

% Smooth outer regions
mesh.x = unique([x_center, x_tip_neg, ...
                 x_tip_pos, x_fixed]);

mesh.x = SmoothMeshLines(mesh.x, ...
    dx_outer, 1.3);

% Transverse mesh near wire
yz_inner = -4:dyz_center:4;

yz_fixed = [-domain_half, -tfsf_half, ...
            -10, -4, 0, 4, 10, ...
            tfsf_half, domain_half];

mesh.y = unique([yz_inner yz_fixed]);
mesh.y = SmoothMeshLines(mesh.y, ...
    dx_outer, 1.3);

mesh.z = mesh.y;

% Ensure probe coordinates lie on grid
x_probe = unique([ ...
    -L/2 + (0.5:0.5:3), ...
    (-L/2+5):probe_spacing:(L/2-5), ...
     L/2 - (0.5:0.5:3), ...
     0]);

mesh.x = unique([mesh.x x_probe]);

% Current contour boundaries
mesh.y = unique([mesh.y, ...
    -probe_halfwidth, probe_halfwidth]);
mesh.z = unique([mesh.z, ...
    -probe_halfwidth, probe_halfwidth]);

% Define FDTD mesh
CSX = DefineRectGrid(CSX, unit, mesh);

fprintf('Mesh: %d x %d x %d cells\n', ...
    numel(mesh.x)-1, ...
    numel(mesh.y)-1, ...
    numel(mesh.z)-1);

%% Thin-wire conductor
% CSX = AddMetal(CSX, 'wire');
% 
% points = [
%     -L/2, L/2;
%      0,   0;
%      0,   0
% ];
% 
% CSX = AddWire(CSX,'wire',10,points,wire_radius);

%% Plane-wave excitation
% Propagation along z, E parallel to wire
k_dir = [0 0 1];
E_dir = [1 0 0];

CSX = AddPlaneWaveExcite(CSX, ...
    'plane_wave', k_dir, E_dir, f0);

CSX = AddBox(CSX, 'plane_wave', 0, ...
    -tfsf_half*[1 1 1], ...
     tfsf_half*[1 1 1]);

%% Current probes along wire
for n = 1:numel(x_probe)
    name = sprintf('I_%03d', n);
    x = x_probe(n);

    CSX = AddProbe(CSX, name, 1, ...
        'NormDir', 0);

    CSX = AddBox(CSX, name, 0, ...
        [x, -probe_halfwidth, -probe_halfwidth], ...
        [x,  probe_halfwidth,  probe_halfwidth]);
end

%% Frequency-domain E-field dumps near tips
tip_dump_size = [8 6 6];

% CSX = AddDump(CSX, 'E_tip_pos', ...
%     'DumpType', 10, ...
%     'FileType', 1, ...
%     'DumpMode', 2, ...
%     'Frequency', f0);

CSX = AddDump(CSX, 'E_tip_pos', ...
    'DumpType', 0, ...
    'FileType', 1, ...
    'DumpMode', 2);

CSX = AddBox(CSX, 'E_tip_pos', 0, ...
    [ L/2, 0, 0] - tip_dump_size, ...
    [ L/2, 0, 0] + tip_dump_size);

% CSX = AddDump(CSX, 'E_tip_neg', ...
%     'DumpType', 10, ...
%     'FileType', 1, ...
%     'DumpMode', 2, ...
%     'Frequency', f0);

CSX = AddDump(CSX, 'E_tip_neg', ...
    'DumpType', 0, ...
    'FileType', 1, ...
    'DumpMode', 2);

CSX = AddBox(CSX, 'E_tip_neg', 0, ...
    [-L/2, 0, 0] - tip_dump_size, ...
    [-L/2, 0, 0] + tip_dump_size);

%% Write and run
xml_file = fullfile(Sim_Path, Sim_CSX);

WriteOpenEMS(xml_file, FDTD, CSX);
CSXGeomPlot(xml_file);

RunOpenEMS(Sim_Path, Sim_CSX);

% %% Read current spectra
% I_fd = complex(zeros(numel(x_probe), ...
%                      numel(f_eval)));
% 
% for n = 1:numel(x_probe)
%     name = sprintf('I_%03d', n);
% 
%     p = ReadUI(name, Sim_Path, f_eval);
% 
%     I_fd(n,:) = reshape(p.FD{1}.val, 1, []);
% end
% 
% %% Current spectrum at wire center
% [~, center_idx] = min(abs(x_probe));
% 
% figure;
% plot(f_eval/1e9, ...
%      abs(I_fd(center_idx,:)), ...
%      'LineWidth', 2);
% xlabel('Frequency (GHz)');
% ylabel('Raw spectral current');
% title('Center current spectrum');
% grid on;
% 
% %% Spatial current distribution
% [~, f_idx] = min(abs(f_eval - 1.35e9));
% 
% I_spatial = abs(I_fd(:,f_idx));
% 
% figure;
% plot(x_probe, I_spatial, 'o-', ...
%     'LineWidth', 1.5);
% 
% xlabel('Position along wire (mm)');
% ylabel('Raw spectral current');
% title(sprintf('Current distribution at %.3f GHz', ...
%     f_eval(f_idx)/1e9));
% grid on;
