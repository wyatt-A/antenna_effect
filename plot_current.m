
clear; clc;

%% Configuration
Sim_Path = 'wire2';

% Must match the probe positions from simulation
L = 100;
probe_spacing = 5;

x_probe = unique([ ...
    -L/2 + (0.5:0.5:3), ...
    (-L/2+5):probe_spacing:(L/2-5), ...
     L/2 - (0.5:0.5:3), ...
     0]);

n_probes = numel(x_probe);

%% Load time-domain current data

for n = 1:n_probes
    name = sprintf('I_%03d', n);

    probe = ReadUI(name, Sim_Path);

    if n == 1
        t = probe.TD{1}.t(:);
        I = zeros(n_probes, numel(t));
    end

    I(n,:) = probe.TD{1}.val(:).';
end

fprintf('Loaded %d probes, %d timepoints\n', ...
    n_probes, numel(t));

fprintf('Time step: %.4f ps\n', ...
    mean(diff(t))*1e12);

fprintf('Total duration: %.2f ns\n', ...
    (t(end)-t(1))*1e9);

%% Plot center current vs time

[~, center_idx] = min(abs(x_probe));

figure;
plot(t*1e9, I(center_idx,:)*1e3, ...
    'LineWidth', 1);

xlabel('Time (ns)');
ylabel('Current (mA)');
title('Center current ring-down');
grid on;

%% Plot all current probes

figure;
imagesc(t*1e9, x_probe, I*1e3);
axis xy;
colorbar;

xlabel('Time (ns)');
ylabel('Position along wire (mm)');
title('Axial current distribution (mA)');

%% Plot maximum current at each position

I_max = max(abs(I), [], 2);

figure;
plot(x_probe, I_max*1e3, 'o-', ...
    'LineWidth', 1.5);

xlabel('Position along wire (mm)');
ylabel('Maximum |I(t)| (mA)');
title('Peak transient current distribution');
grid on;

%% Plot current distribution at selected time

t_target = 5e-9; % seconds
[~, idx] = min(abs(t-t_target));

figure;
plot(x_probe, I(:,idx)*1e3, 'o-', ...
    'LineWidth', 1.5);

xlabel('Position along wire (mm)');
ylabel('Current (mA)');
title(sprintf('Current at t = %.2f ns', ...
    t(idx)*1e9));
grid on;
