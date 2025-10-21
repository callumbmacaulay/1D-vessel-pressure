%% --- 1. SETUP ---
clear; clc;
outputFolder = 'stiff_CT_2048_results';
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

% Parameter ranges
k3_min   = 2e5;  k3_max   = 9e5;   % stiffness (f3)
CT1_min  = 0.3;  CT1_max  = 1.8;   % Windkessel compliance to sweep
Rp1 = 0.7; Rd1 = 1.5;               % fixed Windkessel resistances
p0  = 8.0;                          % mmHg (inlet)
pout = 16.0;                        % mmHg (fixed outlet)

N_total = 2048;                     % total Sobol points
fprintf('Generating %d Sobol points...\n', N_total);

%% --- 2. GENERATE SOBOL SEQUENCE (built-in, no toolbox) ---
% 2 dimensions: [k3, CT1]
p = sobolset(2,'Skip',10,'Leap',20);
sobol_pts = net(p, N_total);

% Scale to desired parameter ranges
k3_vals  = k3_min  + (k3_max  - k3_min)  .* sobol_pts(:,1);
CT1_vals = CT1_min + (CT1_max - CT1_min) .* sobol_pts(:,2);

% Last one is "true", others "candidates"
roles = repmat({'candidate'}, N_total, 1);
roles{end} = 'true';

scatter(k3_vals, CT1_vals, 6, 'filled');
xlabel('k3'); ylabel('CT1'); title('Sobol sequence (scaled)'); grid on;
set(gca,'XScale','linear');   % <- remove the log scaling

%% --- 3. PREPARE LOG ---
outputLog = fullfile(outputFolder, 'sobol_run_log.csv');
fid = fopen(outputLog, 'w');
fprintf(fid, 'Index,Role,k3,CT1,OutputFile\n');
fclose(fid);

fprintf('Starting solver loop for %d runs...\n', N_total);

%% --- 4. MAIN LOOP ---
for i = 1:N_total
    k3   = k3_vals(i);
    CT1  = CT1_vals(i);
    role = roles{i};
    fileID = i; % used by solver for output_%d.2d

    fprintf('Run %03d/%03d: %s | k3 = %.2e | CT1 = %.3f\n', ...
        i, N_total, upper(role(1:3)), k3, CT1);

    % --- Run solver ---
    % Order must match the unchanged C program:
    % ./sor06 k3 Rp1 Rd1 CT1 p0 p_out fileID
    cmd = sprintf('./sor06 %.3e %.2f %.2f %.3f %.2f %.2f %d', ...
        k3, Rp1, Rd1, CT1, p0, pout, fileID);
    system(cmd);

    % --- Move and rename output file ---
    oldFile = sprintf('output_%d.2d', fileID);
    if exist(oldFile, 'file')
        newFile = sprintf('%s_k3_%g_CT1_%.3f.2d', role(1:3), k3, CT1);
        movefile(oldFile, fullfile(outputFolder, newFile));
    else
        warning('No output file found for run %d (k3=%.2e, CT1=%.3f).', i, k3, CT1);
        newFile = 'missing';
    end

    % --- Log results ---
    fid = fopen(outputLog, 'a');
    fprintf(fid, '%d,%s,%.3e,%.3f,%s\n', i, role, k3, CT1, newFile);
    fclose(fid);
end

fprintf('All %d Sobol simulations completed.\n', N_total);
