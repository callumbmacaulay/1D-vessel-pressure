%% --- 1. SETUP ---
clear; clc;
outputFolder = 'sobol_1024_results';
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

% Parameter ranges
k3_min = 2e5; k3_max = 9e5;
pout_min = 15; pout_max = 30;
Rp1 = 0.7; Rd1 = 1.5; CT1 = 1.8; p0 = 8.0;

N_total = 1024;             % total Sobol points
fprintf('Generating %d Sobol points...\n', N_total);

%% --- 2. GENERATE SOBOL SEQUENCE (built-in, no toolbox) ---
p = sobolset(2,'Skip',10,'Leap',20); % 2 dimensions: k3, p_out
sobol_pts = net(p, N_total);

% Scale to desired parameter ranges
k3_vals   = k3_min   + (k3_max   - k3_min)   .* sobol_pts(:,1);
pout_vals = pout_min + (pout_max - pout_min) .* sobol_pts(:,2);

% Last one is "true", others "candidates"
roles = repmat({'candidate'}, N_total, 1);
roles{end} = 'true';

%% --- 3. PREPARE LOG ---
outputLog = fullfile(outputFolder, 'sobol_run_log.csv');
fid = fopen(outputLog, 'w');
fprintf(fid, 'Index,Role,k3,p_out,OutputFile\n');
fclose(fid);

fprintf('Starting solver loop for %d runs...\n', N_total);

%% --- 4. MAIN LOOP ---
for i = 1:N_total
    k3   = k3_vals(i);
    pout = pout_vals(i);
    role = roles{i};
    fileID = i; % used by solver for output_%d.2d

    fprintf('Run %03d/%03d: %s | k3 = %.2e | p_out = %.2f mmHg\n', ...
        i, N_total, upper(role(1:3)), k3, pout);

    % --- Run solver ---
    cmd = sprintf('./sor06 %.3e %.2f %.2f %.2f %.2f %.2f %d', ...
        k3, Rp1, Rd1, CT1, p0, pout, fileID);
    system(cmd);

    % --- Move and rename output file ---
    oldFile = sprintf('output_%d.2d', fileID);
    if exist(oldFile, 'file')
        newFile = sprintf('%s_k3_%g_pout_%.2f.2d', role(1:3), k3, pout);
        movefile(oldFile, fullfile(outputFolder, newFile));
    else
        warning('No output file found for run %d (k3=%.2e, p_out=%.2f).', ...
                i, k3, pout);
        newFile = 'missing';
    end

    % --- Log results ---
    fid = fopen(outputLog, 'a');
    fprintf(fid, '%d,%s,%.3e,%.2f,%s\n', i, role, k3, pout, newFile);
    fclose(fid);
end

fprintf('All %d Sobol simulations completed.\n', N_total);
