%% --- 1. SETUP ---
clear; clc;
outputFolder = 'test_results'; % Use a new folder for this run
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

% Reproducible draw
rng(2);                      
N_cand   = 1000;
N_true   = 10;
N_total  = N_cand + N_true;

%% --- 2. GENERATE A UNIFORM HAMMERSLEY SET ---
hammersley_set = zeros(N_total, 2);
base = 2; % Prime base for the second dimension
for i = 1:N_total
    hammersley_set(i, 1) = i / N_total;
    n = i; inv = 0; p = 1 / base;
    while n > 0
        d = mod(n, base); inv = inv + d * p;
        p = p / base; n = floor(n / base);
    end
    hammersley_set(i, 2) = inv;
end

%% --- 3. SCALE ALL POINTS ---
p0_min = 8;      p0_max = 16;
stiff_min = 75000; stiff_max = 100000;
p0    = p0_min    + (p0_max    - p0_min)    .* hammersley_set(:,1);
stiff = stiff_min + (stiff_max - stiff_min) .* hammersley_set(:,2);

%% --- 4. RANDOMLY ASSIGN ROLES & PREPARE FOR PLOTTING ---
rand_indices = randperm(N_total);
true_indices = rand_indices(1:N_true);
cand_indices = rand_indices(N_true + 1:end);

% Create a "role map" for the loop to use
roles = cell(N_total, 1);
roles(:) = {'candidate'}; % Default all to candidate
roles(true_indices) = {'true'}; % Assign the true roles

%% --- 5. PLOTTING ---
figure;
hold on;
% Plot candidates (use cand_indices to select them from the full p0/stiff arrays)
scatter(p0(cand_indices), stiff(cand_indices), 25, 'b', 'filled', 'DisplayName', 'Candidates');
% Plot true values
scatter(p0(true_indices), stiff(true_indices), 40, 'r', 'filled', 'DisplayName', 'True Values');
hold off;
title('Hammersley Set with Randomly Selected True Values');
xlabel('p0'); ylabel('stiff');
legend('show', 'Location', 'best'); grid on;
axis([p0_min p0_max stiff_min stiff_max]);

%% --- 6. RUN SOLVER LOOP OVER ALL POINTS ---
outputLogFile = fullfile(outputFolder, 'loop_output.txt');

% Create and write header to the new log file
fid = fopen(outputLogFile, 'w');
fprintf(fid, 'Index,Role,Stiffness,P0,Filename\n');
fclose(fid);

fprintf('Starting solver loop for %d points...\n', N_total);
for i = 1:N_total
    stiffVal = round(stiff(i));
    p0Val    = p0(i);
    role     = roles{i}; % Get the role from our map

    fprintf('Processing point %d/%d: Role=%s, Stiff=%d, p0=%.3f\n', ...
            i, N_total, role, stiffVal, p0Val);

    % Clean p0 string for filename
    p0Str = regexprep(sprintf('%.3f', p0Val), {'0+$', '\.$'}, '');

    % ---- 1) Run solver ----
    cmd = sprintf('./sor06 %d 0.7 1.5 1.8 %.3f 16.0 %d%.3f', ...
                  stiffVal, p0Val, stiffVal, p0Val);
    % [status, cmdout] = system(cmd); % Use this for better error checking
    system(cmd);

    % ---- 2) Move produced .2d files ----
    system(sprintf('mv *.2d "%s/"', outputFolder));
    
    % ---- 3) Find and rename the correct output file ----
    cands = dir(fullfile(outputFolder, sprintf('output_%d*.2d', stiffVal)));
    if isempty(cands)
        warning('No output file found for point %d (Stiff=%d, p0=%.3f). Skipping.", i, stiffVal, p0Val');
        continue
    end
    [~, ix] = max([cands.datenum]);
    oldFile = fullfile(outputFolder, cands(ix).name);
    newName = sprintf('%s_%d_%s.2d', role, stiffVal, p0Str);
    newFile = fullfile(outputFolder, newName);
    movefile(oldFile, newFile, 'f');
    
    % ---- 4) Append to the new log file ----
    fid = fopen(outputLogFile, 'a');
    fprintf(fid, '%d,%s,%d,%.6f,%s\n', i, role, stiffVal, p0Val, newName);
    fclose(fid);
end
fprintf('Solver loop finished.\n');