
%%% Script for Analysis of CyCIF tonsil dataset for  "Are low-variance Principal Components of single-cell omics data biologically meaningful, or just noise?"
%%% Aligns TONSILA to TONSILB and writes the warped coordinates 
%%% inputs: TONSILA.ome.tif.tif TONSILB.ome.tif  affine_sweep_results.mat  dataTONSILA.csv dataTONSILB.csv
%%% outputs: dataTONSILA_registered.csv
%%% Data source: https://elifesciences.org/articles/31657#content
%%% Part 1 (commented out, not run): affine parameter sweep. 27 combinations of optimizer settings were tried; the results and the best transform are stored in affine_sweep_results.mat.
%%% Part 2: The script loads the saved transform from affine_sweep_results.mat.
%% Part 1: Affine Parameter Sweep
% Tries 27 combinations of optimizer settings for affine registration
% of Tonsil A → Tonsil B. Saves all transforms + the best one.

%% Input data
% Directory containing Tonsil A/B images and per-cell measurements.
scriptDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(scriptDir, 'Data');
% %%%load the tiff images of tonsils 
% imgA = imread(fullfile(dataDir,'TONSILA.ome.tif.tif'));
% imgB = imread(fullfile(dataDir, 'TONSILB.ome.tif'));
% 
% %% Image preprocessing
% % Convert RGB images to grayscale for intensity-based registration.
%%% Tifs are RGB images
% imgA_grey = im2gray(imgA);
% imgB_grey = im2gray(imgB);
% 
% %%%Load in the per-cell measurements.
% dataA = readtable(fullfile(dataDir, 'dataTONSILA.csv'));
% dataB = readtable(fullfile(dataDir, 'dataTONSILB.csv'));
%%
% Evaluate all 27 combinations of optimizer settings.
% Registration quality is assessed by the Pearson correlation between grayscale Tonsil B and the affine-transformed Tonsil A image.
% The parameter combination yielding the highest correlation is retained.
% The saved .mat file (affine_sweep_results.mat) contains the results and the best-performing transform. 
% Uncomment to re-run with rng(42) to reproduce the sweep

% Best-performing parameters in the final analysis:
% MaximumIterations = 500, InitialRadius = 0.001, GrowthFactor = 1.15 (image correlation = 0.40372).




% 
% 
% % Set random seed for reproducibility.
% rng(42)
% 
% % Optimizer parameters evaluated in the sweep.
% % Maximum number of optimizer iterations.
% MaximumIter = [100,300,500];
% % Initial search radius in transformation-parameter space.
% InitialRad = [0.001, 0.00625, 0.01]; 
% % Growth factor controlling expansion of the optimizer search radius.
% GrowthFactor = [1.05,1.1,1.15];
% 
% % Generate all 27 combinations of optimizer parameters.
% [I, R, G] = ndgrid(1:3, 1:3, 1:3); 
% nTests = numel(I); % 27
% 
% % Build the table with parameter columns + empty results columns
% results = table( ...
%     MaximumIter(I(:))', ...
%     InitialRad(R(:))', ...
%     GrowthFactor(G(:))', ...
%     zeros(nTests, 1), ...       % correlation
%     zeros(nTests, 1), ...       % runtime in seconds
%     'VariableNames', {'Itr', 'Rad', 'Growth', 'Corr', 'Runtime'});

% Store the affine transformation estimated for each parameter combination.
%tforms = cell(nTests, 1);

% for k = 1:nTests 
%     % Configure optimizer for the current parameter combination.
%     itr = results.Itr(k);
%     rad = results.Rad(k);
%     growth = results.Growth(k);
%     tic; % Start timer for performance measurement
%     [optimizer, metric] = imregconfig('multimodal'); %% Use MATLAB's multimodal registration configuration to accommodate differences in intensity distributions between the two serial sections.
%     optimizer.MaximumIterations = itr;
%     optimizer.InitialRadius = rad;
%     optimizer.GrowthFactor = growth;
%     % Estimate affine transformation from Tonsil A to Tonsil B.
%     tform_test = imregtform(imgA_grey, imgB_grey, 'affine', optimizer, metric);
% 
%     % Apply affine to A, render in B's pixel space, measure 2D correlation
%     corr_value = corr2(imgB_grey, imwarp(imgA_grey, tform_test, 'OutputView', imref2d(size(imgB_grey)))); % Calculate correlation. 
%     % Record image correlation and runtime.
%     elapsed = toc; % Stop timer and record elapsed time.
%     results.Corr(k) = corr_value;
%     results.Runtime(k) = elapsed;
%     tforms{k} = tform_test;
%     disp(['[', num2str(k), '/', num2str(nTests), '] Itr:', num2str(itr), ...
%         ' Rad:', num2str(rad), ' Growth:', num2str(growth), ...
%         ' | Corr: ', num2str(corr_value, '%.5f'), ...
%         ' | Time: ', num2str(elapsed, '%.1f'), 's']);
% end
% 
% results_sorted = sortrows(results, 'Corr', 'descend');
% disp(results_sorted);
%%
%% Save affine-registration sweep
% % Identify the transform with the highest image correlation and save the complete sweep results together with the selected transform.

% [~, I] = max(results.Corr);
% trans_afine_best = tforms{I};
% 
% save(fullfile(dataDir, 'affine_sweep_results.mat'), ...
%        'results', 'results_sorted', 'tforms', 'trans_afine_best');


%% Part 2: Register Tonsil A cell coordinates into Tonsil B coordinate space
% Convert Tonsil A stage coordinates to image pixels, apply the selected affine transformation, and convert the transformed points to Tonsil B stage coordinates. 
% Cells outside the Tonsil B field of view are flagged.
% Output: dataTONSILA_registered.csv

%%  Load images, saved transform and per-cell tables (skipped if already in the workspace)
if ~exist('imgA_grey', 'var')
    imgA_grey = im2gray(imread(fullfile(dataDir, 'TONSILA.ome.tif.tif')));
end
if ~exist('imgB_grey', 'var')
    imgB_grey = im2gray(imread(fullfile(dataDir, 'TONSILB.ome.tif')));
end
if ~exist('trans_afine_best', 'var')
  load(fullfile(dataDir, 'affine_sweep_results.mat'), 'trans_afine_best');
end
if ~exist('dataA', 'var')
    dataA = readtable(fullfile(dataDir, 'dataTONSILA.csv'));
end
if ~exist('dataB', 'var')
    dataB = readtable(fullfile(dataDir, 'dataTONSILB.csv'));
end


%% Scale factors (microscope coordinates -> pixels)
% Derive independent X/Y scale factors for each section from the image dimensions and maximum observed stage coordinates.

scaleX_A = size(imgA_grey,2) / max(dataA.Xt);
scaleY_A = size(imgA_grey,1) / max(dataA.Yt);
scaleX_B = size(imgB_grey,2) / max(dataB.Xt);
scaleY_B = size(imgB_grey,1) / max(dataB.Yt);

% Tile dimensions in stage coordinates.
% The CyCIF microscope images tissue in a 7×8 grid of tiles at this spacing
tileW = 1664;
tileH = 1404;



%% Step 1: A coords -> A pixels (flip Y: image Y and stage Yt are inverted)
% Microscope stage Y increases upward; image row index increases downward.
% Subtracting from image height flips the axis to match.
Apx = [dataA.Xt * scaleX_A, size(imgA_grey,1) - dataA.Yt * scaleY_A];

%% Step 2: Apply affine (A pixels -> B pixels).
Bpx = transformPointsForward(trans_afine_best, Apx);


%% Step 3: B pixels -> B coords (un-flip Y).
dataA.Xt_aff     = Bpx(:,1) / scaleX_B;
dataA.Yt_aff     = (size(imgB_grey,1) - Bpx(:,2)) / scaleY_B;


%% Tile assignments after affine transform.
dataA.COL_aff     = floor(dataA.Xt_aff     / tileW) + 1;
dataA.ROW_aff     = floor(dataA.Yt_aff     / tileH) + 1;

%% Identify cells within the registered field of view.
% Cells transformed outside the coordinate extent of Tonsil B are flagged as invalid. 
dataA.valid = dataA.Xt_aff > 0 & dataA.Yt_aff > 0 & ...
    dataA.Xt_aff < max(dataB.Xt) & ...
    dataA.Yt_aff < max(dataB.Yt);

%% Write full table: all original A columns + the new ones appended.
writetable(dataA, fullfile(dataDir, 'dataTONSILA_registered.csv'));
%% Downstream analysis
% The registered per-cell coordinates generated here are used as input to
% downstream correction in "03_Flat-field_correction.Rmd"