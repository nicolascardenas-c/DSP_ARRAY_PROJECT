%-------------------------------------------------------------------------
%% Direction of Arrival Estimation: Bartlett, MVDR and MUSIC
%-------------------------------------------------------------------------
clear; clc; close all;


% Locate the datasets relative to this script.
dataFolder = fileparts(mfilename('fullpath'));
if isempty(dataFolder)
    dataFolder = pwd;
end

files = {'spcom_10sep.mat', 'spcom_50sep.mat'};
methods = {'Bartlett', 'MVDR', 'MUSIC'};
theta = -90:0.01:90;                % Angles tested by the estimators [deg]

% Extract the datasets only if they are not already available.
if ~isfile(fullfile(dataFolder, files{1})) || ...
   ~isfile(fullfile(dataFolder, files{2}))

    zipFile = fullfile(dataFolder, 'spcom.zip');
    if ~isfile(zipFile)
        error('Place spcom.zip or both .mat files next to this script: %s', dataFolder);
    end
    unzip(zipFile, dataFolder);
end

% Parameters used in the two sensitivity experiments.
snapshotCounts = [5 10 20];       % Number of snapshots to retain
snrLevels = [30 20 10 0];         % Added-noise SNR [dB]
nTrials = 100;                    % Independent noise realizations for Monte Carlo
tolerance = 2;                    % All DoAs must be within 2 deg for success
minProminence = 3;                % Peak must stand out by at least 3 dB

% Reference DoAs inferred from the datasets (not provided as labels).
referenceAngles = {[20 30], [-20 30]};

%-------------------------------------------------------------------------
%% MODEL AND VARIABLES
%-------------------------------------------------------------------------
% Received data:  X = A_true*S + Noise
%
% M     : Number of antennas (5)
% N     : Number of snapshots (20 in the original datasets)
% K     : Number of sources (2, known from the dataset)
% Delta : Antenna spacing in wavelengths (0.5)
%
% X     : Received signals [M x N]; each column is one snapshot
% A     : Steering vectors [M x number of candidate angles]
% R     : Sample spatial covariance [M x M], estimated as X*X'/N
%
% All methods use the same covariance R. MVDR has no regularization.
% The source directions are estimated from peaks of the angular spectra.

for i = 1:numel(files)

    %---------------------------------------------------------------------
    %% 1. Load data and construct the array model
    %---------------------------------------------------------------------
    data = load(fullfile(dataFolder, files{i}));
    X = data.X;
    [M, N] = size(X);
    K = double(data.num_sources);
    Delta = double(data.Delta);

    m = (0:M-1).';
    A = exp(1j*2*pi*Delta * (m*sind(theta)));
    % For a uniform linear array, an arrival angle determines the phase
    % difference across antennas. Every column of A is one tested angle.

    R = X*X'/N;
    % This is the sample spatial covariance. It summarizes how the
    % antennas vary together, using all N available snapshots.

    fprintf('\n==================================================\n');
    fprintf('%s\n', files{i});
    fprintf('==================================================\n');

    %---------------------------------------------------------------------
    %% PART I - ORIGINAL DATA
    %---------------------------------------------------------------------
    % Apply the three methods to the dataset without modifying X.
    % First inspect the spatial spectra, then summarize the eigenvalues.
     
    % Compute the three spatial spectra and covariance eigenvalues.
    [P, eigenvalues] = doaSpectra(R, A, K);

    % Plot Bartlett, MVDR and MUSIC on the same angular axis.
    figure('Color','w', 'Name',[files{i} ' - Original spectra']);
    plotSpectra(theta, P);
    title('Original data');
    legend(methods, 'Location','best');
    
    % Print the estimated source angles for each method.
    fprintf('\nOriginal data:\n');
    printDoAs(theta, P, K, methods, minProminence);

    % Plot the eigenvalues to inspect the signal and noise subspaces.
    figure('Color','w', 'Name',[files{i} ' - Original eigenvalues']);
    plotEigenvalues(eigenvalues, K);
    title('Original covariance eigenvalues');

    %---------------------------------------------------------------------
    %% PART II - NUMBER OF SNAPSHOTS
    %---------------------------------------------------------------------
    % Keep the first 5, 10 or 20 columns of the ORIGINAL observations.
    % We change only the number of observations used to estimate R.
    % With N = M = 5, R may be poorly conditioned; MVDR is sensitive
    % because it solves a linear system involving this matrix.

    % One subplot per snapshot count; save eigenvalues for comparison.
    figure('Color','w', 'Name',[files{i} ' - Snapshots spectra']);
    eigSnapshots = zeros(M, numel(snapshotCounts));

    for c = 1:numel(snapshotCounts)
        % Keep only the first N_short observations and recompute R.
        N_short = snapshotCounts(c);
        X_short = X(:,1:N_short);
        R_short = X_short*X_short'/N_short;

        % Reapply the same three methods to the reduced dataset.
        [P, eigenvalues] = doaSpectra(R_short, A, K);
        eigSnapshots(:,c) = eigenvalues;  % Save results for the summary figure

        % Plot the resulting spectra for this observation length.
        subplot(1,numel(snapshotCounts),c);
        plotSpectra(theta, P);
        title(sprintf('N = %d snapshots', N_short));
        if c == 1
            legend(methods, 'Location','best');
        end

        % Print estimated angles to check how they change with N.
        fprintf('\nN = %d snapshots:\n', N_short);
        printDoAs(theta, P, K, methods, minProminence);
    end

    
    % Compare covariance eigenvalues for all three snapshot counts.
    figure('Color','w', 'Name',[files{i} ' - Snapshots eigenvalues']);
    plotEigenvalues(eigSnapshots, K);
    title('Covariance eigenvalues - number of snapshots');
    legend(compose('N = %d', snapshotCounts), 'Location','best');

    %---------------------------------------------------------------------
    %% PART III - ADDING NOISE
    %---------------------------------------------------------------------
    % Start from the original X and add spatially white complex Gaussian
    % noise. One noise realization is reused at every SNR for the PLOTS:
    % only its amplitude changes, making the comparison easier to see.
    % SNR is defined relative to the power of the original X, which
    % already contains a very small amount of noise.

    % Generate one reproducible Gaussian noise matrix for this dataset.
    rng(i);                          % Reproducible visual comparison
    W = (randn(M,N) + 1j*randn(M,N))/sqrt(2);
    signalPower = mean(abs(X(:)).^2);

    % Compare the spatial spectra at different noise levels.
    figure('Color','w', 'Name',[files{i} ' - Noise spectra']);
    eigNoise = zeros(M, numel(snrLevels));

    for c = 1:numel(snrLevels)
        snrDB = snrLevels(c);

        % Scale the noise by SNR, then recompute the spatial covariance.
        % W has unit power on average (before this scaling).
        noiseStd = sqrt(signalPower / 10^(snrDB/10));
        X_noisy = X + noiseStd*W;
        R_noisy = X_noisy*X_noisy'/N;

        % Estimate the three spectra using the noisy observations.
        [P, eigenvalues] = doaSpectra(R_noisy, A, K);
        eigNoise(:,c) = eigenvalues;  % Save results for the summary figure

        % Plot the spectra corresponding to the current SNR.
        subplot(2,2,c);
        plotSpectra(theta, P);
        title(sprintf('SNR = %d dB', snrDB));
        if c == 1
            legend(methods, 'Location','best');
        end

        % Print the DoA estimates to compare their sensitivity to noise.
        fprintf('\nSNR = %d dB:\n', snrDB);
        printDoAs(theta, P, K, methods, minProminence);
    end

    % Eigenvalue summary for all noise levels in ONE plot.
    figure('Color','w', 'Name',[files{i} ' - Noise eigenvalues']);
    plotEigenvalues(eigNoise, K);
    title('Covariance eigenvalues - added noise');
    legend(compose('%d dB', snrLevels), 'Location','best');

    %---------------------------------------------------------------------
    %% PART IV - MONTE CARLO (NOISE ROBUSTNESS)
    %---------------------------------------------------------------------
    % The spectra above describe ONE noise realization. To assess whether
    % the estimates are reliable, repeat the same noise experiment with
    % 100 independent noise realizations and summarize the ANGLE errors.
    %
    % Success rate: percentage of trials in which ALL K estimated angles
    %               differ by at most 'tolerance' degrees from reference.
    % RMSE: angular error over trials with K significant detected peaks.
    %       Incorrect estimates are included; trials with fewer than K
    %       peaks count as failures but cannot contribute to the RMSE.

    % Use a different reproducible seed and a coarser angular grid.
    rng(100+i);
    thetaMC = -90:0.1:90;           % Smaller grid to speed up Monte Carlo
    A_MC = exp(1j*2*pi*Delta * (m*sind(thetaMC)));

    % Accumulate successful trials, valid detections and squared errors.
    nSuccess = zeros(numel(snrLevels), numel(methods));
    nDetected = zeros(numel(snrLevels), numel(methods));
    squaredError = zeros(numel(snrLevels), numel(methods));

    % Repeat the experiment with new noise in every Monte Carlo trial.
    for trial = 1:nTrials
        W = (randn(M,N) + 1j*randn(M,N))/sqrt(2);
        % New noise in every trial, but the SAME W is used across SNRs.

        for c = 1:numel(snrLevels)
            noiseStd = sqrt(signalPower / 10^(snrLevels(c)/10));
            X_noisy = X + noiseStd*W;
            R_noisy = X_noisy*X_noisy'/N;

            % Estimate the spectra independently for this realization.
            [P, ~] = doaSpectra(R_noisy, A_MC, K);

            for j = 1:numel(methods)
                % Extract up to K DoA candidates from this method's spectrum.
                estimates = peakAngles(thetaMC, P(j,:), K, minProminence);

                % Compute angular errors only when K peaks are detected.
                if numel(estimates) == K
                    errors = estimates - referenceAngles{i};
                    nDetected(c,j) = nDetected(c,j) + 1;
                    squaredError(c,j) = squaredError(c,j) + sum(errors.^2);

                    % Count a success only if every DoA is within tolerance.
                    if all(abs(errors) <= tolerance)
                        nSuccess(c,j) = nSuccess(c,j) + 1;
                    end
                end
            end
        end
    end

    % Print success rate and conditional RMSE for every SNR and method.
    fprintf('\nMONTE CARLO - %s (%d trials)\n', files{i}, nTrials);
    fprintf('%7s  %-9s  %10s  %13s\n', ...
        'SNR', 'Method', 'Success %', 'RMSE [deg]');

    for c = 1:numel(snrLevels)
        for j = 1:numel(methods)
            successRate = 100*nSuccess(c,j)/nTrials;
            rmse = NaN;

            % RMSE is undefined if no trial produced K detected peaks.
            if nDetected(c,j) > 0
                rmse = sqrt(squaredError(c,j)/(K*nDetected(c,j)));
            end

            fprintf('%5d dB  %-9s  %9.1f%%  %13.2f\n', ...
                snrLevels(c), methods{j}, successRate, rmse);
        end
    end
end

%-------------------------------------------------------------------------
%% FUNCTIONS
%-------------------------------------------------------------------------
% These functions avoid repeating the same estimation, plotting and
% peak-detection steps. Only the input data change between experiments.

function [P, eigenvalues] = doaSpectra(R, A, K)
% Compute the three DoA spectra from one covariance matrix.
% P has one row per method and one column per candidate angle.

    M = size(R,1);

    % Bartlett: scan each angle with a conventional, non-adaptive
    % beamformer. The result is the spatial output power.
    P_bartlett = real(sum(conj(A).*(R*A),1))/M;

    % MVDR: choose weights to minimize output power while preserving
    % unit gain towards the tested angle (no diagonal regularization).
    % R\A solves R*Y=A, rather than explicitly computing inv(R)*A.
    Y = R\A;
    P_mvdr = 1 ./ max(real(sum(conj(A).*Y,1)), realmin);

    % MUSIC: separate signal and noise subspaces using the eigenvectors.
    % For a correct DoA, A is almost orthogonal to the noise subspace,
    % giving a small denominator and a peak in the pseudospectrum.
    [V,D] = eig(R);
    [eigenvalues, order] = sort(real(diag(D)), 'descend');
    E_noise = V(:,order(K+1:end));   % M-K smallest eigenvalue directions
    P_music = 1 ./ max(sum(abs(E_noise'*A).^2,1), realmin);

    % Stack the three spectra in the order given by 'methods'.
    P = [P_bartlett; P_mvdr; P_music];
end

function plotSpectra(theta, P)
% Normalize each spectrum separately to 0 dB for visual comparison.
% MUSIC is a pseudospectrum: its peak height is not received source power.

    % Plot one normalized curve per method on the same angular axis.
    hold on;
    for j = 1:size(P,1)
        p = P(j,:);
        p_dB = 10*log10(max(p/max(p),1e-12));
        plot(theta, p_dB, 'LineWidth',1.2);
    end
    xlabel('Angle [degrees]');
    ylabel('Normalized spectrum [dB]');
    xlim([-90 90]);
    ylim([-50 2]);
    grid on;
    hold off;
end

function plotEigenvalues(eigenvalues, K)
% Plot the eigenvalues of one or several covariance matrices.
% Each COLUMN contains the M eigenvalues of one experimental condition,
% already ordered from largest to smallest by doaSpectra().
% This allows comparing different snapshot counts or SNRs in one plot.

    eigenvalues = max(eigenvalues,realmin);  % Avoid log(0)
    M = size(eigenvalues,1);
    semilogy(1:M, eigenvalues, '-o', 'LineWidth',1.2);
    xlabel('Eigenvalue index');
    ylabel('Eigenvalue (log scale)');
    xticks(1:M);
    grid on;

    % K dominant eigenvalues describe the signal subspace in our model;
    % the remaining M-K are associated with the noise subspace.
    hold on;
    xline(K+0.5, 'k:', 'HandleVisibility','off');
    hold off;
end

function printDoAs(theta, P, K, methods, minProminence)
% Print at most K significant DoA candidates for each method.
% If only one strong peak exists, do not invent a second source angle.

    % Apply the same peak-selection rule to every estimator.
    for j = 1:numel(methods)
        estimates = peakAngles(theta, P(j,:), K, minProminence);
        fprintf('  %-9s: %s deg (%d/%d peaks)\n', ...
            methods{j}, mat2str(estimates,4), numel(estimates), K);
    end
end

function estimates = peakAngles(theta, p, K, minProminence)
% Find local maxima without Signal Processing Toolbox.
% Use peak PROMINENCE rather than a fixed fraction of the highest peak:
% a valid MUSIC peak can be much lower than another MUSIC peak.
%
% Prominence is the height of a peak above the higher surrounding valley.
% We discard shallow sidelobes and keep up to K strongest remaining peaks.
% This detection rule is practical, not part of the estimator itself.

    % Normalize to dB and first find every local maximum.
    p_dB = 10*log10(max(p/max(p),1e-12));
    peaks = find(p_dB(2:end-1) > p_dB(1:end-2) & ...
                 p_dB(2:end-1) >= p_dB(3:end)) + 1;

    keep = false(size(peaks));
    for k = 1:numel(peaks)
        q = peaks(k);

        % On either side, consider the valley up to the next higher peak
        % (or the edge of the angular grid if there is no higher peak).
        left = find(p_dB(1:q-1) > p_dB(q), 1, 'last');
        if isempty(left)
            left = 1;
        end

        right = find(p_dB(q+1:end) > p_dB(q), 1, 'first');
        if isempty(right)
            right = numel(p_dB);
        else
            right = q + right;
        end

        leftValley = min(p_dB(left:q));
        rightValley = min(p_dB(q:right));
        prominence = p_dB(q) - max(leftValley,rightValley);
        keep(k) = prominence >= minProminence;
    end

    % Keep sufficiently prominent peaks, ranked by spectral height.
    peaks = peaks(keep);
    [~,order] = sort(p(peaks),'descend');
    peaks = peaks(order(1:min(K,numel(order))));
    estimates = sort(theta(peaks));
end
