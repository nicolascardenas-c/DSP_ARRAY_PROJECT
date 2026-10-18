%% Direction of Arrival Estimation: Bartlett, MVDR and MUSIC
clear; clc; close all;

% Load the two datasets (spcom.zip must be in the current folder)
unzip('spcom.zip');
files = {'spcom_10sep.mat', 'spcom_50sep.mat'};
theta = -90:0.01:90;                    % Candidate angles [degrees]
methods = {'Bartlett', 'MVDR', 'MUSIC'};

figure('Color', 'w');

%% MODEL AND VARIABLES
%
% Model: X = A_true*S + Noise
% Goal: estimate the Directions of Arrival (DoA) of K sources.
%
% M       : Number of antennas (5)
% N       : Number of snapshots (20)
% K       : Number of sources (2)
% Delta   : Antenna spacing in wavelengths (0.5)
%
% X       : Received signals [M x N]
% A_true  : Array response at the true DoAs [M x K] (unknown)
% S       : Transmitted source signals [K x N] (unknown)
% Noise   : Additive noise [M x N] (unknown)
%
% theta   : Candidate angles [-90, 90] degrees
% A       : Steering vectors for all candidate angles
% R       : Estimated spatial covariance, R = X*X'/N

for i = 1:2
    %% 1. Load data and construct the array model
    data = load(files{i});
    X = data.X;                        % M antennas x N snapshots
    [M, N] = size(X);
    K = data.num_sources;     % Number of sources (= 2)
    Delta = data.Delta;       % Antenna spacing in wavelengths

    m = (0:M-1).';                     % Antenna indices
    A = exp(1j*2*pi*Delta * (m * sind(theta))); % Steering vectors
    R = X*X'/N;                       % Sample spatial covariance

    %% 2. Bartlett (classical beamformer)
    P_bartlett = real(sum(conj(A) .* (R*A), 1)) / M;

    %% 3. MVDR (Capon)
    P_mvdr = 1 ./ real(sum(conj(A) .* (R\A), 1));

    %% 4. MUSIC
    [V, D] = eig(R);
    [~, idx] = sort(real(diag(D)), 'ascend');
    E_noise = V(:, idx(1:M-K));       % Eigenvectors of noise subspace
    P_music = 1 ./ sum(abs(E_noise' * A).^2, 1);

    %% 5. Find peaks and compare spectra
    spectra = {P_bartlett, P_mvdr, P_music};
    
    subplot(2, 1, i);
    hold on;
    fprintf('\n%s\n', files{i});
    
    for j = 1:3
        P = spectra{j};
    
        % Detect peaks above 10% of the maximum
        [peakValues, peakIdx] = findpeaks(P,'MinPeakHeight', 0.10*max(P));
    
        % Select up to K strongest peaks
        [~, order] = sort(peakValues, 'descend');
        peakIdx = peakIdx(order(1:min(K, numel(order))));
    
        estimates = sort(theta(peakIdx));
    
        % Print estimated angles
        fprintf('%-9s: ', methods{j});
        fprintf('%.2f  ', estimates);
        fprintf('degrees\n');
    
        % Normalize spectrum to 0 dB
        P_dB = 10*log10(max(P/max(P), 1e-12));
        plot(theta, P_dB, 'LineWidth', 1.2);
    end
    
    title(strrep(files{i}, '.mat', ''));
    xlabel('Angle [degrees]');
    ylabel('Normalized spectrum [dB]');
    legend(methods, 'Location', 'best');
    grid on;
    xlim([-90 90]);
    ylim([-60 2]);
    hold off;
    end