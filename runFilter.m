%% Code to run Monte Carlo simulations to test how the CF, MCF and MEKF w
addpath('tools')
rng(1) % For reproducibility

% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

%% Do 100 MC simulations and pre-allocate space for errors
nMC = 100;
e_MC = zeros(nMC,3);

for iMC = 1:nMC
    % Simulate data with or without outliers
    options.outliers = 1;
    options.percOutliers = 5;
    options.magnOutliers = 1;
    data = simulateInertialData(options);
    settings = data.settings;
    accGyrMag = data.accGyrMag;
    
    % Settings for filter:
    % Flag to estimate gyroscope bias as in Madgwick
    settings.estGyrBias = 0; 
    % Flag to estimate the local magnetic field as in
    settings.estimateMagneticField = 0;
    settings.init_q_nb = [1;0;0;0]; % Initial orientation
        
    % Run filter
    q = oriEst(accGyrMag, settings);
    
    % Compute RMSE in Euler angles
    e = quat2euler( qMult( q, qInv(data.groundTruth.qnb) ));
    e_MC(iMC,:) = rms(e);
end

%% Compute resulting RMSE over all Monte Carlo simulations
sqrt(mean(e_MC.^2))

