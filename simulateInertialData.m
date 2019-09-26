function data = simulateInertialData(options)
% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

%% Simulate data for rotating sensor
nSamples = 800;  
ni = nSamples / 4;
T = 1/10;
rotSpeed = 2*pi/ni/T; % rotate around each axis in 100 samples
gyr = [zeros(ni,3) ; rotSpeed*ones(ni,1), zeros(ni,2) ; ...
    zeros(ni,1), rotSpeed*ones(ni,1), zeros(ni,1) ; ...
    zeros(ni,2), rotSpeed*ones(ni,1)]; 

% Make data set that is 10 times longer
gyr = repmat(gyr,[10 1]);

% Compute orientations based on angular velocities
nSamples = length(gyr);
q_init = [1;0;0;0];
qnb = zeros(nSamples,4);
qnb(1,:) = q_init';
for iSample = 2:nSamples+1
    qnb(iSample,:) = qRight( expq(T/2*gyr(iSample-1,:))) * qnb(iSample-1,:)';
end

%% Make accelerometer and magnetometer data
g = [0;0;-9.82];
% d = -67*pi/180; % Dip angle the Netherlands
d = 0*pi/180; % Dip angle equator
mn = [cos(d) ; 0 ; sin(d)];
acc = zeros(nSamples,3);
mag = zeros(nSamples,3);
for iSample = 1:nSamples
    acc(iSample,:) = - quat2rmat(qnb(iSample,:))' * g;
    mag(iSample,:) = quat2rmat(qnb(iSample,:))' * mn;
end

%% Add noise and disturb magfield
sigmaGyr = 5*pi/180;
sigmaAcc = -g(3)/100;
sigmaMag = 0.01;
gyr = gyr + sigmaGyr*randn(nSamples,3);
acc = acc + sigmaAcc*randn(nSamples,3);
mag = mag + sigmaMag*randn(nSamples,3);

if options.outliers
    threshold = 1-options.percOutliers/100;
    magnOutliers = options.magnOutliers;
    % Add outliers with same relative magnitude
    u = rand(nSamples,1);
    acc(u > threshold,:) = acc(u > threshold,:) + ...
        magnOutliers*g(3)*randn(length(acc(u > threshold,1)),3);
    u = rand(nSamples,1);
    mag(u > threshold,:) = mag(u > threshold,:) + ...
        magnOutliers*randn(length(mag(u > threshold,1)),3);
end

%% Save data
t = (0:length(gyr)-1)'*T;
accGyrMag = [acc,gyr,mag];
groundTruth.qnb = qnb;

settings.sigmaAcc = sigmaAcc*ones(3,1);
settings.sigmaGyr = sigmaGyr*ones(3,1);
settings.sigmaMag = sigmaMag*ones(3,1);
settings.g = g;
settings.mn = mn;
settings.T = T;

data.accGyrMag = accGyrMag;
data.settings = settings;
data.groundTruth = groundTruth; 

end