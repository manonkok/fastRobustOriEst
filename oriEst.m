function [q,meanTime] = oriEst(accGyrMag, settings)
% Runs the fast and robust orientation estimation algorithm
% using the accGyrMag data. Settings are used to pass on the 
% local earth magnetic field, sampling time and the initial orientation

% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

N = length(accGyrMag); % Number of data points
q = zeros(N,4); % Pre-allocate orientation estimates
q(1,:) = settings.init_q_nb; % Initial orientation
gyrBias = zeros(3,1); % Gyroscope bias
mn = settings.mn; % Local magnetic field vector

for i = 1:N 
    % Extract data
    acc = accGyrMag(i,1:3);
    gyr = accGyrMag(i,4:6);
    mag = accGyrMag(i,7:9);
    
    % Do filter update for this time step
    [q(i+1,:),gyrBias, mn] = filterUpdate(acc, gyr, mag, ...
        q(i,:), gyrBias, mn, settings);
end

end

function [q, gyrBias, mn] = filterUpdate(acc, gyr, mag, ...
    q, gyrBias, mn, settings)

    T = settings.T; % Sampling time
    
    % Local magnetic field vector 
    b_x = mn(1); b_z = mn(3);
    
    % Compute beta
    gyroMeasError = settings.sigmaGyr(1,1);
    beta = sqrt(3) * gyroMeasError;
    
    if settings.estGyrBias
        % Estimate gyroscope bias
        w_bx = gyrBias(1); w_by = gyrBias(2); w_bz = gyrBias(3);
        zeta = sqrt(3) * 1E-2;
    end

    %% Extract variables
    w_x = gyr(1); w_y = gyr(2); w_z = gyr(3); % Gyr data
    a_x = acc(1); a_y = acc(2); a_z = acc(3); % Acc data
    m_x = mag(1); m_y = mag(2); m_z = mag(3); % Mag data
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4); % Quaternion

    %% Auxiliary variables to avoid repeated calculations
    q0_half = 0.5 * q0;
    q1_half = 0.5 * q1;
    q2_half = 0.5 * q2;
    q3_half = 0.5 * q3;
    q0_two = 2 * q0;
    q1_two = 2 * q1;
    q2_two = 2 * q2;
    b_x_two = 2 * b_x;
    b_z_two = 2 * b_z;
    q0q2 = q0 * q2; 
    q1q3 = q1 * q3;

    %% Normalise the accelerometer measurement
    n = sqrt(a_x^2 + a_y^2 + a_z^2);
    a_x = a_x / n;
    a_y = a_y / n;
    a_z = a_z / n;

    %% Normalise the magnetometer measurement
    n = sqrt(m_x^2 + m_y^2 + m_z^2);
    m_x = m_x / n;
    m_y = m_y / n;
    m_z = m_z / n;

    %% Compute the objective function and Jacobian for acc
    fbara_1 = q1_two * q3 - q0_two * q2;
    fbara_2 = q0_two * q1 + q2_two * q3;
    fbara_3 = 1 - q1_two * q1 - q2_two * q2;

    %% Compute the objective function and Jacobian for mag
    fbarm_1 = b_x_two * (0.5 - q2^2 - q3^2) + b_z_two * (q1q3 - q0q2);
    fbarm_2 = b_x_two * (q1 * q2 - q0 * q3) + b_z_two * (q0 * q1 + q2 * q3) ;
    fbarm_3 = b_x_two * (q0q2 + q1q3) + b_z_two * (0.5 - q1^2 - q2^2) ;

    %% Compute the gradient (matrix multiplication)
    wHatAccMag_1 = fbara_2 * a_z - fbara_3 * a_y + fbarm_2 * m_z - fbarm_3 * m_y;
    wHatAccMag_2 = fbara_3 * a_x - fbara_1 * a_z + fbarm_3 * m_x - fbarm_1 * m_z;
    wHatAccMag_3 = fbara_1 * a_y - fbara_2 * a_x + fbarm_1 * m_y - fbarm_2 * m_x;

    %% Normalise the gradient
    n = sqrt(wHatAccMag_1^2 + wHatAccMag_2^2 + wHatAccMag_3^2);
    wHatAccMag_1 = wHatAccMag_1 / n;
    wHatAccMag_2 = wHatAccMag_2 / n;
    wHatAccMag_3 = wHatAccMag_3 / n;
    
    %% Compensate for gyr bias
    if settings.estGyrBias
        % compute and remove the gyroscope baises
        w_bx = w_bx + wHatAccMag_1 * T * zeta;
        w_by = w_by + wHatAccMag_2 * T * zeta;
        w_bz = w_bz + wHatAccMag_3 * T * zeta;
        w_x = w_x - w_bx;
        w_y = w_y - w_by;
        w_z = w_z - w_bz;
    end
    
    %% Compute the estimated angular velocity and compute S(q) w/2
    wHat_x = T * (w_x - beta * wHatAccMag_1);
    wHat_y = T * (w_y - beta * wHatAccMag_2);
    wHat_z = T * (w_z - beta * wHatAccMag_3);
    dq0 = -q1_half * wHat_x - q2_half * wHat_y - q3_half * wHat_z;
    dq1 = q0_half * wHat_x + q2_half * wHat_z - q3_half * wHat_y;
    dq2 = q0_half * wHat_y - q1_half * wHat_z + q3_half * wHat_x;
    dq3 = q0_half * wHat_z + q1_half * wHat_y - q2_half * wHat_x;

    %% Integrate the estimated quaternion derivative
    q0 = q0 + dq0; 
    q1 = q1 + dq1;
    q2 = q2 + dq2;
    q3 = q3 + dq3;

    %% Normalise quaternion
    n = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
    q0 = q0 / n; q1 = q1 / n; q2 = q2 / n; q3 = q3 / n;

    %% Return updated quaternion
    q = [q0 ; q1 ; q2 ; q3];
    
    if settings.estGyrBias
        % and optionally also return the estimated gyroscope bias
        gyrBias = [w_bx ; w_by ; w_bz];
    end
    
    %% Compute magnetic field in navigation frame
    if settings.estimateMagneticField
        m_x_two = 2 * m_x;
        m_y_two = 2 * m_y;
        m_z_two = 2 * m_z;
        % Auxiliary variables
        q0q1 = q0 * q1; 
        q0q2 = q0 * q2;
        q0q3 = q0 * q3;
        q2q3 = q2 * q3;
        q1q2 = q1 * q2;
        q2q3 = q1 * q3;
        
        % Compute local magnetic field
        h_x = m_x_two * (0.5 - q2^2 - q3^2) + ...
            m_y_two * (q1q2 - q0q3) + ...
            m_z_two * (q2q3 + q0q2);
        h_y = m_x_two * (q1q2 + q0q3) + ...
            m_y_two * (0.5 - q1^2 - q3^2) + ...
            m_z_two * (q2q3 - q0q1);
        h_z = m_x_two * (q2q3 - q0q2) + ...
            m_y_two * (q2q3 + q0q1) + ...
            m_z_two * (0.5 - q1^2 - q2^2);

        % and rotate to get the horizontal component only in the x-direction
        b_x = sqrt(h_x^2 + h_y^2);
        b_z = h_z;
        mn = [b_x, 0, b_z];
    end
end