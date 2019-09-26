function q = qMult(q1, q2) 
% Computes the quaternion multiplication of q1 and q2
% Assumes that q1 and q2 are Nx4. For single quaternions, 
% both 1x4 and 4x1 are accepted.
% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

if any(size(q1) == 1)
    if size(q2,1) == 4
        q = qLeft(q1) * q2;
    else
        q = qLeft(q1) * q2';
    end
else
    N = size(q2,1);
    q1L = qLeft(q1);
    q = zeros(N,4);
    for i = 1:N
        q(i,:) = q1L(:,:,i) * q2(i,:)';
    end
end