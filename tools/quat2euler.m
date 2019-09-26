function e = quat2euler(q)
% e = quat2euler(q) -- Converts quaternions into Euler angles
% Assumes that q = Nx4. For a single quaternion, both 1x4 and 4x1 are
% accepted.
% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

if any(size(q) == 1)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
else
    q0 = q(:,1); q1 = q(:,2); q2 = q(:,3); q3 = q(:,4);
end

e = 180/pi*[atan2( (2*q2.*q3 - 2*q0.*q1) , (2*q0.^2 + 2*q3.^2 - 1 ) ) , ...
    - real(asin(2*q1.*q3 + 2*q0.*q2)) , ...
    atan2( (2*q1.*q2 - 2*q0.*q3) , (2*q0.^2 + 2*q1.^2 -1 ) ) ];
end
