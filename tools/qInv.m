function q = qInv(q)
% Computes the quaternion inverse of q
% Assumes that q = Nx4. For a single quaternion, both 1x4 and 4x1 are
% accepted.
% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

if any(size(q) == 1)
    q(2:4) = -q(2:4);
else
    q = [q(:,1) -q(:,2) -q(:,3) -q(:,4)];
end

end