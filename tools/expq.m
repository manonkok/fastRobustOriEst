function eq = expq(phi)
% eq = expq(phi) -- Computes the quaternion exponential of the vector phi
% Assumes that a single 3-dimensional vector is used as input argument
% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

if any(size(phi) == 1)
    mag_phi = norm(phi);
    norm_phi = phi./(mag_phi + (mag_phi == 0));
    if size(norm_phi) == [1 3]
        norm_phi = norm_phi';
    end

    eq = [cos(mag_phi) ; norm_phi*sin(mag_phi)];
    if eq(1) < 0
        eq = - eq;
    end
else
    mag_phi = sqrt(phi(:,1).^2 + phi(:,2).^2 + phi(:,3).^2);
    norm_phi = phi./(mag_phi + (mag_phi == 0));
    eq = [cos(mag_phi) , norm_phi*sin(mag_phi)];
    eq(eq(:,1)<=0,:) = -eq(eq(:,1)<=0,:);
end