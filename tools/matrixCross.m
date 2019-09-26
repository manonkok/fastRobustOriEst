function M = matrixCross(v)
% Computes  matrix cross product of vector v: M = [v x]
% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

M = [0,-v(3),v(2);...
    v(3),0,-v(1);...
    -v(2),v(1),0];

end

