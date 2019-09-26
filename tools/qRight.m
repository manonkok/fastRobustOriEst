function qR = qRight(q)
% Computes the right multiplication matrices of a quaternion.
% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

if any(size(q) == 1)
    if size(q,1) == 1
        q = q';
    end
    qR = [q(1) , -q(2:4)' ; q(2:4) , q(1)*eye(3)-matrixCross(q(2:4))];
else
    N = size(q,1); % Determine N
    for i = 1:N
        qR(:,:,i) = [q(i,1) , -q(i,2:4) ; q(i,2:4)', ...
            q(i,1)*eye(3)-matrixCross(q(i,2:4)')];
    end
end