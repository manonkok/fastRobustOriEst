function pL = qLeft(q)
% Computes the left multiplication matrix of a quaternions.
% Copyright (C) 2019 by Manon Kok and Thomas B. Schon.

% Determine left matrix multiplication matrix
if any(size(q) == 1)
    if size(q) == [1 4]
        q = q';
    end
    pL = [q(1) , -q(2:4)' ; q(2:4) , q(1)*eye(3)+matrixCross(q(2:4))];
else
    N = size(q,1); % Determine N
    for i = 1:N
        pL(:,:,i) = [q(i,1) , -q(i,2:4) ; q(i,2:4)', ...
            q(i,1)*eye(3)+matrixCross(q(i,2:4)')];
    end
end