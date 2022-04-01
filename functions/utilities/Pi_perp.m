function P = Pi_perp(A)
    % P = PI_PERP(A) returns the projection operator that projects the row
    % space of a matrix onto the row space of the orthogonal complement of A
        P = eye(size(A,2))-A'*pinv(A*A')*A;
