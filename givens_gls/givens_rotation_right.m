function [ A, G ] = givens_rotation_right( A, i, xj, yj )
    %GIVENS Perform givens rotation on row j from the right side of A
    % [A,G] = givens_rotation_right( A, xi, yi, j ) returns
    % A, The matrix after applying the rotation zeroing A(yi,j)
    % G, The givens rotation
    % basd on C.C.Paige. Fast numerically stable computations for
    % generalized linear least squares problems  
    [A,G] = givens_rotation(A',i, xj,yj);
    A = A'; G=G';
end

