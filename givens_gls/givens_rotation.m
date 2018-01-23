function [ A, G ] = givens_rotation( A, j, xi, yi )
    %GIVENS Perform givens rotation on column j
    % [A,G] = givens_rotation( A, xi, yi, j ) returns
    % A, The matrix after applying the rotation zeroing A(yi,j)
    % G, The givens rotation
    % basd on C.C.Paige. Fast numerically stable computations for
    % generalized linear least squares problems   
    G = givens(A(xi,j),A(yi,j));
    A([xi,yi],:) = G * A([xi,yi],:);
    if nargout==2
        I = eye(size(A,1));
        I([xi,yi],:) = G' * I([xi,yi],:);
        G = I;
    end
end

