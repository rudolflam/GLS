function [R,Q]=house_rq(A)
    % Householder reflections for RQ decomposition using QR.
    % [R,Q] = house_rq(A) returns
    % R, the upper triangular factor, and
    % Q, the orthogonal matrix.   
    R = flipud(A)';
    
    H = @(u,x) x - u*(u'*x);
    [m,n] = size(R);
    U = zeros(m,n);
    
    for j = 1:min(m,n)
        u = house_gen(R(j:m,j));
        U(j:m,j) = u;
        R(j:m,j:n) = H(u,R(j:m,j:n));
        R(j+1:m,j) = 0;
    end
    Q = house_apply(U, eye(min(m,n)));
    
    Q = flipud(Q')';
    R = flipud(flipud(R)');
end
