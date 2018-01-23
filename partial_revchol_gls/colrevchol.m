function [ R ] = colrevchol( S, t )
    %COLREVCHOL column wise reverse cholesky factorization
    % [ R ] = colrevchol( A, t ) returns
    % R, The (right most k columns of) upper triangular factor of A = R*R'
    k=3;
    n = size(S,1);
    L = zeros(size(S));
    for j=n:-1:t
        S(j,j) = sqrt(S(j,j));
        L(j,j) = S(j,j);
        a = S(j,j);
        S(j,1:j-1) = S(j,1:j-1)/a;
        L(j,1:j-1) = S(j,1:j-1);
        for k=1:j-1
           for i=k:j-1
               S(i,k)=S(i,k) - S(j,k)*S(j,i) ;
           end
        end
    end
    R=L';
end

