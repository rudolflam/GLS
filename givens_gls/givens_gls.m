function [ x_hat , mu_hat] = givens_gls( A, Sigma, y )
    %GIVENS_GLS GLS using givens rotations
    % [x_hat] = givens_gls( A, Sigma, y ) returns
    % x_hat, the estimate of y = A x + v ~(0, Sigma) 
    % mu_hat, the estimate for mean of noise distribution
    % basd on C.C.Paige. Fast numerically stable computations for
    % generalized linear least squares problems
    
    B = chol(Sigma,'lower');
    [x_hat, mu_hat] = givens_gls_(A, B, y);
end

