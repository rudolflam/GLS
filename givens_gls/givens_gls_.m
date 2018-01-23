function [ x_hat , mu_hat] = givens_gls_( A, B, y )
    %GIVENS_GLS GLS using givens rotations
    % [x_hat] = givens_gls( A, Sigma, y ) returns
    % x_hat, the estimate of y = A x + v ~(0, Sigma) 
    % mu_hat, the estimate for mean of noise distribution
    % basd on C.C.Paige. Fast numerically stable computations for
    % generalized linear least squares problems
    A_ = [y A B];
    n = size(A,2);
    m = size(B,2);
    Gt = eye(m);
    P = eye(m+n+1);
    for i=1:n+1
       yi = i;
       for j=n+1:-1:n+2-i
          [A_,Q] = givens_rotation(A_, j, yi+1, yi);
          %Gt = Q' * Gt;

          [A_,Q] = givens_rotation_right(A_, yi, n+yi+1,n+yi+2);
          %P = P*Q';
          yi = yi-1;
       end
    end
    r = A_(:, size(A_,2)-n);
    r = r(end-n:end,:);
    L = [r A_(2:end, 2:1+n)];
    z = A_(2:end,1);
    hat = L\z;
    mu_hat = hat(1,:);
    x_hat = hat(2:end, :);
end

