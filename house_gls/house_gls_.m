function [ x_hat ] = house_gls_( A, B, y )
    %HOUSE_GLS GLS using householder transformations
    % [x_hat] = house_gls( A, Sigma, y ) returns
    % x_hat, the estimate of y = A x + v ~(0, Sigma) where B*B' = Sigma
    [U, R] = house_qr_partial(A);
    QtB = house_apply_transpose(U, B);
    [S,~] = house_rq(QtB);
    [m,n] = size(R);
    %Q = house_apply(U, eye(m));
    Qty = house_apply_transpose(U,y);
    U2cols = S(:,  end-(m-n)+1:end);
    % [R U2cols] is a triangular matrix
    x_z = [R U2cols] \ Qty;
    x_hat = x_z(1:n);
end

