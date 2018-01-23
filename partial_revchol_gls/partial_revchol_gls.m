function [ x_hat ] = partial_revchol_gls( A, Sigma, y )
    %PARTIAL_REVCHOL_GLS GLS using reverse partial cholesky on Q'SigmaQ
    % where A = QR 
    % [x_hat] = partial_revchol_gls( A, Sigma, y ) returns
    % x_hat, the estimate of y = A x + v ~(0, Sigma) 
    
    [~,n] = size(A);
    [U, R] = house_qr_partial(A);
    QtS = house_apply_transpose(U, Sigma);
    QtSQ = house_apply_transpose(U, QtS')';
    
    U2cols = colrevchol(QtSQ,n+1);
    
    Qty = house_apply_transpose(U,y);
    x_z = [R U2cols(:,n+1:end)] \ Qty;
    x_hat = x_z(1:n);
end

