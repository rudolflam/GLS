function [ x_hat ] = house_gls( A, Sigma, y )
    %HOUSE_GLS GLS using householder transformations
    % [x_hat] = house_gls( A, Sigma, y ) returns
    % x_hat, the estimate of y = A x + v ~(0, Sigma) 
    
    B = chol(Sigma, 'lower');
    x_hat = house_gls_(A,B,y);
end

