
iters = 1100;
t = zeros(3,iters);
for n = 1:iters
    A = randi(10,5,3);
    y = randi(10,5,1);
    Sigma = randi(10,5,5);
    Sigma = Sigma'*Sigma;
    tic;
    house_gls(A, Sigma, y);
    t(1,n) = toc;
    tic;
    givens_gls(A, Sigma, y);
    t(2,n) = toc;
    tic;
    partial_revchol_gls(A, Sigma, y);
    t(3,n) = toc;
    
end
t1 = zeros(1,1000);
t2 = zeros(1,1000);
t3 = zeros(1,1000);
for i=101:1100
    t1(i-100)=t(1,i);
    t2(i-100)=t(2,i);
    t3(i-100)=t(3,i);
end
plot(t1);
hold on
plot(t2);
plot(t3);
title('Computation times of different GLS implementations per trial');
xlabel('trial');
ylabel('time(s)');
legend('householder', 'givens', 'reverse cholesky');
hold off