function [x2_hat,alpha2] = lmmse_x(r2_hat, gamma2, q2_hat, delta2, A, y)
    r2_hat = r2_hat';
    q2_hat = q2_hat';
    
    [m, N] = size(A);
    
    
    %Sigma = gamma2*A*A' + delta2*eye(m);
    %invSigma = inv(Sigma);
    %Gauss = mvnpdf(y-q2_hat,A*r2_hat,Sigma);
    %Gauss = 1;
    
    term1 = inv(gamma2*eye(N) + delta2*A'*A);
    term2 = (gamma2*r2_hat + delta2 * A' * (y-q2_hat));

    x2_hat =  term1 * term2;
    alpha2 = (gamma2/N) * trace(term1);
    x2_hat = reshape(x2_hat, [1,N]);

    %term1 = (r2_hat + gamma2*A'*invSigma*(y-q2_hat-A*r2_hat));
    %x2_hat = Gauss * term1;
    %dx2dr2 = Gauss * (((term1 - r2_hat)/gamma2)*term1' + ...
    %                   eye(N) - gamma2*A'*invSigma*A);    
    %alpha2 = (1/N) * trace(dx2dr2);
    %x2_hat = reshape(x2_hat, [1,N]);
end