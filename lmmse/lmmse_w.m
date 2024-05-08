function [w2_hat,alpha2] = lmmse_w(r2_hat, gamma2, q2_hat, delta2, A, y)

    r2_hat = r2_hat';
    q2_hat = q2_hat';    
    
    m = size(A,1);  

    Q = inv(A*A');
    term1 = inv(delta2*eye(m) + gamma2*Q);
    term2 = (delta2*q2_hat + gamma2 * Q * (y-A*r2_hat));
    w2_hat =  term1 * term2;
    alpha2 = (delta2/m) * trace(term1);
    w2_hat = reshape(w2_hat, [1,m]);
end