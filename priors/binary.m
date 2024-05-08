function [x1, alpha1] = binary(r1,gamma1)
    K = size(r1,1);

    x1 = tanh(r1.*repmat(gamma1,K,1));
    
    alpha1 = gamma1.*mean(1-x1.^2);
    
    alpha1 = max(alpha1,1e-6);    
end
