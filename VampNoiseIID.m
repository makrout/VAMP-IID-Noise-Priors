function [x_est, nrmses] = VampNoiseIID(A, y, x, prior_x, prior_w, Nbr_iter, damp)

    % inline definition of the clip function
    clip = @(x,a,b) min(max(x,a),b);

    [m, N] =  size(A);

    % variables for the estimation of x
    r1_hat = zeros(1,N);
    x1_hat = zeros(1,N);
    gamma1 = 1;

    r1_hat_old = zeros(1,N);
    gamma1_old = 1;
    r2_hat_old = zeros(1,N);
    gamma2_old = 1;


    % variables for the estimation of w
    q1_hat = zeros(1,m);
    delta1 = 1;

    q1_hat_old = zeros(1,m);
    delta1_old = 1;
    q2_hat_old = zeros(1,m);
    delta2_old = 1;

    % initialize convergence criterion
    diff_old = inf;
    
    nrmses = [];
    
    for t = 1:Nbr_iter
        %% x estimation
        % denoising
        if strcmp(prior_x.name, 'bg')
            [x1_hat, alpha1] = bg(r1_hat,gamma1,N, 1, prior_x.rho);
        elseif strcmp(prior_x.name, 'binary')
            [x1_hat, alpha1] = binary(r1_hat,gamma1);
        end
        % extrinsic information computation
        eta1 = gamma1/alpha1;
        gamma2 = eta1 - gamma1;
        gamma2 = clip(gamma2, 10^(-11),10^(11));
        r2_hat = (eta1*x1_hat - gamma1*r1_hat)/gamma2;
        
        %% damping
        r2_hat = ((1-damp)*r2_hat_old*diag(gamma2_old) ...
                     + damp*r2_hat*diag(gamma2)) ...
                        * diag(1./((1-damp)*gamma2_old + damp*gamma2));
        r2_hat_old = r2_hat;
        gamma2 = (1-damp)*gamma2_old + damp*gamma2;
        gamma2_old = gamma2;

        %% w estimation
        % denoising
        if strcmp(prior_w.name, 'bg')
            [w1_hat, alpha1_] = bg(q1_hat,delta1,m, 1, prior_w.rho);
        elseif strcmp(prior_w.name, 'binary')
            [w1_hat, alpha1_] = binary(q1_hat,delta1);
        end
        % extrinsic information computation
        beta1 = delta1/alpha1_;
        delta2 = beta1 - delta1;
        delta2 = clip(delta2, 10^(-11),10^(11));
        q2_hat = (beta1*w1_hat - delta1*q1_hat)/delta2;
        
        %% damping
        q2_hat = ((1-damp)*q2_hat_old*diag(delta2_old) ...
                     + damp*q2_hat*diag(delta2)) ...
                        * diag(1./((1-damp)*delta2_old + damp*delta2));
        q2_hat_old = q2_hat;
        delta2 = (1-damp)*delta2_old + damp*delta2;
        delta2_old = delta2;

        %% x LLMSE step
        [x2_hat, alpha2] = lmmse_x(r2_hat, gamma2, q2_hat, delta2, A, y);
        eta2 = gamma2/alpha2;
        gamma1 = eta2-gamma2;
        gamma1 = clip(gamma1, 10^(-11),10^(11));
        r1_hat = (eta2*x2_hat - gamma2*r2_hat)/gamma1;
        x_est = reshape(x1_hat, [N,1]);

        %% damping
        r1_hat = ((1-damp)*r1_hat_old*diag(gamma1_old) ...
                     + damp*r1_hat*diag(gamma1)) ...
                        * diag(1./((1-damp)*gamma1_old + damp*gamma1));
        r1_hat_old = r1_hat;
        gamma1 = (1-damp)*gamma1_old + damp*gamma1;
        gamma1_old = gamma1;

        %% w LLMSE step
        [w2_hat, alpha2_] = lmmse_w(r2_hat, gamma2, q2_hat, delta2, A, y);
        beta2 = delta2/alpha2_;
        delta1 = beta2-delta2;
        delta1 = clip(delta1, 10^(-11),10^(11));
        q1_hat = (beta2*w2_hat - delta2*q2_hat)/delta1;


        %% damping
        q1_hat = ((1-damp)*q1_hat_old*diag(delta1_old) ...
                     + damp*q1_hat*diag(delta1)) ...
                        * diag(1./((1-damp)*delta1_old + damp*delta1));
                    
        % compute stopping criterion before updating the variables
        diff = mean(abs(q1_hat-q1_hat_old), 'all');
        
        q1_hat_old = q1_hat;
        delta1 = (1-damp)*delta1_old + damp*delta1;
        delta1_old = delta1;

        % compute error
        nrmse = sqrt(mean((A*x - A*x_est).^2/mean((A*x).^2), 'all'));
        nrmses = [nrmses, nrmse];
        
        % break the algoithm when the update between two iterations increase
        if (t>1) && (diff > diff_old)
            break;
        end
        diff_old = diff;
    end
end