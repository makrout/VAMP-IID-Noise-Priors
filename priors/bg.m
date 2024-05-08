function [x_hat, alpha] = bg(r_hat,gamma,N, signal_var,rho)

    mu_1_r = r_hat;
    N_1_r = 2/gamma;

    
    mu_2_x = zeros(1,N);
    N_2_x= 2*signal_var;
    
    mu_r_x = (N_2_x*mu_1_r + N_1_r*mu_2_x)/(N_1_r + N_2_x);            
    N_r_x = (N_1_r*N_2_x)/(N_1_r + N_2_x);
    Ax_1 = (rho*mu_r_x./(sqrt(N_1_r + N_2_x))).*exp(-(abs(mu_1_r - mu_2_x).^2)/(N_1_r + N_2_x));
    
    Ax_2 = (((N_r_x/2 + mu_r_x.^2)*rho)./sqrt(N_1_r + N_2_x)).*...
                                                   exp(-(abs(mu_1_r - mu_2_x).^2)./(N_1_r + N_2_x));
                                                                                     
    Bx_1 = ((1 - rho)/sqrt(N_1_r)).*exp(-(abs(mu_1_r).^2)/(N_1_r)) + ...
                                (rho/sqrt(N_1_r + N_2_x)).*exp(-(abs(mu_1_r - mu_2_x).^2)/(N_1_r + N_2_x));
   
    x_hat = Ax_1./Bx_1;
    
    Q_x = Ax_2./Bx_1 - x_hat.^2;
        
    alpha = gamma*mean(Q_x');
end