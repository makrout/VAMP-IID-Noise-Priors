close all;
clear all;
clc;

% add folder paths
addpath('priors')  
addpath('lmmse') 

%% initialize parameters
N = 200;
m = 100;
Nbr_iter = 150;
damp = 0.8;
fprintf(1,"--> Problem dimensions:\n N: %d \n m: %d \n", N,m);

%% signal and noise priors
% two priors supported: 'bg' (Bernoulli-Gaussian) and 'binary'
% choose each prior by commenting the undesired one

% prior for x: choose between 'binary' and 'bg'
%prior_x.name = 'binary';
prior_x.name = 'bg';

% prior for w: choose between 'binary' and 'bg'
%prior_w.name = 'binary';
prior_w.name = 'bg';

% generate the noise w according to its prior
if strcmp(prior_w.name, 'bg')
    prior_w.rho = 0.05; % percentage of non-0
    w = zeros(m, 1);
    n_nonzeros_ = ceil(prior_w.rho * m);
    non_zero_indices_ = randperm(m, n_nonzeros_);
    w(non_zero_indices_) = randn(n_nonzeros_, 1);
elseif strcmp(prior_w.name, 'binary')
    w = 2*randi(2,m,1)-3; 
end

% generate the signal x according to its prior
if strcmp(prior_x.name, 'bg')
    % generate a sparse vector x with sparsity 1-rho
    prior_x.rho = 0.05; % percentage of non-0
    x = zeros(N, 1);
    n_nonzeros = ceil(prior_x.rho * N);
    non_zero_indices = randperm(N, n_nonzeros);
    x(non_zero_indices) = randn(n_nonzeros, 1);

elseif strcmp(prior_x.name, 'binary')
    % generate a binary vector x {-1,1}    
    x = 2*randi(2,N,1)-3;
else
    disp('unknown prior');
end

fprintf(1,"--> Priors:\n signal: %s \n noise: %s \n", prior_x.name, prior_w.name);

% generate a sensing matrix A
A = randn(m,N);

% construct the observation
y = A*x + w;
snr = 10*log10(norm(A*x)^2/norm(w)^2);

% run VAMP with i.i.d. noise priors
[x_est_VAMP, nrmses] = VampNoiseIID(A, y, x, prior_x, prior_w, Nbr_iter, damp);

if strcmp(prior_x.name, 'binary')
    x_est_VAMP = round(x);
end
nrmse = sqrt(mean((A*x - A*x_est_VAMP).^2/mean((A*x).^2), 'all'));
fprintf(1,'\n Final nrmse = %f \n', nrmse);

%% Plotting
% the true signal vs the estimate
figure(1)
hold on
stem(x,'b')
stem(x_est_VAMP,'--r')
grid on;
legend('true', 'estimated')
title('True vs. estimated signals');

% the nrmse over iterations
figure(2)
loglog(1:length(nrmses), nrmses)
grid on;
xlabel('Iterations')
ylabel('NRMSE')
title('NRMSE vs Iterations')