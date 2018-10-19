%% Double integrator playground for Information Theoretic MPC
% Model definitions and parameters 

A = [1 1; 0 1];
B = [0; 1];

f = @(x,v) A*x + B*v;

sigma_v = .5; %input covariance

%cost matricies
Q = [1 0; 0 1];
R = 10; 

N = 20; % number of timesteps to run the simulation
K = 1000; % number of samples
T = 20; % MPC time horizon
alpha = .1;
lambda = 1;
gamma = lambda*(1-alpha);
x0 = [10; 0];

%% generate the initial control sequence via QP 
[U_init] = solve_ftmpc(f, Q, R, x0, T, 2, 1);

%% plot and forward propagate U to double check
X_traj = zeros(2, T+1);
X_traj(:,1) = x0;

for i = 1:T
    X_traj(:, i+1) = f(X_traj(:,i), U_init(i));
end

figure
plot(X_traj(1,:));
hold on
plot(X_traj(2,:));
legend('x1', 'x2');

figure
plot(U_init);
%% Run the Sample Based Information Theoretic MPC
X_traj_it = zeros(2, T+1);
X_traj_it(:,1) = x0;
X_traj_mpc = zeros(2, T+1);
X_traj_mpc(:,1) = x0;

U = U_init;
U_itmpc = zeros(1,N);
U_mpc = zeros(1,N);
V_itmpc = zeros(1,N);
for i = 1:N
    x = X_traj_it(:,i);
    %generate a set of trajectory perturbations (Epsilon)
    Epsilon = normrnd(0, sigma_v^2, K, T);
    threshold = floor((1-alpha)*K);
    %generate trajectories to sample from 
    V = cat(1, repmat(U, threshold,1) + Epsilon(1:threshold,:), Epsilon(threshold+1:end,:));
    S = zeros(1, K);
    % TODO: compute the traj costs
    for k = 1:K
        x_temp = x;
        for t = 1:T
            x = f(x, V(k,t));
            S(k) = S(k) + x'*Q*x + gamma*U(t)*inv(sigma_v)*V(k,t);
        end
    end
    %compute the sample weights
    S = S - min(S);
    W = exp(-(1/lambda)*S);
    W = W / sum(W);
    %now apply the weights
    %currently using unsmoothed inputs
    U = U + sum(diag(W)*Epsilon,1);
    U_itmpc(i) = U(1);
    %leftshift
    U = circshift(U, -1);
    U(T) = 0;
    
    %forward propagate the model subject to actuation noise
    e = normrnd(0, sigma_v^2);
    v = U_itmpc(i) + e;
    V_itmpc(i) = v;
    X_traj_it(:,i+1) = f(X_traj_it(:,i), v);
    %in parallel run normal MPC
    [U_temp] = solve_ftmpc(f, Q, R, X_traj_mpc(:,i), T, 2, 1);
    U_mpc(i) = U_temp(1);
    X_traj_mpc(:,i+1) = f(X_traj_mpc(:,i), U_mpc(i)+e);
end

figure
plot(X_traj_it(1,:));
hold on
plot(X_traj_it(2,:));
legend('x1', 'x2');
title("Information Theoretic MPC")

figure
plot(X_traj_mpc(1,:));
hold on
plot(X_traj_mpc(2,:));
legend('x1', 'x2');
title("Standard MPC")

figure
plot(U_itmpc);
hold on
plot(U_mpc);
legend('Information Theoretic', 'Standard');