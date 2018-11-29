%% double integrator system definition
k = -0.5;  % Friction
delta_t = 1;  % 1s sampling time
A = [1 delta_t;
    0 1 - k];
b = [0;
    delta_t];

C = [1 0];
%% system observer
L = place(A', C', [-0.5 -0.1]);
L = L';
disp(L);
disp(A - L*C);

%% Testing just normal MPC with Gaussian disturbance
xL = [-10; -3];
xU = [10; 10];
uL = -5;
uU = 5;

M = 40;
N = 5;
sigma = 0.5;
Q = [1 0; 
    0 1];
P = Q;
R = 1;

x0 = [5; 0];
x0hat = [5; 0];
[feas, xOpt, uOpt, xhat, predErr] = MPC_stochastic(A, b, C,L, P,...
                                                   x0, x0hat, M, N, sigma, Q, R,...
                                                   xL, xU, uL, uU);
plot(xOpt(1,:));

%% Testing unscented MPC
xL = [-10; -3];
xU = [10; 5];
uL = -5;
uU = 5;

M = 10;
N = 5;
Q = [5 0; 
    0 5];
P = Q;
R = 1;

cov_f = [0.1 0; 
        0 0.1];
cov_h = 0.2;

x0 = [5; 0];
x0hat = [5; 0];
x0cov = zeros(2);

[feas, xOpt, uOpt, xhat, predErr] = UKF_MPC_linear(A, b, C, P, x0, x0hat, x0cov, M, N,...
                                                      Q, R, xL, xU, uL, uU, cov_f, cov_h);                                                 