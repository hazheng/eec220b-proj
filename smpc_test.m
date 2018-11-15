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
    0 0];
P = Q;
R = 1;

x0 = [5; 0];
x0hat = [5; 0];

[feas, xOpt, uOpt, xhat, predErr] = MPC_stochastic(A, b, C, L, P,...
                                                   x0, x0hat, M, N, sigma, Q, R,...
                                                   xL, xU, uL, uU);
plot(xOpt(1,:));

%% Testing unscented MPC

%% function definitions
function [feas, xOpt, uOpt, JOpt] = solve_cftoc(A, B, P, Q, R, N, x0, ...
                                          xL, xU, uL, uU, bf, Af, target_x)
      nu = size(B, 2);
      nx = size(A, 2);
      x = sdpvar(nx, N+1);  % include x0 in this
      u = sdpvar(nu, N);
      
      cost = (x(:,N+1) - target_x)' * P * (x(:,N+1) - target_x);
      for i =[1:N]
          cost = cost + (x(:,i) - target_x)'*Q*(x(:,i)-target_x) + u(:,i)'*R*u(:,i);
          % cost = cost + u(:,i)' * R * u(:,i);
      end
      constraints = [xL <= x(:,1) <= xU; 
          uL <= u(:,1) <= uU];
      index = 2;
      for i = 1:N-1
          constraints = vertcat(constraints, x(:,index) == A*x(:,index-1) + B*u(:,index-1),...
          xL <= x(:,index) <= xU,...
          uL <= u(:,index) <= uU);
          index = index + 1;
      end
      constraints = vertcat(constraints, x(:,N+1) == A * x(:,N) + B * u(:,N));
      constraints = vertcat(constraints, x(:,1) == x0);
      if isempty(Af)
          constraints = vertcat(constraints, x(:,N+1) == bf);
      else
          constraints = vertcat(constraints, Af * x(:,N+1) <= bf);
      end
      
      options = sdpsettings('verbose', 0);
      diag = optimize(constraints, cost, options);
      if diag.problem == 0
          feas = 1;
          xOpt = value(x);
          uOpt = value(u);
          JOpt = value(cost);
      else
          feas = 0;
          xOpt = [];
          uOpt = [];
          JOpt = inf;
      end
end


% make a conservative choice of terminal set: basically, the terminal 
% set for each iteration is the set that you need to reach that you could
% get to the final point going the minimum possible speed; ie. at time 20s,
% the terminal set is the set at which you could reach the terminal point
% going at (distance_to_final_point)/(tf - 20) m/s, assuming no
% disturbances (not persistently feasible, but provides a reasonable
% terminal set). 
function [feas, xOpt, uOpt, xhat, predErr] = MPC_stochastic(A, b, C, L, P, x0, x0hat, M, N, sigma,...
                                                      Q, R, xL, xU, uL, uU)
    predErr = zeros(2, M-N+1);
    pred_trajs = zeros(2, N, M);
    feas = zeros(1, M);
    xOpt = zeros(2, M+1);
    
    % observer
    xhat = zeros(2, M+1);
    xhat(:,1) = x0hat;
    
    xOpt(:,1) = x0;
    uOpt = zeros(1, M);
    xf = 0;
    % M is tfinal
    for i=1:M
        tx = xhat(1,i) - (xOpt(1,i) - xf)/(M - i); % for now using exact state
        disp(tx);
        % disp(tx);
        Af = [1 0; 0 0]; 
        % Af = [];
        bf = [tx; 0];
        
        [f, xo, uo, jo] = solve_cftoc(A, b, P, Q, R, N, xOpt(:,i), ...  % for now using exact state
                                          xL, xU, uL, uU, bf, Af,...
                                          xf);
        if f == 0
            feas(i:end) = 0;
            return
        end
        
        feas(i) = 1;
        pred_trajs(:,:,i) = xo(:,2:end);
        uOpt(i) = uo(1);
        xOpt(:,i+1) = A * xOpt(:,i) + b * uOpt(i) + randn(2,1) * sigma;
        xhat(:,i+1) = A * xhat(:,i) + b * uOpt(i) + L * (C*xOpt(:,i) - C * xhat(:,i));
        if i > N && i <= M - N
            predErr(1,i-N) = norm(pred_trajs(1,:,i-N) - xOpt(1,i-N+1:i), 2);
            predErr(2,i-N) = norm(pred_trajs(2,:,i-N) - xOpt(2,i-N+1:i), 2);
        end
    end
end


function [feas, xOpt, uOpt, xhat, predErr] = UKF_MPC(A, b, C, L, P, x0, x0hat, M, N, sigma,...
                                                      Q, R, xL, xU, uL, uU)
    % Procedure:
    % At each timestep, make a measurement
    % Pass this measurement through the ukf to obtain the state estimate
    % pass this state estimate into the stochastic cftoc solver
    % stochastic MPC solver constraints are generated according to
    % https://arxiv.org/pdf/1709.01201.pdf, using the propagation 
    % of the covariances to generate chance constraints
    % uses a robust time horizon to avoid the covariance 'blowing up' 
    return;
end


% f is dynamics function
% does one predict step
% implementation of https://arxiv.org/pdf/1709.01201.pdf
function [mean_pred, cov_pred] = propagate_mean_cov(x, cov_x, f, u, nx, cov_f)
    alpha = 1e-3;
    beta = 2;
    kappa = 0;
    
    L = 2 * nx + 1;
    lambda = alpha^2 * (L + kappa) - L;
    wm0 = lambda/(L + lambda);
    wc0 = lambda/(L + lambda) + (1 - alpha^2 + beta);
    
    % Doing this with sdpvars might get hairy...
    sqrt_cov = chol((L + lambda) * cov_x, 'nocheck');
    
    sigma_pts = zeros(nx, 2*L+1);
    sigma_pts(:,1) = x;
    for i=2:L
        sigma_pts(:,i) = x + sqrt_cov(i-1,:)';
    end
    for i=L+1:2*L
        sigma_pts(:,i) = x + sqrt_cov(i-1-L,:)';
    end
    ws = repmat(1/(2 * (L + lambda)), 1, 2 * L);
    sigma_pts_prop = zeros(nx, 2*L+1);
    sigma_pts_prop(:,1) = f(sigma_pts(:,1), u);
    mean_pred = wm0 * sigma_pts_prop(:,1);
    cov_pred = cov_f + wc0 * (sigma_pts_prop(:,1)-mean_pred)*...
        (sigma_pts_prop(:,1)-mean_pred)';
    
    for i=2:2*L+1
        sigma_pts_prop(:,i) = f(sigma_pts(:,i), u);
        mean_pred = mean_pred + ws(i) * sigma_pts_prop(:,i);
        cov_pred = cov_pred + ws(i) * (sigma_pts_prop(:,i) - mean_pred)*...
            (sigma_pts_prop(:,i) - mean_pred)';
    end
    
end
