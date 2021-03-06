%% Test out energy heuristics / cost models for the high level energy scheduler
% Mock strategy over the course of a day
% sim params 
N_day = 16;
N = N_day*3;
dt = 8/N_day; %time unit is in hours, 
soc_max = 5000; %units in W/h
pout_max = 50000; % units in W
% solar = zeros(N,1);
solar = 1400*sin(linspace(pi/6, 5*pi/6, N_day));
solar_cloudy = 400 * sin(linspace(pi/6, 5*pi/6, N_day));

solar = repmat(solar,1,6);
solar_cloudy = repmat(solar, 1, 6);

%account for evening and morning charge
energy_evening = 800*4;
solar(N_day+1) = solar(N_day+1) + energy_evening/dt;
solar(2*N_day+1) = solar(2*N_day+1)+ energy_evening/dt;

solar_cloudy(N_day+1) = solar_cloudy(N_day+1) + energy_evening/(2 * dt);
solar_cloudy(2*N_day+1) = solar_cloudy(2*N_day+1)+ energy_evening/(2 * dt);

soc_init = .95*soc_max;
%%
% compute the mean power draw
% state just consists of SOC 
soc = sdpvar(N+1, 1);
u = sdpvar(N,1); 
soc(1) = soc_init;

energy_total = soc_init + sum(solar)*dt;
power_avg = energy_total/24;

% assign(u, power_avg/pout_max*ones(N, 1));
c = [soc_max >= soc >= 0, 1 >= u >= 0];

for i = 1:N
    c = [c, soc(i+1) == soc(i) - dt*pout_max*u(i) + solar(i)*dt];
end

obj = -sum(u.^(1/3));

ops = sdpsettings('fmincon.maxiter', 10000);
optimize(c, obj, ops); 

u = double(u);
soc = double(soc);

figure
subplot(2,1,1)
plot(soc)
subplot(2,1,2)
plot(u*pout_max)
%% Testing the Unscented MPC on this global optimization
f = @(soc, u, i) soc - dt*pout_max*u + solar(i) * dt;
f_cloudy = @(soc, u, i) soc - dt*pout_max*u + solar_cloudy(i) * dt;

h = @(soc) soc; % assume fine measurement for now

cov_f = 500;  
cov_h = 100;  % measurement can be off by 100 Wh (2% of total SoC)
nx = 1;
nu = 1;

% unused, legacy parameters to make stuff fit
P = 1;
Q = 1;
R = 1;

x0cov = 10;  % could be off by 10W
x0hat = .9 * soc_max;
x0 = .9 * soc_max;

M = N_day * 3;  % test for 3 days simulation horizon
N = N_day * 3;

uL = 0;
uU = 1;

xU = soc_max;
xL = 0.05 * soc_max;
noise_seq = randn(M+1,1) * 0;
%%
[feas, xOpt, uOpt, xhat, pred_traslidjs, pred_u] = UKF_MPC_tv(nx, nu, f, f_cloudy, h, P, x0, x0hat, x0cov, M, N,...
                                      Q, R, xL, xU, uL, uU, cov_f, cov_h, noise_seq);
%%
times = linspace(1, 24, 48);                                  
figure(1);
% subplot(2, 1, 1);
% hold on;
% for i = 1:30
%     plot(times+i, pred_trajs(:,:,i));
% end
% title("predicted trajectories UKF MPC");
% hold off;
subplot(2, 1, 1);
plot(times, xOpt(1:end-1));
hold on;
title("xOpt UKF MPC");
subplot(2, 1, 2);
plot(times, uOpt*pout_max);
hold on;
title("uOpt UKF MPC");
xlabel("Hours of Race");
%% Testing the Naive MPC 
f = @(soc, u, i) soc - dt*pout_max*u + solar(i) * dt;
f_cloudy = @(soc, u, i) soc - dt*pout_max*u + solar_cloudy(i) * dt;

h = @(soc) soc; % assume fine measurement for now

cov_f = 500;
cov_h = 100;  % measurement can be off by 100 Wh (2% of total SoC)
nx = 1;
nu = 1;

% unused, legacy parameters to make stuff fit
P = 1;
Q = 1;
R = 1;

x0cov = 10;  % could be off by 10W
x0hat = .9 * soc_max;
x0 = .9 * soc_max;

M = N_day * 3;  % test for 3 days simulation horizon
N = N_day * 3;

uL = 0;
uU = 1;

xU = soc_max;
xL = 0.05 * soc_max;

[feas, xOpt_naive, uOpt_naive, pred_trajs_naive, pred_u_naive] = MPC_Naive(nx, nu, f, f_cloudy, h, P, x0, M, N,...
                                      Q, R, xL, xU, uL, uU, cov_f, noise_seq);
%%
times = linspace(1, 24, 48);                                  
figure(1);
% subplot(3, 1, 1);
% hold on;
% for i = 1:30
%     plot(times+i, pred_trajs_naive(:,:,i));
% end
% title("pred trajs naive mpc");
subplot(2, 1, 1);
plot(times, xOpt_naive(1:end-1));
title("State of Charge Over Race");
legend({"UKF", "Naive"});
subplot(2, 1, 2);
plot(times, uOpt_naive*pout_max);
legend({"UKF", "Naive"});
title("Throttle Input");

%% Saving this bs code for some reason
        % Fixing the transitions so can compare with Naive
%         if mod(i, 20) == 0
%             if cur_model == 1
%                 cur_model = 0;
%                 f_tv = @(x, u, time) f_cloudy(x, u, time+ i);
%             else
%                 cur_model = 1;
%                 f_tv = @(x, u, time) f(x, u, time + i);
%             end
%         elseif cur_model == 1
%             f_tv = @(x, u, time) f(x, u, time + i);
%         else
%             f_tv = @(x, u, time) f_cloudy(x, u, time+ i);
%         end
% 
%         prob = rand();
%         if prob > 0.8 && cur_model == 1  % switching to cloudy 
%             cur_model = 0;
%             f_tv = @(x, u, time) f_cloudy(x, u, time+ i);
%         elseif prob > 0.8 && cur_model == 0  % switching to sunny
%             cur_model = 1;
%             f_tv = @(x, u, time) f(x, u, time + i);  
%         elseif cur_model == 1
%             f_tv = @(x, u, time) f(x, u, time + i);
%         else
%             f_tv = @(x, u, time) f_cloudy(x, u, time+ i);
%         end
%% functions
% need to change the damn mpc formulation to account for time-varying
% solar, so here's all the functions again

% this assumes time-varying f (with constant covariance for now), but time-invariant h
function [feas, xOpt, uOpt, xhat, pred_trajs, pred_u] = UKF_MPC_tv(nx, nu, f, f_cloudy, h, P, x0, x0hat, x0cov, M, N,...
                                                      Q, R, xL, xU, uL, uU, cov_f, cov_h, f_noise_seq)
    tr = 10;
    H = [-1];
    
    g = [-xL];
    
    p = [0.95]; 
    
%     predErr = zeros(nx, M-N+1);
    pred_trajs = zeros(nx, N, M);
    pred_u = zeros(nu, N, M);
    feas = zeros(1, M);
    xOpt = zeros(nx, M+1);
    
    % UKF observer
    xhat = zeros(nx, M+1);
    xhat(:,1) = x0hat;
    xCovs = zeros(nx, nx,M+1);
    xCovs(:,:,1) = x0cov;

    
    xOpt(:,1) = x0;
    uOpt = zeros(1, M);
    xf = 0;
    cur_model = 1;  % currently using sunny model
    % M is tfinal
    for i=1:M
        disp(i);
        f_tv = @(x, u, time) f(x, u, time + i);
        
        horz_length = min(N, M + 1 - i);
        [f_prob, xo, uo, jo] = solve_ucftoc_timevarying(nx, nu, P, xhat(:,i), xCovs(:,:,i), horz_length,...
                                       Q, R, H, g, p, f_tv, cov_f, tr, uL, uU,...
                                       0);
%         disp(xo);
        % disp(uo(1));
        if f_prob == 0
            feas(i:end) = 0;
            return
        end
        
        feas(i) = 1;
        pred_trajs(:,:,i) = [xo(:,2:end) zeros(1, N - horz_length)];
        pred_u(:,:,i) = [uo zeros(1, N - horz_length)];
        uOpt(i) = uo(1);
        
        if isempty(f_noise_seq)
            noise_f = cov_f * randn(nx, 1);
        else
            noise_f = f_noise_seq(i);
        end
        
        noise_h = cov_h * randn();
        
        if cur_model == 0
            xOpt(:,i+1) = f_cloudy(xOpt(:,i), uOpt(i), i) + noise_f;  % only allowing additive noise for now
        else
            xOpt(:,i+1) = f(xOpt(:,i), uOpt(i), i) + noise_f;  % only allowing additive noise for now
        end
        [x_pred, cov_x_pred, sigma_pts_prop, wm0, wc0, ws] = propagate_mean_cov_tv(xhat(:,i), xCovs(:,:,i), f, uOpt(i), nx, cov_f, i);
        [y_pred, cov_y_pred, cov_xy_pred] = generate_output_prediction(sigma_pts_prop, x_pred, wm0, wc0, ws, h, cov_h);
        y_meas = h(xOpt(:,i+1)) + noise_h;
        [x_update, cov_update] = ukf_update(x_pred, y_pred, cov_x_pred, cov_y_pred, cov_xy_pred, y_meas);
       
        xhat(:,i+1) = x_update;
        xCovs(:,:,i+1) = cov_update;
    end    
end

function [feas, xOpt, uOpt, JOpt] = solve_ucftoc_timevarying(nx, nu, P, x0, x_cov, N,...
                                                 Q, R, H, g, p, f, cov_f, tr, uL, uU,...
                                                 use_constraints)
    nc = size(H, 1);  % number of rows of H = number of constraints per timestep
    
    x = sdpvar(nx, N+1);  % this is the mean, or expected state value
    u = sdpvar(nu, N);
    
    % have to do shenanigans with declaring arrays of sdpvars
    % because of https://yalmip.github.io/naninmodel/
    % covs_tr = zeros(nx, nx, tr+1);  % only track covs up to tr
    covs_tr = x_cov;
    
    cost = 0;
    for i = [1:N]
        cost = cost - u(i)^(1/3);
    end
    
    
    constraints = [uL <= u(:,1) <= uU];
    index = 2;
    for i = 1:N-1
      constraints = [constraints; uL <= u(:,index) <= uU];
      index = index + 1;
    end
    
    constraints = [constraints; x(:,1) == x0];
    for i=2:N+1
        [mean_pred, cov_pred, sigma_pts_prop, wm0, wc0, ws] = propagate_mean_cov_tv(x(:,i-1), covs_tr(:,:,i-1), f, u(:,i-1), nx, cov_f, i-1);
        constraints = [constraints; x(:,i) == mean_pred];
        if i <= tr
            % covs_tr(:,:,i) = cov_pred;
            covs_tr = cat(3, covs_tr, cov_pred);
%             cost = cost - (x(:,i) - 3 * sqrtm(cov_pred));
        else
            % covs_tr(:,:,i) = covs_tr(:,:,i-1);
            covs_tr = cat(3, covs_tr, covs_tr(:,:,i-1));
        end
        % expectation + variance padding formulation
        % cost = cost + 0.5 * exp(x(:,i) - 5000) + exp(500 + 0.5 * covs_tr(:,:,i) - x(:,i));
        
        % Expectation Formulation
        cost = cost + exp(0.01 * (x(:,i) - 5000)) + exp(0.01 * (-x(:,i)));
        
        % reformulation as a probabilistic constraint based on covariance
        % reformulation as per https://arxiv.org/pdf/1709.01201.pdf
        if use_constraints == 1
            for j=1:nc
                constraints = [constraints;
                    (quantile(p(j)) * sqrtm(H(j,:) * covs_tr(:,:,i) * H(j,:)')+...
                    H(j,:) * x(:,i)) <= g(j)];
            end     
        end
    end
        
    options = sdpsettings('verbose', 0, 'fmincon.maxiter', 1000000);
%     options = sdpsettings('verbose', 0);
    diag = optimize(constraints, cost, options);
%     diag = optimize(constraints, cost);
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

function [mean_pred, cov_pred, sigma_pts_prop, wm0, wc0, ws] = propagate_mean_cov_tv(x, cov_x, f_tv, u, nx, cov_f, time)
    alpha = 1e-3;
    beta = 2;
    kappa = 0;
    
    L = nx;
    lambda = alpha^2 * (L + kappa) - L;
    wm0 = lambda/(L + lambda);
    wc0 = lambda/(L + lambda) + (1 - alpha^2 + beta);
    
    % Doing this with sdpvars might get hairy...
    if all(class(cov_x) == 'double') && all(cov_x(:) == 0)
        sqrt_cov = zeros(nx);
    else
        sqrt_cov = cholesky((L + lambda) * cov_x);
    end
    
    % unfortunately you have issues with
    % assigning SDPVARS to arrays of zeros []. 
    % sigma_pts = zeros(nx, 2*L+1);
    sigma_pts = x;
    for i=2:L+1
        % sigma_pts(:,i) = x + sqrt_cov(i,:)';
        sigma_pts = [sigma_pts, x+sqrt_cov(i-1,:)'];
    end
    for i=L+2:(2*L + 1)
        % sigma_pts(:,i) = x + sqrt_cov(i-L,:)';
        sigma_pts = [sigma_pts, x-sqrt_cov(i-1-L,:)'];
    end
    ws = repmat(1/(2 * (L + lambda)), 1, 2 * L);
    % sigma_pts_prop = zeros(nx, 2*L+1);
    % sigma_pts_prop(:,1) = f(sigma_pts(:,1), u);
    sigma_pts_prop = f_tv(sigma_pts(:,1), u, time);
    
    mean_pred = wm0 * sigma_pts_prop(:,1);
    
    for i=2:2*L+1
        % sigma_pts_prop(:,i) = f(sigma_pts(:,i), u);
        sigma_pts_prop = [sigma_pts_prop, f_tv(sigma_pts(:,i),u, time)];
        mean_pred = mean_pred + ws(i-1) * sigma_pts_prop(:,i);
        
    end
    
    cov_pred = cov_f + wc0 * (sigma_pts_prop(:,1)-mean_pred)*...
        (sigma_pts_prop(:,1)-mean_pred)';
    for i = 2:2*L+1
        cov_pred = cov_pred + ws(i-1) * (sigma_pts_prop(:,i) - mean_pred)*...
            (sigma_pts_prop(:,i) - mean_pred)';
    end
end

function [feas, xOpt, uOpt, JOpt] = solve_cftoc_timevarying(nx, nu, x0, N,...
                                                            xL, xU, f, uL, uU)
    % nc = size(A, 1);  % number of rows of A = number of constraints per timestep
    
    x = sdpvar(nx, N+1);
    u = sdpvar(nu, N);
  
    cost = 0;
    for i = [1:N]
        cost = cost - u(i)^(1/3);
    end
    constraints = [];
    for i = 1:N
      constraints = [constraints; uL <= u(:,i) <= uU];
    end
    
    constraints = [constraints; x(:,1) == x0];
    for i=2:N+1
        constraints = [constraints; 
            x(:,i) == f(x(:,i-1), u(:,i-1), i-1);
            xL <= x(:,i) <= xU];
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

% assume perfect knowledge
function [feas, xOpt, uOpt, pred_trajs, pred_u] = MPC_Naive(nx, nu, f, f_cloudy, h, P, x0, M, N,...
                                                      Q, R, xL, xU, uL, uU, cov_f, f_noise_seq)
    pred_trajs = zeros(nx, N, M);
    pred_u = zeros(nu, N, M);
    feas = zeros(1, M);
    xOpt = zeros(nx, M+1);
        
    xOpt(:,1) = x0;
    uOpt = zeros(1, M);
    cur_model = 1;  % currently using sunny model
    % M is tfinal
    for i=1:M
        disp(i);
        f_tv = @(x, u, time) f(x, u, time + i);
        
        horz_length = min(N, M + 1 - i);
        [f_prob, xo, uo, jo] = solve_cftoc_timevarying(nx, nu, xOpt(:,i), horz_length, xL, xU, f_tv, uL, uU);
%         disp(xo);
        % disp(uo(1));
        if f_prob == 0
            feas(i:end) = 0;
            return
        end
        
        feas(i) = 1;
        pred_trajs(:,:,i) = [xo(:,2:end) zeros(1, N - horz_length)];
        pred_u(:,:,i) = [uo zeros(1, N - horz_length)];
        uOpt(i) = uo(1);
        if isempty(f_noise_seq)
            noise_f = cov_f * randn(nx, 1);
        else
            noise_f = f_noise_seq(i);
        end
        
        if cur_model == 0
            xOpt(:,i+1) = f_cloudy(xOpt(:,i), uOpt(i), i) + noise_f;  % only allowing additive noise for now
        else
            xOpt(:,i+1) = f(xOpt(:,i), uOpt(i), i) + noise_f;  % only allowing additive noise for now
        end
    end
end