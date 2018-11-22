% sigma is a 2x1. Assume no covariance for now.
% Procedure:
% At each timestep, make a measurement
% Pass this measurement through the ukf to obtain the state estimate
% pass this state estimate into the stochastic cftoc solver
% stochastic MPC solver constraints are generated according to
% https://arxiv.org/pdf/1709.01201.pdf, using the propagation 
% of the covariances to generate chance constraints
% uses a robust time horizon to avoid the covariance 'blowing up' 

% assume we have a u from MPC, and we now apply it to the system
% then we make a measurement of the true system y
function [feas, xOpt, uOpt, xhat, predErr] = UKF_MPC_linear(A, b, C, P, x0, x0hat, x0cov, M, N,...
                                                      Q, R, xL, xU, uL, uU, cov_f, cov_h)
    
    nx = size(A, 2);
    nu = size(b, 2);
    f = @(x,u) A*x + b*u;
    h = @(x) C*x;
    tr = 2;
    H = [1 0;
        0 1;
        -1 0;
        0 -1];
    g = [xU;
        -xL];
    p = [0.3;
        0.3;
        0.3;
        0.3];
    
    predErr = zeros(2, M-N+1);
    pred_trajs = zeros(2, N, M);
    feas = zeros(1, M);
    xOpt = zeros(2, M+1);
    
    % UKF observer
    xhat = zeros(2, M+1);
    xhat(:,1) = x0hat;
    xCovs = zeros(2, 2,M+1);
    xCovs(:,:,1) = x0cov;

    
    xOpt(:,1) = x0;
    uOpt = zeros(1, M);
    xf = 0;
    
    % M is tfinal
    for i=1:M
        disp(i);
%         tx = xhat(1,i) - (xOpt(1,i) - xf)/(M - i); % for now using exact state
%         disp(tx);
%         % disp(tx);
%         Af = [1 0; 0 0]; 
%         % Af = [];
%         bf = [tx; 0];
        [f_prob, xo, uo, jo] = solve_ucftoc(A, b, P, xhat(:,i), xCovs(:,:,i), N,...
                                       Q, R, H, g, p, f, cov_f, tr, uL, uU,...
                                       xf);
                                   
        if f_prob == 0
            feas(i:end) = 0;
            return
        end
        
        feas(i) = 1;
        pred_trajs(:,:,i) = xo(:,2:end);
        uOpt(i) = uo(1);
        if i == 2
            disp("hello");
        end

        xOpt(:,i+1) = A * xOpt(:,i) + b * uOpt(i) + cov_f * randn(2,1);
        [x_pred, cov_x_pred, sigma_pts_prop, wm0, wc0, ws] = propagate_mean_cov(xhat(:,i), xCovs(:,:,i), f, uOpt(i), nx, cov_f);
        [y_pred, cov_y_pred, cov_xy_pred] = generate_output_prediction(sigma_pts_prop, x_pred, wm0, wc0, ws, h, cov_h);
        y_meas = h(xOpt(:,i)) + cov_h * randn();
        [x_update, cov_update] = ukf_update(x_pred, y_pred, cov_x_pred, cov_y_pred, cov_xy_pred, y_meas);
       
        xhat(:,i+1) = x_update;
        xCovs(:,:,i+1) = cov_update;
        
        if i > N && i <= M - N
            predErr(1,i-N) = norm(pred_trajs(1,:,i-N) - xOpt(1,i-N+1:i), 2);
            predErr(2,i-N) = norm(pred_trajs(2,:,i-N) - xOpt(2,i-N+1:i), 2);
        end
    end    
end
