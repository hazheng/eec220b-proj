% assuming linearity for now, can be extended to a nonlinear function later
% P(H * x <= g) >= p
function [feas, xOpt, uOpt, JOpt] = solve_ucftoc(A, b, P, x0, x_cov, N,...
                                                 Q, R, H, g, p, f, cov_f, tr, uL, uU,...
                                                 target_x)
    nu = size(b, 2);
    nx = size(A, 2);
    nc = size(H, 1);  % number of rows of H = number of constraints per timestep
    
    x = sdpvar(nx, N+1);  % this is the mean, or expected state value
    u = sdpvar(nu, N);
    
    % have to do shenanigans with declaring arrays of sdpvars
    % because of https://yalmip.github.io/naninmodel/
    % covs_tr = zeros(nx, nx, tr+1);  % only track covs up to tr
    covs_tr = x_cov;
    
    cost = (x(:,N+1) - target_x)' * P * (x(:,N+1) - target_x);
    for i =[1:N]
      cost = cost + (x(:,i) - target_x)'*Q*(x(:,i)-target_x) + u(:,i)'*R*u(:,i);
    end
    
    
    constraints = [uL <= u(:,1) <= uU];
    index = 2;
    for i = 1:N-1
      constraints = [constraints; uL <= u(:,index) <= uU];
      index = index + 1;
    end
    
    constraints = [constraints; x(:,1) == x0];
    for i=2:N+1
        [mean_pred, cov_pred, sigma_pts_prop, wm0, wc0, ws] = propagate_mean_cov(x(:,i-1), covs_tr(:,:,i-1), f, u(:,i-1), nx, cov_f);
%         sdisplay(mean_pred);
        constraints = [constraints; x(:,i) == mean_pred];
        if i <= tr
            % covs_tr(:,:,i) = cov_pred;
            covs_tr = cat(3, covs_tr, cov_pred);
        else
            % covs_tr(:,:,i) = covs_tr(:,:,i-1);
            covs_tr = cat(3, covs_tr, covs_tr(:,:,i-1));
        end
        
        % reformulation as a probabilistic constraint based on covariance
        % reformulation as per https://arxiv.org/pdf/1709.01201.pdf
        for j=1:nc
            constraints = [constraints;
                (quantile(p(j)) * sqrtm(H(j,:) * covs_tr(:,:,i) * H(j,:)')+...
                H(j,:) * x(:,i)) <= g(j)];
        end            
    end
        
    options = sdpsettings('verbose', 0);
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