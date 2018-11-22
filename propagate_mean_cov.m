% f is dynamics function
% does one predict step
% implementation of https://arxiv.org/pdf/1709.01201.pdf
function [mean_pred, cov_pred, sigma_pts_prop, wm0, wc0, ws] = propagate_mean_cov(x, cov_x, f, u, nx, cov_f)
    alpha = 1e-3;
    beta = 2;
    kappa = 0;
    
    L = nx;
    lambda = alpha^2 * (L + kappa) - L;
    wm0 = lambda/(L + lambda);
    wc0 = lambda/(L + lambda) + (1 - alpha^2 + beta);
    
    % Doing this with sdpvars might get hairy...
    if all(class(cov_x) == 'double') && all(cov_x(:) == 0)
        sqrt_cov = zeros(2);
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
        sigma_pts = [sigma_pts, x - sqrt_cov(i-1-L,:)'];
    end
    ws = repmat(1/(2 * (L + lambda)), 1, 2 * L);
    % sigma_pts_prop = zeros(nx, 2*L+1);
    % sigma_pts_prop(:,1) = f(sigma_pts(:,1), u);
    sigma_pts_prop = f(sigma_pts(:,1), u);
    mean_pred = wm0 * sigma_pts_prop(:,1);
    
    for i=2:2*L+1
        % sigma_pts_prop(:,i) = f(sigma_pts(:,i), u);
        sigma_pts_prop = [sigma_pts_prop, f(sigma_pts(:,i),u)];
        mean_pred = mean_pred + ws(i-1) * sigma_pts_prop(:,i);
        
    end
    
    cov_pred = cov_f + wc0 * (sigma_pts_prop(:,1)-mean_pred)*...
        (sigma_pts_prop(:,1)-mean_pred)';
    for i = 2:2*L
        cov_pred = cov_pred + ws(i-1) * (sigma_pts_prop(:,i) - mean_pred)*...
            (sigma_pts_prop(:,i) - mean_pred)';
    end
end