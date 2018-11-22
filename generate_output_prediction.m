% generates output prediction and covariances based on ukf
function [y_pred, cov_y_pred, cov_xy] = generate_output_prediction(sigma_pts_prop, x_pred, wm0, wc0, ws, h, cov_h)
    L = size(sigma_pts_prop, 2);
    L = (L - 1)/2;
    sigma_pts_output = zeros(1, size(sigma_pts_prop,2));
    for i = 1:2*L+1
        sigma_pts_output(:,i) = h(sigma_pts_prop(:,i));
    end
    y_pred = sigma_pts_output(:,1) * wm0;
    
    for i = 2:2*L+1
        y_pred = y_pred + ws(i-1) * sigma_pts_output(:,i);
    end
    
    cov_y_pred = cov_h + wc0 * (sigma_pts_output(:,1)-y_pred)*...
        (sigma_pts_output(:,1)-y_pred)';
    cov_xy = wc0 * (sigma_pts_prop(:,1) - x_pred) * (sigma_pts_output(:,1) - y_pred)';
    for i = 2:2*L+1
        cov_y_pred = cov_y_pred + ws(i-1) * (sigma_pts_output(:,i) - y_pred)*...
            (sigma_pts_output(:,i) - y_pred)';
        cov_xy = cov_xy + ws(i-1) * (sigma_pts_prop(:,i) - x_pred)*...
            (sigma_pts_output(:,i) - y_pred)'; 
    end
end
