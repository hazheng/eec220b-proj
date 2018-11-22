function [x_new, cov_new] = ukf_update(x_pred, y_pred, cov_x, cov_y, cov_xy, y)
    K = cov_xy/cov_y;
    x_new = x_pred + K * (y - y_pred);
    cov_new = cov_x - K * cov_y * K';
end
