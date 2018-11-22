% make a conservative choice of terminal set: basically, the terminal 
% set for each iteration is the set that you need to reach that you could
% get to the final point going the minimum possible speed; ie. at time 20s,
% the terminal set is the set at which you could reach the terminal point
% going at (distance_to_final_point)/(tf - 20) m/s, assuming no
% disturbances (not persistently feasible, but provides a reasonable
% terminal set). 
function [feas, xOpt, uOpt, xhat, predErr] = MPC_stochastic_naive(A, b, C, L, P, x0, x0hat, M, N,...
                                                      Q, R, xL, xU, uL, uU, f_cov, h_cov)
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