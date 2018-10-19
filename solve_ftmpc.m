function [uOpt] = solve_ftmpc(f, Q, R, x0, T, nx, nu)
%SOLVE_FTMPC solve and unconstrained finite time MPC
% f - state evolution function
% Q,R - cost matricies
% x0 - initial state
% T - time cost
% nx - dimention of the state
% nu - dimention of the input
X = sdpvar(nx, T+1);
U = sdpvar(nu, T);

% obj = l(X(:,T+1), zeros(nu));
% c = [X(:,1) == x0];
% 
% for i = 1:T
%     obj = obj + l(X(:,i),U(:,i));
%     c = [c, X(:,i+1) == f(X(:,i), U(i))];
% end

obj = X(:, T+1)'*Q*X(:, T+1);
for i = 1:T
    obj = obj + X(:,i)'*Q*X(:,i) + U(:,i)'*R*U(:,i);
end

c = [];

for i = 1:T
    c = [c, X(:, i+1) == f(X(:,i), U(i))];
end

c = [c, X(:,1) == x0];
optimize(c, obj);

uOpt = double(U);

end

