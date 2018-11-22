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