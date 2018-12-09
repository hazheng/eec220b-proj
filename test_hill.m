%% Open Loop test to Verify Dynamics
pmax = 5000; %max power generated by the motor
f_max = 250;
cd = .25; 

v_init = 20;

t = linspace(1,20,201);
dt = .1;

x = zeros(51, 3);
x(1,:) = [v_init, 0, 3600];

u = 1;
for i = 1:201
    x(i+1,:) = hill_model(x(i,:), u, dt);
end

figure
plot(x(:,2), x(:,1))

%% Generate optimize using YALMIP 
n = 200;
dt = .1;
x0 = [10, 0, 0];

dyn = @(x,u) hill_model(x, u, dt);

X = sdpvar(n+1, 3);
U = sdpvar(n, 1);

obj = X(end,3)^2;
% ops = sdpsettings('solver','fmincon');

c = [X(1,:) == x0, X(end,2) >= 200, 0<= U <=1];

for i = 1:n 
    c = [c, X(i+1,:) == dyn(X(i,:), U(i))];
end

optimize(c, obj);

U = double(U);
X = double(X);
