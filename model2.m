%% Take 2
%% System parameter
m = 295;        % Mass (kg)
g = 9.8;        % Gravitational constant (m/s^2)
rho = 1.225;    % Density of air (kg/m^3)
Cd = 0.21;      % Drag coefficient (unitless)
Af = 0.98;      % Frontal area (m^2)
Crr = 0.003;    % Coefficient of rolling resistance (unitless)
Q = 4800*3600;  % Battery capacity (J)
vmax = 30;      % Speed limit (m/s) (~65mph)
vmin = -10;     % Speed limit in reverse (m/s)
x0 = 0;         % Start position (m)
q0 = Q;         % Start charge (J)

%% Simulation parameters
Ts = 3600;      % Timestep (s)
M = 24;         % 24hr * 60min/hr * 1step/60min
N = 8;          % 8hr * 60min/hr * 1step/60min

%% Vehicle dynamics
Faero = @(v) 0.5*Cd*Af*rho*v^2;     % Aerodynamic drag
Froll = @(v) Crr*m*g*v;             % Rolling resistance
Fgrav = @(theta) m*g*sin(theta);   % Force from gravity
Ftot = @(v,theta) Faero(abs(v)) + Froll(abs(v));% + Fgrav(theta);

%% System dynamics
theta = @(x) 0;%sin(x/100);        % Sinusoidal slope
sunin = @(t) 0;%600*(1+sin(t/60)); % Solar input

%% Model 
x = sdpvar(1,N+1);
q = sdpvar(1,N+1);
v = sdpvar(1,N);
constr = [0 <= q <= Q,  ... % Charge between 0 and max
    0 <= v <= vmax];     % Speed limits, can be turned into soft constraint
for j = 1:N
    constr = [constr, x(j+1) == x(j) + v(j)*Ts];    % Longitudinal motion
end

%% MPC
xOpt = NaN(1,M+1);
qOpt = NaN(1,M+1);
vOpt = NaN(1,M);
xOpt(1) = x0;
qOpt(1) = q0;
for i = 1:M-1
    x = sdpvar(1,M+1-i);
    q = sdpvar(1,M+1-i);
    v = sdpvar(1,M-i);
    constr = [0 <= q <= Q,  ... % Charge between 0 and max
        0 <= v <= vmax,     ... % Speed limits, can be turned into soft constraint
        x(1) == xOpt(i),    ... % Initial position
        q(1) == qOpt(i)];   ... % Initial charge
    for j = 1:M-i
       constr = [constr,                ...
           x(j+1) == x(j) + v(j)*Ts,    ...
           q(j+1) == q(j) - Ftot(v(j),theta(x(j)))*v(j)*Ts];
%            q(j+1) == q(j) - Ftot(v(j),theta(x(j)))*abs(v(j))*Ts];% + sunin(i+j-2)*Ts];   % Battery charge
    end
    disp(['Solving iteration ' num2str(i) '...'])
    J = -x(end);
    options = sdpsettings('verbose',0);
    debug = optimize(constr,J,options);
    xOpti = value(x);
    qOpti = value(q);
    vOpti = value(v);
    xOpt(i+1) = xOpti(2);
    qOpt(i+1) = qOpti(2);
    vOpt(i) = vOpti(1);
end
xOpt(i:end) = xOpti;
qOpt(i:end) = qOpti;
vOpt(i:end) = vOpti;
