%% Take 3
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
track = [1000 1000 1000 1000;
            5    0   -5    0];
segs = size(track,2);

%% Simulation parameters
Ts = 3600;      % Timestep (s)
M = 24;         % 24hr * 60min/hr * 1step/60min
N = 8;          % 8hr * 60min/hr * 1step/60min

%% Vehicle dynamics
Faero = @(v) 0.5*Cd*Af*rho*v^2;     % Aerodynamic drag
Froll = @(v) Crr*m*g*v;             % Rolling resistance
Fgrav = @(theta) m*g*sind(theta);   % Force from gravity
Ftot = @(v,theta) Faero(abs(v)) + Froll(abs(v)) + Fgrav(theta);

%% System dynamics
sunin = @(t) 0;%600*(1+sin(t/60)); % Solar input

% Model
v = sdpvar(1,segs);
segt = sdpvar(1,segs);
q = sdpvar(1,segs+1);
stepsize = 100000;
iter = Q/stepsize/10;

options = sdpsettings('verbose',0);
for i = 1:iter
    qlap = i*stepsize;
    
    constr = [0 <= q,   ...
        0 <= segt,      ...
        0 <= v <= vmax  ...
        segt.*v == track(1,:)];
    J = sum(segt);
    for j = 1:segs
        constr = [constr,q(j+1) == q(j) - Ftot(v(j),track(2,j))*track(1,j)];
    end
    
    disp(['Solving iteration ' num2str(i) ' of ' num2str(iter) '...']);
    optimize([constr, q(1) == qlap],J,options);
    vi = value(v)
    segti = value(segt)
    qi = value(q)
    Ji = value(J)
    segti.*vi
    subplot(3,1,1);
    plot(1:segs,value(v)); hold on;
    xlabel('Track segment'); ylabel('Velocity (m/s)');
    subplot(2,1,2);
    plot(1:segs+1,value(q)); hold on;
    xlabel('[Start of] track segment'); ylabel('Charge (J)');
end

