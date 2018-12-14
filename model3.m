%% Lap-based controller
%% System parameters
m = 295;        % Mass (kg)
g = 9.8;        % Gravitational constant (m/s^2)
rho = 1.225;    % Density of air (kg/m^3)
Cd = 0.21;      % Drag coefficient (unitless)
Af = 0.98;      % Frontal area (m^2)
Crr = 0.003;    % Coefficient of rolling resistance (unitless)
Q = 4800*3600;  % Battery capacity (J)
vmax = 30;      % Maximum speed (m/s) (~65mph)
vmin = 0;       % Minimum speed (m/s)
amax = 1;       % Maximum accel (m/s^2)
amin = -5;      % Maximum decel (m/s^2) (~0.5g)
seglengths  = [1000 500 1000 2000]; % Track segment lengths
segslopes   = [   5   0    5   -5]; % Track segment slopes
segsplims   = [  30  30   30   30]; % Track segment speed limits
tracklength = sum(seglengths);
segs = numel(seglengths);
if (abs(sum(seglengths.*sind(segslopes))) > 1e-6)  ...
    || any(segsplims > vmax)                      ...
    || any(segsplims < vmin)
    error('Invalid track');
end
posslopes   = [];
possplims    = [];
for i = 1:segs
    posslopes = [posslopes segslopes(j)*ones(1,seglengths(i))];
    possplims = [possplims segsplims(j)*ones(1,seglengths(i))];
end

%% Vehicle dynamics
Faero = @(v) 0.5*Cd*Af*rho*v^2;     % Aerodynamic drag
Froll = @(v) Crr*m*g*v;             % Rolling resistance
Fgrav = @(theta) m*g*sind(theta);   % Force from gravity
Ftot = @(v,theta) Faero(abs(v)) + Froll(abs(v)) + Fgrav(theta);

%% Solar dynamics
sunin = @() max(0,min(800,600*(1+0.5*randn)));  % Solar input as a function of position?
sunin = @() 0;

%% Model
v = sdpvar(1,segs);
q = sdpvar(1,segs+1);
options = sdpsettings('verbose',0);

stepsize = 100000;
iter = Q/stepsize/10;

constr = [0 <= q <= Q];
J = sum(seglengths./v);
for j = 1:segs
    constr = [constr,vmin <= v(j) <= segsplims(j),    ...
        q(j+1) == q(j) - Ftot(v(j),segslopes(j))*seglengths(j) + sunin()*seglengths(j)/v(j)];
end

for i = 1:iter
    qlap = i*stepsize;
    disp(['Solving iteration ' num2str(i) ' of ' num2str(iter) '...']);
    optimize([constr, q(1) == qlap],J,options);
    vi = value(v)
    qi = value(q)
    Ji = value(J)
    vs = [];
    qs = [];
    for j = 1:segs
        vs = [vs vi(j)*ones(1,seglengths(j))];
        qs = [qs linspace(qi(j),qi(j+1),seglengths(j))];
    end
    subplot(2,1,1);
    plot(1:tracklength,vs); hold on;
    xlabel('Position (m)'); ylabel('Velocity (m/s)');
    subplot(2,1,2);
    plot(1:tracklength,qs); hold on;
    xlabel('Position(m)'); ylabel('Charge (J)');
end

