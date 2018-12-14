%% Take 4
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
% track = [1000 1000 1000 1000;
%             5    0   -5    0];
% segs = size(track,2);
% tracklength = sum(track(1,:));
% seglength = tracklength/segs;

%% System dynamics
Faero = @(v) 0.5*Cd*Af*rho*v.^2;     % Aerodynamic drag
Froll = @(v) Crr*m*g*v;             % Rolling resistance
Fgrav = @(theta) m*g*sind(theta);   % Force from gravity
Ftot = @(v,theta) Faero(abs(v)) + Froll(abs(v)) + Fgrav(theta);
sunin = @(t) 400+400*(1+sin(t/60)); % Solar input

%% Simulation parameters
q0 = 25;               % Starting SoC (%)
q0 = Q*min(q0/100,1);   % Max charge (J)
Ts = 60;                % Timestep (s)
Tpred = 15*60;          % Prediction horizon (s)
Trace = 60*60;          % Total race time (s)
N = Tpred/Ts;
M = Trace/Ts;

%% Model
x = sdpvar(1,N+1);
v = sdpvar(1,N);
q = sdpvar(1,N+1);
xol = NaN(M-N,N+1);
vol = NaN(M-N,N);
qol = NaN(M-N,N+1);
xact = NaN(1,M+1);
vact = NaN(1,M);
qact = NaN(1,M+1);
Jact = NaN(1,M-N+1);
xact(1) = x0;
qact(1) = q0;
qtarg = @(k) q0*(1-(k-1)/M);
options = sdpsettings('verbose',0);
baseconstr = [0 <= q,   ...
        0 <= v <= vmax];
for j = 1:N
    baseconstr = [baseconstr,x(j+1) == x(j) + v(j)*Ts];
end
for i = 1:M-N+1
    constr = [baseconstr,   ...
        x(1) == xact(i),    ...
        q(1) == qact(i)];
    J = norm(q-qtarg(i:i+N));
    for j = 1:N
        constr = [constr,q(j+1) == q(j) - Ftot(v(j),0)*v(j)*Ts + sunin((i+j)*Ts)*Ts];%track(2,rem(x(j),tracklength)/seglength+1))*v(j)*Ts];
    end
    disp(['Solving iteration ' num2str(i) ' of ' num2str(M-N+1) '...']);
    optimize(constr,J,options);
    xi = value(x);
    vi = value(v);
    qi = value(q);
    xol(i,:) = xi;
    vol(i,:) = vi;
    qol(i,:) = qi;
    xact(i+1:i+N) = xi(2:end);
    vact(i:i+N-1) = vi;
    qact(i+1:i+N) = qi(2:end);
    Jact(i) = value(J);
end

%% Plot
subplot(3,1,1);
plot(0:M,xact); hold on;
xlabel('Time'); ylabel('Position (m)');
subplot(3,1,2);
plot(0:M-1,vact); hold on;
xlabel('Time'); ylabel('Velocity (m/s)');
subplot(3,1,3);
plot(0:M,qact,'g','linewidth',1.5); hold on; plot(0:M,qtarg(1:M+1),'r--');
xlabel('Time'); ylabel('Charge (J)');
yyaxis right;
plot(0:M-1,sunin((0:M-1)*Ts)); plot(0:M-1,Ftot(vact,0).*vact); plot(0:M-1,diff(qact)/Ts);
ylabel('Power (W)'); legend('Charge','Target charge','Solar in','Power out','Net power');