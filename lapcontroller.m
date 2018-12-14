%% Lap-based controller
function [vi,qi,Ji] = lapcontroller(car,g,rho,tracksegs,track,q0)

%% System parameters
m       = car.m;    % Mass (kg)
Cd      = car.Cd;   % Drag coefficient (unitless)
Af      = car.Af;   % Frontal area (m^2)
Crr     = car.Crr;  % Coefficient of rolling resistance (unitless)
Q       = car.Q;    % Battery capacity (J)
vmax    = car.vmax; % Maximum speed (m/s) (~65mph)
vmin    = car.vmin; % Minimum speed (m/s)
amax    = car.amax; % Maximum accel (m/s^2)
amin    = car.amin; % Maximum decel (m/s^2) (~0.5g)
seglengths  = tracksegs(1,:);  % Track segment lengths
segslopes   = tracksegs(2,:);  % Track segment slopes
segsplims   = tracksegs(3,:);  % Track segment speed limits
tracklength = track.length;
segs        = track.segs;

%% Vehicle dynamics
Faero = @(v) 0.5*Cd*Af*rho*v^2;     % Aerodynamic drag
Froll = @(v) Crr*m*g*v;             % Rolling resistance
Fgrav = @(theta) m*g*sind(theta);   % Force from gravity
Ftot = @(v,theta) Faero(abs(v)) + Froll(abs(v)) + Fgrav(theta);

%% Model
v = sdpvar(1,segs);
q = sdpvar(1,segs+1);
options = sdpsettings('verbose',0);

constr = [0 <= q <= Q];
J = sum(tracksegs(1,:)./v);
for j = 1:segs
    constr = [constr,vmin <= v(j) <= min(vmax,segsplims(j)),    ...
        q(j+1) == q(j) - Ftot(v(j),segslopes(j))*seglengths(j)];
end

optimize([constr, q(1) == q0],J,options);
vi = value(v);
qi = value(q);
Ji = value(J);

