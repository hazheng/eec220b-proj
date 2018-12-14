%% Inner Controller Test Harness

%% Car parameters
car.m       = 295;          % Mass (kg)
car.Cd      = 0.21;         % Drag coefficient (unitless)
car.Af      = 0.98;         % Frontal area (m^2)
car.Crr     = 0.003;        % Coefficient of rolling resistance (unitless)
car.Q       = 4800*3600;    % Battery capacity (J)
car.vmax    = 30;           % Maximum speed (m/s) (~65mph)
car.vmin    = 1;            % Minimum speed (m/s)
car.amax    = 1;            % Maximum accel (m/s^2)
car.amin    = -5;           % Maximum decel (m/s^2) (~0.5g)

%% Track parameters
tracksegs   = [1000 500 1000 2000;   ... % Segment lengths (m)
                  5   0    5   -5;   ... % Segment slopes (deg)
                 30  30   30   30];      % Segment speed limits (e.g. curves)
track.segs  = size(tracksegs,2);
track.length = sum(tracksegs(1,:));

%% Environment parameters
g   = 9.8;      % Gravitational acceleration (m/s^2)
rho = 1.225;    % Density of air (kg/m^3)

%% Arithmetic test parameters
q0start = 0.004*car.Q;
q0step  = 0.006*car.Q;
q0end   = 0.100*car.Q;
iter    = floor((q0end-q0start)/q0step)+1;
lgd     = [];

for i = 1:iter
    q0 = q0start + (i-1)*q0step;
    disp(['Solving iteration ' num2str(i) ' of ' num2str(iter) '...']);
    [v,q,J] = lapcontroller(car,g,rho,tracksegs,track,q0);
    plotoverdist(tracksegs(1,:),track,v,q,car.Q);
    lgd = [lgd {['q0 = ' num2str(100*q0/car.Q) '% SoC']}];
end
subplot(2,1,1); title('Arithmetic test progression'); legend(lgd);

%% Geometric test parameters
q0      = 0.001*car.Q;
q0fact  = 1.5;
q0end   = 0.100*car.Q;
lgd     = [];

while q0 < q0end
    disp(['Solving ' num2str(100*q0/q0end) '% of maximum starting SoC...']);
    [v,q,J] = lapcontroller(car,g,rho,tracksegs,track,q0);
    plotoverdist(tracksegs(1,:),track,v,q,car.Q);
    q0 = q0fact*q0;
    lgd = [lgd {['q0 = ' num2str(100*q0/car.Q) '% SoC']}];
end
subplot(2,1,1); title('Geometric test progression'); legend(lgd);
