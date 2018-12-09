%% Test out energy heuristics / cost models for the high level energy scheduler
% Mock strategy over the course of a day
% sim params 
N_day = 32;
N = N_day*3
dt = 8/N_day; %time unit is in hours, 
soc_max = 5000; %units in W/h
pout_max = 50000; % units in W
% solar = zeros(N,1);
solar = 1400*sin(linspace(pi/6, 5*pi/6, N_day));
solar = repmat(solar,1,3);
%account for evening and morning charge
energy_evening = 800*4;
solar(N_day+1) = solar(N_day+1) + energy_evening/dt;
solar(2*N_day+1) = solar(2*N_day+1)+ energy_evening/dt;

soc_init = .95*soc_max;
% compute the mean power draw
% state just consists of SOC 
soc = sdpvar(N+1, 1);
u = sdpvar(N,1); 
soc(1) = soc_init;

energy_total = soc_init + sum(solar)*dt;
power_avg = energy_total/24;

% assign(u, power_avg/pout_max*ones(N, 1));
c = [soc_max >= soc >= 0, 1 >= u >= 0];

for i = 1:N
    c = [c, soc(i+1) == soc(i) - dt*pout_max*u(i) + solar(i)*dt];
end

obj = -sum(u.^(1/3));

ops = sdpsettings('fmincon.maxiter', 10000);
optimize(c, obj, ops); 

u = double(u);
soc = double(soc);

figure
subplot(2,1,1)
plot(soc)
subplot(2,1,2)
plot(u*pout_max)