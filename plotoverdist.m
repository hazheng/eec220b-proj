function [] = plotoverdist(seglengths,track,v,q,Q)
vs = [];
qs = [];
for i = 1:track.segs
    vs = [vs v(i)*ones(1,seglengths(i))];
    qs = [qs linspace(q(i),q(i+1),seglengths(i))];
end
subplot(2,1,1);
plot(1:track.length,vs); hold on;
xlabel('Position (m)'); ylabel('Velocity (m/s)');
subplot(2,1,2);
plot(1:track.length,100*qs/Q); hold on;
xlabel('Position(m)'); ylabel('% Charge');