% Reproducing Figure 9 from [Sète & Zur, The transport of images method:
... computing all zeros of harmonic mappings by continuation]
    
% Plot specifications
LW = 1;
ms = 5;

%Rhie's function for n = 7, r = 0.7, epsilon = 0.1
n = 7;
r = 0.7;
epsilon = 0.1;
fun = rhiefun_old(n,r,epsilon);

% applying tiroots
caus_pts = 2^8;
theta = pi/50;
[zer, ~, caus, tsteps_list, ~, ~, steps_failed] = tiroots_ray(fun, theta, caus_pts);
tp = pts_from_tsteps(tsteps_list);

% Create Figure 9
subplot(2,2,[1, 2]);
axis equal
axis([0 .3 -0.03 .03]);
hold on
plotcpcell(caus,'k','LineWidth', LW);
hold on
plotcp(tp, 'ko', 'MarkerFaceColor', 'r', 'Markersize', ms);

xticks([0 0.05 0.1 0.15 0.2 0.25 0.3])
yticks([-0.03 0 0.03])
box on
hold off

subplot(2,2,3);
axis equal
axis([.012 .016 .0005 .0015]);
hold on
plotcpcell(caus,'k','LineWidth', LW);
plotcp(tp, 'ko', 'MarkerFaceColor', 'r', 'Markersize', ms);

xticks([0.012 .016])
yticks([0.0005  0.0015])

box on
hold off

subplot(2,2,4);
axis equal
axis([.0574 .059 0.0035 .0039]);
hold on
plotcpcell(caus,'k','LineWidth', LW);
plotcp(tp, 'ko', 'MarkerFaceColor', 'r', 'Markersize', ms);

xticks([.0574 .059])
yticks([0.0035 .0039])

box on
hold off