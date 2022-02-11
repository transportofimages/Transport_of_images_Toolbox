% Reproducing Figure 11 and Table 3 from [Sète & Zur, The transport of
...images method: computing all zeros of harmonic mappings by continuation]


% 1. Mao-Petters-Witt function for n = 100, r = 0.94
n = 100;
r = 0.94;
fun = mpwfun(n,r);

% tiroots
rng(1);
caus_pts = 2^8;
tic; % time measurement
[zer_mpw,~,~,~,iter] = tiroots(fun, 30, caus_pts);
time_mpw = toc;
maxres_mpw = max(abs(fun.f(zer_mpw)));

% Printing data as in Table 3
disp('1. Mao-Petters-Witt function');
disp(['Number of zeros: ', num2str(numel(zer_mpw))]);
disp(['Maximal residual: ', num2str(maxres_mpw)]);
disp(['Computation time: ', num2str(time_mpw), ' secs.']);


% 2. Rhie's function for n = 25, r = 0.9, epsilon = 0.4
n = 25;
r = 0.9;
epsilon = 0.4;
fun = rhiefun_old(n,r,epsilon);
poles = roots(fun.sden);

% tiroots
rng(1);
caus_pts = 2^8;
tic; % time measurement
[zer_rhie,crit,caus,~] = tiroots(fun, 30, caus_pts);
time_rhie = toc;
maxres_rhie = max(abs(fun.f(zer_rhie)));

% Printing data as in Table 3
disp('2. Rhie''s function');
disp(['Number of zeros: ', num2str(numel(zer_rhie))]);
disp(['Maximal residual: ', num2str(maxres_rhie)]);
disp(['Computation time: ', num2str(time_rhie), ' secs.']);


% Create Figure 11

% Plot specifications
ms = 8;

figure(1);
plotcpcell(caus,'k');
hold on
axis equal
hold off

figure(2);
plot_phase(fun.f,1.3);
hold on
plotcpcell(crit,'k');
plotcp(zer_rhie, 'ko', 'MarkerFaceColor', 'k','Markersize',ms);
plotcp(poles, 'ko', 'MarkerFaceColor', 'w','Markersize',ms);
xticks([-1 0 1])
yticks([-1 0 1])
hold off