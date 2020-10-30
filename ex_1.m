% Reproduction of Figure 7 and Table 1 from [Sète & Zur, The transport of 
...images method: computing all zeros of harmonic mappings by continuation]

% 1. Wilmshurst's polynomial

% construction of the function
n = 3;
fun = wilmshurstpoly(n);

% Applying tiroots
caus_pts = 2^8;
theta = pi/50;
[zer_wilmshurst, ~, caus_wilmshurst, ~, iter_wilmshurst, steps_wilmshurst, steps_failed_wilmshurst] = tiroots_ray(fun, theta, caus_pts, "poly");

% maximal residual
maxres_wilmshurst = max(abs(fun.f(zer_wilmshurst)));


% 2. The Mao-Petters-Witt function

% construction of the function
n = 7;
r = 0.7;
fun = mpwfun(n,r);

% Applying tiroots
caus_pts = 2^8;
theta = pi/50;
[zer_mpw, ~, caus_mpw, ~, iter_mpw, steps_mpw, steps_failed_mpw] = tiroots_ray(fun, theta, caus_pts);

% maximal residual
maxres_mpw = max(abs(fun.f(zer_mpw)));


% 3. Rhie's function
n = 7;
r = 0.7;
epsilon = 0.1;
fun = rhiefun(n,r,epsilon);

% Applying tiroots
caus_pts = 2^8;
theta = pi/50;
[zer_rhie, ~, caus_rhie, ~, iter_rhie, steps_rhie, steps_failed_rhie] = tiroots_ray(fun, theta, caus_pts);

% maximal residual
maxres_rhie = max(abs(fun.f(zer_rhie)));

% 4. f(z) = z^2 + 1/conj(z) + 1/conj(z+1) + 2log|z|
fun = harmonicRat([1 0 0], 1, [2 1], [1 1 0], 2, 0);

% Applying tiroots
caus_pts = 2^8;
theta = pi/50;
[zer_log,~,~,~,iter_log,steps_log,steps_failed_log] = tiroots_ray(fun,theta,caus_pts);

% maximal residual
maxres_log = max(abs(fun.f(zer_log)));

disp('1. f(z) = z^2 + 1/conj(z) + 1/conj(z+1) + 2log|z|')
disp(['Number of zeros: ', num2str(numel(zer_log))]);
disp(['Maximal residual: ', num2str(maxres_log)]);
disp(['harm. Newton iterations: ', num2str(iter_log)]);
disp(['Number of transport steps: ', num2str(steps_log - steps_failed_log)]);
disp(['Refinements of the transport path: ', num2str(steps_failed_log)]);


% Printing data as in Table 1
disp('2. Wilmshurts''s function');
disp(['Number of zeros: ', num2str(numel(zer_wilmshurst))]);
disp(['Maximal residual: ', num2str(maxres_wilmshurst)]);
disp(['Harm. Newton iterations: ', num2str(iter_wilmshurst)]);
disp(['Number of transport steps: ', num2str(steps_wilmshurst - steps_failed_wilmshurst)]);
disp(['Refinements of the transport path: ', num2str(steps_failed_wilmshurst)]);

disp('3. Mao-Petters-Witt function');
disp(['Number of zeros: ', num2str(numel(zer_mpw))]);
disp(['Maximal residual: ', num2str(maxres_mpw)]);
disp(['Harm. Newton iterations: ', num2str(iter_mpw)]);
disp(['Number of transport steps: ', num2str(steps_mpw - steps_failed_mpw)]);
disp(['Refinements of the transport path: ', num2str(steps_failed_mpw)]);

disp('4. Rhies''s function')
disp(['Number of zeros: ', num2str(numel(zer_rhie))]);
disp(['Maximal residual: ', num2str(maxres_rhie)]);
disp(['harm. Newton iterations: ', num2str(iter_rhie)]);
disp(['Number of transport steps: ', num2str(steps_rhie - steps_failed_rhie)]);
disp(['Refinements of the transport path: ', num2str(steps_failed_rhie)]);



% Plot the caustics as in Figure 7
subplot(1,3,1);
plotcpcell(caus_wilmshurst,'k');
axis square

subplot(1,3,2);
plotcpcell(caus_mpw,'k');
axis square

subplot(1,3,3);
plotcpcell(caus_rhie,'k');
axis square

