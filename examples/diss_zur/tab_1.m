% This m-file reproduces the results in Table 1 from
%
% [Zur, Jan “On the Zeros of Harmonic Mappings: Analysis, Computation and 
... Application". PhD thesis. TU Berlin, 2022].
% 
% Attention: Chebfun and rootsb are required to run this m-file.
%

rng(1);

tic
n = 3;
fun = wilmshurstpoly(n);
zer_wilmshurst = tiroots(fun, [], [], "poly");
time_wilmshurst = toc;
maxres_wilmshurst = max(abs(fun.f(zer_wilmshurst)));

tic
n = 7;
r = 0.7;
fun = mpwfun(n,r);
zer_mpw = tiroots(fun);
time_mpw = toc;
maxres_mpw = max(abs(fun.f(zer_mpw)));


tic
n = 8;
r = 0.7;
epsilon = 0.2;
fun = rhiefun(n,r,epsilon);
zer_rhie = tiroots(fun);
time_rhie = toc;
maxres_rhie = max(abs(fun.f(zer_rhie)));


% Printing data as in Table 1
disp('Wilmshurts''s polynomial');
disp(['Number of zeros: ', num2str(numel(zer_wilmshurst))]);
disp(['Maximum residual: ', num2str(maxres_wilmshurst)]);
disp(['Time (ms): ', num2str(ceil(1000*time_wilmshurst))]);

disp('Mao-Petters-Witt function');
disp(['Number of zeros: ', num2str(numel(zer_mpw))]);
disp(['Maximum residual: ', num2str(maxres_mpw)]);
disp(['Time (ms): ', num2str(ceil(1000*time_mpw))]);

disp('Rhies''s function')
disp(['Number of zeros: ', num2str(numel(zer_rhie))]);
disp(['Maximum residual: ', num2str(maxres_rhie)]);
disp(['Time (ms): ', num2str(ceil(1000*time_rhie))]);


