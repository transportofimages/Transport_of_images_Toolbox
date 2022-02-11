% Reproducing Table 4 from [Sète & Zur, The transport of images method:
... computing all zeros of harmonic mappings by continuation]
% Attention: this m-file requires the 
% Chebfun-toolbox, 
% which can be downloaded from www.chebfun.org, and
% rootsb,
% which can be downloaded from
% https://de.mathworks.com/matlabcentral/fileexchange/44084-computing-common-zeros-of-two-bivariate-functions
%



% 1. Mao-Petters-Witt function for n = 100, r = 0.94
n = 100;
r = 0.94;

% tiroots
rng(1);
caus_pts = 2^8;
tic; % time measurement
fun = mpwfun(n,r);
[zer_mpw,~,~,~,iter] = tiroots(fun, 30, caus_pts);
time_mpw = toc;

% Printing data as in Table 4
disp('Mao-Petters-Witt function: Transport of images');
disp(['Number of zeros: ', num2str(numel(zer_mpw))]);
disp(['Computation time: ', num2str(time_mpw), ' secs.']);

tic
f = @(z) z.*conj(z.^n - r^n) - conj(z.^(n-1));
f1 = @(x,y)real(f(x+1i*y));
f2 = @(x,y)imag(f(x+1i*y));
[x,y] = rootsb(f1,f2,1.2*[-1, 1, -1, 1]);
rootsbzer = x+1i*y;

time_rootsb = toc;
% Printing data as in Table 3

disp('Mao-Petters-Witt function: Rootsb');
disp(['Number of zeros: ', num2str(numel(rootsbzer))]);
disp(['Computation time: ', num2str(time_rootsb), ' secs.']);

tic
f = @(z) z.*conj(z.^n - r^n) - conj(z.^(n-1));
chebf = chebfun2(f, 1.2*[-1, 1, -1, 1], 'vectorize');
chebzer = roots(chebf);

time_chebfun = toc;

disp('Mao-Petters-Witt function: Chebfun');
disp(['Number of zeros: ', num2str(numel(chebzer))]);
disp(['Computation time: ', num2str(time_chebfun), ' secs.']);


% 2. Rhie's function for n = 25, r = 0.9, epsilon = 0.4
n = 25;
r = 0.9;
epsilon = 0.4;



% tiroots
rng(1);
caus_pts = 2^8;
tic; % time measurement
fun = rhiefun_old(n,r,epsilon);
[zer_rhie,crit,caus,~] = tiroots(fun, 30, caus_pts);
time_rhie = toc;

% Printing data as in Table 3
disp('Rhie''s function: Transport of images');
disp(['Number of zeros: ', num2str(numel(zer_rhie))]);
disp(['Computation time: ', num2str(time_rhie), ' secs.']);

tic
f = @(z) (abs(z).^2 - 1).*conj(z.^n) + (epsilon - abs(z).^2)*r^n;
f1 = @(x,y)real(f(x+1i*y));
f2 = @(x,y)imag(f(x+1i*y));
[x,y] = rootsb(f1,f2,1.2*[-1, 1, -1, 1]);
rootsbzer = x+1i*y;

time_rootsb = toc;
% Printing data as in Table 3

disp('Rhie''s function: Rootsb');
disp(['Number of zeros: ', num2str(numel(rootsbzer))]);
disp(['Computation time: ', num2str(time_rootsb), ' secs.']);

tic
f = @(z) (abs(z).^2 - 1).*conj(z.^n) + (epsilon - abs(z).^2)*r^n;
chebf = chebfun2(f, 1.2*[-1, 1, -1, 1], 'vectorize');
chebzer = roots(chebf);
time_chebfun = toc;

disp('Rhie''s function: Chebfun');
disp(['Number of zeros: ', num2str(numel(chebzer))]);
disp(['Computation time: ', num2str(time_chebfun), ' secs.']);

