% This m-file reproduces the results in Table 4 from
%
% [Zur, Jan “On the Zeros of Harmonic Mappings: Analysis, Computation and 
... Application". PhD thesis. TU Berlin, 2022].
% 
% Attention: Chebfun and rootsb are required to run this m-file.
%
% 1. Chebfun 
% https://www.chebfun.org
%
% 2. rootsb
% https://www.mathworks.com/matlabcentral/fileexchange/
... 44084-computing-common-zeros-of-two-bivariate-functions
%



% 1. Mao-Petters-Witt function for n = 100, r = 0.94
n = 100;
r = 0.94;

disp('Mao-Petters-Witt function, n = 100, rho = 0.94');
disp(' ');

% tiroots
rng(1);
tic; % time measurement
fun = mpwfun(n,r);
[zer_mpw,~,~,~,iter] = tiroots(fun);
time_mpw = toc;


% Printing data as in Table 4
disp('tiroots');
disp(['Number of zeros: ', num2str(numel(zer_mpw))]);
disp(['Maximum residual: ', num2str(max(abs(fun.f(zer_mpw))))]);
disp(['Computation time: ', num2str(time_mpw), ' secs.']);
disp(' ');

% Chebfun2
tic
f = @(z) z.*conj(z.^n - r^n) - conj(z.^(n-1));
f1 = @(x,y)real(f(x+1i*y));
f2 = @(x,y)imag(f(x+1i*y));
chebf1 = chebfun2(f1,1.2*[-1, 1, -1, 1],'vectorize');
chebf2 = chebfun2(f2,1.2*[-1, 1, -1, 1],'vectorize');
rr = roots(chebf1,chebf2);
cheb_zer = rr(:,1)+1i*rr(:,2);
time_chebfun = toc;

disp('Chebfun');
disp(['Number of zeros: ', num2str(numel(cheb_zer))]);
disp(['Maximum residual: ', num2str(max(abs(fun.f(cheb_zer))))]);
disp(['Computation time: ', num2str(time_chebfun), ' secs.']);
disp(' ');

% rootsb
tic
f = @(z) z.*conj(z.^n - r^n) - conj(z.^(n-1));
f1 = @(x,y)real(f(x+1i*y));
f2 = @(x,y)imag(f(x+1i*y));
[x,y] = rootsb(f1,f2,1.2*[-1, 1, -1, 1]);
rootsbzer = x+1i*y;

time_rootsb = toc;

% Printing data as in Table 4
disp('rootsb');
disp(['Number of zeros: ', num2str(numel(rootsbzer))]);
disp(['Maximum residual: ', num2str(max(abs(fun.f(rootsbzer))))]);
disp(['Computation time: ', num2str(time_rootsb), ' secs.']);
disp(' ');


% 2. Rhie's function for n = 25, r = 0.8, epsilon = 0.4
n = 25;
r = 0.8;
epsilon = 0.4;


disp('Rhie''s function, n = 25, rho = 0.8, eps = 0.4');
disp(' ');


% tiroots
rng(1);
tic;
fun = rhiefun(n,r,epsilon);
[zer_rhie,crit,caus,~] = tiroots(fun);
time_rhie = toc;

% Printing data as in Table 4
disp('tiroots');
disp(['Number of zeros: ', num2str(numel(zer_rhie))]);
disp(['Maximum residual: ', num2str(max(abs(fun.f(zer_rhie))))]);
disp(['Computation time: ', num2str(time_rhie), ' secs.']);
disp(' ');

% Chebfun2
tic
f = @(z) (abs(z).^2 - 1).*conj(z.^n) + (epsilon - abs(z).^2)*r^n;
f1 = @(x,y)real(f(x+1i*y));
f2 = @(x,y)imag(f(x+1i*y));
chebf1 = chebfun2(f1,1.2*[-1, 1, -1, 1],'vectorize');
chebf2 = chebfun2(f2,1.2*[-1, 1, -1, 1],'vectorize');
rr = roots(chebf1,chebf2);
cheb_zer = rr(:,1)+1i*rr(:,2);
time_chebfun = toc;

% Printing data as in Table 4
disp('Rhie''s function: Chebfun');
disp(['Number of zeros: ', num2str(numel(cheb_zer))]);
disp(['Maximum residual: ', num2str(max(abs(fun.f(cheb_zer))))]);
disp(['Computation time: ', num2str(time_chebfun), ' secs.']);
disp(' ');

% rootsb
tic
f = @(z) (abs(z).^2 - 1).*conj(z.^n) + (epsilon - abs(z).^2)*r^n;
f1 = @(x,y)real(f(x+1i*y));
f2 = @(x,y)imag(f(x+1i*y));
[x,y] = rootsb(f1,f2,1.2*[-1, 1, -1, 1]);
rootsbzer = x+1i*y;

time_rootsb = toc;

% Printing data as in Table 4
disp('Rhie''s function: rootsb');
disp(['Number of zeros: ', num2str(numel(rootsbzer))]);
disp(['Maximum residual: ', num2str(max(abs(fun.f(rootsbzer))))]);
disp(['Computation time: ', num2str(time_rootsb), ' secs.']);
disp(' ');