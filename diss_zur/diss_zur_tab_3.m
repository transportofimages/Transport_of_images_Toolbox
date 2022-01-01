% This m-file reproduces the results in Table 3 from
%
% [Zur, Jan “On the zeros of harmonic mappings: theory, computation and 
... application". Doctoral Thesis. TU Berlin, 2022].
% 
% Attention: Chebfun and roots are required to run this m-file.
%
% 1. Chebfun 
% https://www.chebfun.org
%
% 2. rootsb
% https://www.mathworks.com/matlabcentral/fileexchange/
... 44084-computing-common-zeros-of-two-bivariate-functions
%

rng(1);

minn = 11;
maxn = 20;
r = 0.75;
epsilon = 0.25;

for n = minn : 1 : maxn
    disp(' ');
    disp(['n=', num2str(n)]);
    disp(' ');
    
    % tiroots
    tic;
    fun = rhiefun(n,r,epsilon);
    tizer = tiroots(fun);
    tires = max(abs(fun.f(tizer)));
    time_ti = toc;
    disp(' ');
    disp('tiroots');
    disp(['Number of zeros: ', num2str(numel(tizer))]);
    disp(['Computation time: ', num2str(time_ti), ' secs.']);
    disp(['Maximum residual: ', num2str(tires)]);

    % Chebfun2 - marching squares
    tic;
    f = @(z) (abs(z).^2 - 1).*conj(z.^n) + (epsilon - abs(z).^2)*r^n;
    f1 = @(x,y)real(f(x+1i*y));
    f2 = @(x,y)imag(f(x+1i*y));
    chebf1 = chebfun2(f1,1.2*[-1, 1, -1, 1],'vectorize');
    chebf2 = chebfun2(f2,1.2*[-1, 1, -1, 1],'vectorize');
    rr = roots(chebf1,chebf2,'marchingsquares');
    cheb_ms_zer = rr(:,1)+1i*rr(:,2);
    cheb_ms_res = max(abs(f(cheb_ms_zer)));
    time_chebfun_ms = toc;
    disp(' ');
    disp('Chebfun2 marching squares');
    disp(['Number of zeros: ', num2str(numel(cheb_ms_zer))]);
    disp(['Computation time: ', num2str(time_chebfun_ms), ' secs.']);
    disp(['Maximum residual: ', num2str(cheb_ms_res)]);

    % Chebfun2 - resultant
    tic;
    f = @(z) (abs(z).^2 - 1).*conj(z.^n) + (epsilon - abs(z).^2)*r^n;
    f1 = @(x,y)real(f(x+1i*y));
    f2 = @(x,y)imag(f(x+1i*y));
    chebf1 = chebfun2(f1,1.2*[-1, 1, -1, 1],'vectorize');
    chebf2 = chebfun2(f2,1.2*[-1, 1, -1, 1],'vectorize');
    rr = roots(chebf1,chebf2,'resultant');
    cheb_res_zer = rr(:,1)+1i*rr(:,2);
    cheb_res_res = max(abs(f(cheb_res_zer)));
    time_chebfun_res = toc;
    disp(' ');
    disp('Chebfun2 resultant');
    disp(['Number of zeros: ', num2str(numel(cheb_res_zer))]);
    disp(['Computation time: ', num2str(time_chebfun_res), ' secs.']);
    disp(['Maximum residual: ', num2str(cheb_res_res)]);

    % rootsb
    tic;
    f = @(z) (abs(z).^2 - 1).*conj(z.^n) + (epsilon - abs(z).^2)*r^n;
    f1 = @(x,y)real(f(x+1i*y));
    f2 = @(x,y)imag(f(x+1i*y));
    [x,y] = rootsb(f1,f2,1.2*[-1, 1, -1, 1]);
    rootsb_zer = x+1i*y;
    rootsb_res = max(abs(f(rootsb_zer)));
    time_rootsb = toc;
    disp(' ');
    disp('rootsb');
    disp(['Number of zeros: ', num2str(numel(rootsb_zer))]);
    disp(['Computation time: ', num2str(time_rootsb), ' secs.']);
    disp(['Maximum residual: ', num2str(rootsb_res)]);    
end