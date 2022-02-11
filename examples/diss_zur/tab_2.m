% This m-file reproduces the results in Table 2 from
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

rng(1);

minn = 11;
maxn = 20;
r = 0.75;


for n = minn :1: maxn
    disp(' ');
    disp(['n=', num2str(n)]);
    disp(' ');
    
    % tiroots
    tic;
    fun = mpwfun(n,r);
    tizer = tiroots(fun);
    tires = max(abs(fun.f(tizer)));
    time_ti = toc;
    disp(' ');
    disp('tiroots');
    disp(['Number of zeros: ', num2str(numel(tizer))]);
    disp(['Maximum residual: ', num2str(tires)]);
    disp(['Computation time: ', num2str(time_ti), ' secs.']);

    % Chebfun2 - marching squares
    tic;
    f = @(z) z.*conj(z.^n - r^n) - conj(z.^(n-1));
    f1 = @(x,y)real(f(x+1i*y));
    f2 = @(x,y)imag(f(x+1i*y));
    chebf1 = chebfun2(f1,1.2*[-1, 1, -1, 1],'vectorize');
    chebf2 = chebfun2(f2,1.2*[-1, 1, -1, 1],'vectorize');
    rr = roots(chebf1,chebf2,'marchingsquares');
    cheb_ms_zer = rr(:,1)+1i*rr(:,2);
    cheb_ms_res = max(abs(f(cheb_ms_zer)));
    cheb_ms_res_orig = max(abs(fun.f(cheb_ms_zer)));
    time_chebfun_ms = toc;
    disp(' ');
    disp('Chebfun2 marching squares');
    disp(['Number of zeros: ', num2str(numel(cheb_ms_zer))]);
    disp(['Maximum residual (auxiliary function): ', num2str(cheb_ms_res)]);
    disp(['Maximum residual (original function): ', num2str(cheb_ms_res_orig)]);
    disp(['Computation time: ', num2str(time_chebfun_ms), ' secs.']);
    
    % Chebfun2 - resultant
    tic;
    f = @(z) z.*conj(z.^n - r^n) - conj(z.^(n-1));
    f1 = @(x,y)real(f(x+1i*y));
    f2 = @(x,y)imag(f(x+1i*y));
    chebf1 = chebfun2(f1,1.2*[-1, 1, -1, 1],'vectorize');
    chebf2 = chebfun2(f2,1.2*[-1, 1, -1, 1],'vectorize');
    rr = roots(chebf1,chebf2,'resultant');
    cheb_res_zer = rr(:,1)+1i*rr(:,2);
    cheb_res_res = max(abs(f(cheb_res_zer)));
    cheb_res_res_orig = max(abs(fun.f(cheb_res_zer)));
    time_chebfun_res = toc;
    disp(' ');
    disp('Chebfun2 resultant');
    disp(['Number of zeros: ', num2str(numel(cheb_res_zer))]);
    disp(['Maximum residual (auxiliary function): ', num2str(cheb_res_res)]);
    disp(['Maximum residual (original function): ', num2str(cheb_res_res_orig)]);
    disp(['Computation time: ', num2str(time_chebfun_res), ' secs.']);
    
    % rootsb
    tic;
    f = @(z) z.*conj(z.^n - r^n) - conj(z.^(n-1));
    f1 = @(x,y)real(f(x+1i*y));
    f2 = @(x,y)imag(f(x+1i*y));
    [x,y] = rootsb(f1,f2,1.2*[-1, 1, -1, 1]);
    rootsb_zer = x+1i*y;
    rootsb_res = max(abs(f(rootsb_zer)));
    rootsb_res_orig = max(abs(fun.f(rootsb_zer)));
    time_rootsb = toc;
    disp(' ');
    disp('rootsb');
    disp(['Number of zeros: ', num2str(numel(rootsb_zer))]);
    disp(['Maximum residual (auxiliary function): ', num2str(rootsb_res)]);  
    disp(['Maximum residual (original function): ', num2str(rootsb_res_orig)]);
    disp(['Computation time: ', num2str(time_rootsb), ' secs.']);     
end