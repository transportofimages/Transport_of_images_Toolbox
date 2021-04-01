% Reproducing Figure 10 from [Sète & Zur, The transport of images method:
... computing all zeros of harmonic mappings by continuation]
% Attention: this m-file requires the 
% Chebfun-toolbox, 
% which can be downloaded from www.chebfun.org, and
% rootsb,
% which can be downloaded from
% https://de.mathworks.com/matlabcentral/fileexchange/44084-computing-common-zeros-of-two-bivariate-functions
%
% Note, that running this file (with sample size 50) may take up to 
% 25 minutes.

tStart = tic;

rng(1)

% Turn Chebfun2 warning off:
warning('off', 'CHEBFUN:CHEBFUN2:constructor:vectorize')

caustic_pts = 2^8;
samples = 50;
number_of_rays = 10;
minn = 3;
maxn = 12;

time_transport = zeros(maxn,1);
transport_wrong = zeros(maxn,1);

time_chebfun = zeros(maxn,1);
chebfun_wrong = zeros(maxn,1);

time_rootsb = zeros(maxn,1);
rootsb_wrong = zeros(maxn,1);

%% 1. Mao-Petters-Witt function

disp('*************************');
disp('Mao-Petters-Witt function');

for n = minn:maxn
    % randomized radii
    rcrit = sqrt((n-2)/n) * (2/(n-2))^(1/n);
    r = 0.7 + (rcrit - 0.7)* rand(1,samples);
    
    tic
    for k = 1:samples
        %% MPW as f(z) = z - conj(r(z)), n, r(k),
        fun = mpwfun(n, r(k));
        [zer,~,~,rays_used,iter] = tiroots(fun, number_of_rays, caustic_pts);
        if size(zer) ~= 3*n + 1
            transport_wrong(n) = transport_wrong(n) + 1;
        end
    end
    time_transport(n) = toc;
    
    %% Chebfun2
    R = 1.2; % assumes the zeros to be in 1.2*[-1,1,-1,1]
    tic
    for k = 1:samples
        f = @(z) z.*conj(z.^n - r(k)^n) - conj(z.^(n-1));
        chebf = chebfun2(f, R*[-1, 1, -1, 1], 'vectorize');
        chebzer = roots(chebf);
        if size(chebzer) ~= 3*n + 1
            chebfun_wrong(n) = chebfun_wrong(n) + 1;
        end
    end
    time_chebfun(n) = toc;
    
     %% Rootsb
    R = 1.2; % assumes the zeros to be in 1.2*[-1,1,-1,1]
    tic
    for k = 1:samples
        f = @(z) z.*conj(z.^n - r(k)^n) - conj(z.^(n-1));
        f1 = @(x,y)real(f(x+1i*y));
        f2 = @(x,y)imag(f(x+1i*y));
        [x,y] = rootsb(f1,f2,1.2*[-1, 1, -1, 1]);
        rootsbzer = x+1i*y;
        if size(rootsbzer) ~= 3*n + 1
            rootsb_wrong(n) = rootsb_wrong(n) + 1;
        end
    end
    time_rootsb(n) = toc;
end

%% Plotting
set(0,'defaultTextInterpreter','latex');
figure(1)
plot(minn:maxn,time_transport(minn:maxn)/samples,'-dk','MarkerFaceColor','k','MarkerSize',10,'Linewidth',3);
hold on
plot(minn:maxn,time_chebfun(minn:maxn)/samples,':ok','MarkerFaceColor','k','MarkerSize',10,'Linewidth',3);
plot(minn:maxn,time_rootsb(minn:maxn)/samples,'--vk','MarkerFaceColor','k','MarkerSize',10,'Linewidth',3);
legend('Transport of images', 'Chebfun2','Rootsb','Location','NORTHWEST')
xlabel('$n$')
ylabel('avg. time in secs.')
set(gca,'FontSize',18)
hold off
title('Mao-Petters-Witt function')

%% Wrong results
for k = minn:maxn
    if transport_wrong(k) > 0
        disp(['Transport of images for n = ', num2str(k), ': ', num2str(transport_wrong(k)), '/', num2str(samples), ' wrong']);
    end
    if chebfun_wrong(k) > 0
        disp(['Chebfun2 for n = ', num2str(k), ': ', num2str(chebfun_wrong(k)), '/', num2str(samples), ' wrong']);
    end
    if rootsb_wrong(k) > 0
        disp(['rootsb for n = ', num2str(k), ': ', num2str(rootsb_wrong(k)), '/', num2str(samples), ' wrong']);
    end
end

%% End Mao-Petters-Witt function

rng(1)

time_transport = zeros(maxn,1);
transport_wrong = zeros(maxn,1);

time_chebfun = zeros(maxn,1);
chebfun_wrong = zeros(maxn,1);

%% 2. Rhie's function

disp('*************************');
disp('Rhie''s function');

epsilon(samples) = 0;
for n = minn:maxn
    % randomized radii
    rcrit = sqrt((n-2)/n) * (2/(n-2))^(1/n);
    r = 0.7 + (rcrit - 0.7)* rand(1,samples);
    for jj = 1: samples
        epsstar = rhie_eps(n,r(jj));
        epsilon(jj) = epsstar/2 + rand(1,1)*epsstar/4;
    end
    
    tic
    for k = 1:samples
        %% Rhie as f(z) = z - conj(r(z)), with n, r(k), epsilon(k)
        fun = rhiefun(n, r(k), epsilon(k));
        [zer,~,~,rays_used,iter] = tiroots(fun, number_of_rays, caustic_pts);
        
        if size(zer) ~= 5*n
            transport_wrong(n) = transport_wrong(n) + 1;
        end
    end
    time_transport(n) = toc;
    
    %% Chebfun2
    R = 1.2; % assumes the zeros to be in 1.2*[-1,1,-1,1]
    tic
    for k = 1:samples
        f = @(z) (abs(z).^2 - 1).*conj(z.^n) + (epsilon(k) - abs(z).^2)*r(k)^n;
        chebf = chebfun2(f, R*[-1, 1, -1, 1], 'vectorize');
        chebzer = roots(chebf);
        if size(chebzer) ~= 5*n
            chebfun_wrong(n) = chebfun_wrong(n) + 1;
        end
    end
    time_chebfun(n) = toc;
    
    %% Rootsb
    tic
    for k = 1:samples
        f = @(z) (abs(z).^2 - 1).*conj(z.^n) + (epsilon(k) - abs(z).^2)*r(k)^n;
        f1 = @(x,y)real(f(x+1i*y));
        f2 = @(x,y)imag(f(x+1i*y));
        [x,y] = rootsb(f1,f2,1.2*[-1, 1, -1, 1]);
        rootsbzer = x+1i*y;
        if size(rootsbzer) ~= 5*n
            rootsb_wrong(n) = rootsb_wrong(n) + 1;
        end
    end
    time_rootsb(n) = toc;
end

% Turn Chebfun2 warning off:
warning('on', 'CHEBFUN:CHEBFUN2:constructor:vectorize')

%% Plotting
set(0,'defaultTextInterpreter','latex');
figure(2)
plot(minn:maxn,time_transport(minn:maxn)/samples,'-dk','MarkerFaceColor','k','MarkerSize',10,'Linewidth',3);
hold on
plot(minn:maxn,time_chebfun(minn:maxn)/samples,':ok','MarkerFaceColor','k','MarkerSize',10,'Linewidth',3);
plot(minn:maxn,time_rootsb(minn:maxn)/samples,'--vk','MarkerFaceColor','k','MarkerSize',10,'Linewidth',3);
legend('Transport of images', 'Chebfun2','Rootsb','Location','NORTHWEST')
xlabel('$n$')
ylabel('avg. time in secs.')
set(gca,'FontSize',18)
hold off
title('Rhie''s function')

%% Wrong results
disp('*************************');
disp('Rhie''s function');
for k = minn:maxn
    if transport_wrong(k) > 0
        disp(['Transport of images for n = ', num2str(k), ': ', num2str(transport_wrong(k)), '/', num2str(samples), ' wrong']);
    end
    if chebfun_wrong(k) > 0
        disp(['Chebfun2 for n = ', num2str(k), ': ', num2str(chebfun_wrong(k)), '/', num2str(samples), ' wrong']);
    end
    if rootsb_wrong(k) > 0
        disp(['rootsb for n = ', num2str(k), ': ', num2str(rootsb_wrong(k)), '/', num2str(samples), ' wrong']);
    end
end

%% End Rhie's function

% Turn Chebfun2 warning on:
warning('on', 'CHEBFUN:CHEBFUN2:constructor:vectorize')
toc(tStart)