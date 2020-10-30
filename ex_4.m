% Reproducing Table 2 from [Sète & Zur, The transport of images method:
... computing all zeros of harmonic mappings by continuation]
    
caustic_pts = 2^8;
samples = 50;
number_of_rays = 10;
minn = 3;
maxn = 12;

% 1. Mao-Petters-Witt function

% seed for rng
rng(1)
for n = minn:maxn
    % randomized radii
    rcrit = sqrt((n-2)/n) * (2/(n-2))^(1/n);
    r = 0.7 + (rcrit - 0.7)* rand(1,samples);
    
    time_caustics = 0;
    rays = 0;
    newton_iter = 0;
    tic;
    for k = 1:samples
        %% MPW as f(z) = z - conj(r(z)), n poles, radius r(k)
        fun = mpwfun(n, r(k));
        % Computate critical curves and caustics:
        tcaus = tic;
        [crit, caus] = critical_curves(fun, caustic_pts);
        time_caus = toc(tcaus);
        time_caustics = time_caustics + time_caus;
        
        [zer,~,~,rays_used,iter] = tiroots(fun, number_of_rays, {crit, caus});
        rays = rays + rays_used;
        newton_iter = newton_iter + iter;
    end
    time_transport = toc;
    disp('********************');
    disp(['Mao-Petters-Witt function with n = ',num2str(n), ' with ', num2str(samples), ' runs']);
    disp(['Number of zeros: ', num2str(numel(zer))]);
    disp(['avg. number of harm. Newton iterations: ', num2str(newton_iter/samples)]);
    disp(['avg. number of restarts: ', num2str(rays/samples - 1)]);
    disp(['rel. computation time for the caustics: ', num2str(time_caustics/time_transport)]);
end

% 2. Rhie's function

% see for rng
rng(1)

epsilon(samples) = 0;

for n = minn:maxn
    % randomized radii and residues
    rcrit = sqrt((n-2)/n) * (2/(n-2))^(1/n);
    r = 0.7 + (rcrit - 0.7)* rand(1,samples);
    for jj = 1:samples
        epsstar = rhie_eps(n,r(jj));
        epsilon(jj) = epsstar/2 + rand(1,1)*epsstar/4;
    end
    
    time_caustics = 0;
    rays = 0;
    newton_iter = 0;
    tic;
    for k = 1:samples
        %% Rhie as f(z) = z - conj(r(z)), n poles, radius r(k), epsilon = 0.1
        fun = rhiefun(n,r(k),epsilon(k));
        % Compute critical curves and caustics:
        tcaus = tic;
        [crit, caus] = critical_curves(fun, caustic_pts);
        time_caus = toc(tcaus);
        time_caustics = time_caustics + time_caus;
        
        [zer,~,~,rays_used,iter] = tiroots(fun, number_of_rays, {crit, caus});
        rays = rays + rays_used;
        newton_iter = newton_iter+ iter;
    end
    time_transport = toc;
    disp('********************');
    disp(['Rhie''s function with n = ',num2str(n), ' with ', num2str(samples), ' runs']);
    disp(['Number of zeros: ', num2str(numel(zer))]);
    disp(['avg. number of harm. Newton iterations: ', num2str(newton_iter/samples)]);
    disp(['avg. number of restarts: ', num2str(rays/samples - 1)]);
    disp(['rel. computation time for the caustics: ', num2str(time_caustics/time_transport)]);
end
