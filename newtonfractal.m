function h = newtonfractal(f, dh, dg, Zgrid, zer, tol, maxit, plotzer, plotcrit, tolsep)
%NEWTONFRACTAL   Plot basins of attractions in the harmonic Newton method.
%   NEWTONFRACTAL(f, dh, dg, Zgrid, zer)
%   f(z) = h(z) + conj(g(z)), dh = h', dg = g', zer are the zeros of f,
%   Zgrid is the grid on which to evaluate the function.
%
%   NEWTONFRACTAL(f, dh, dg, Zgrid, zer, tol, maxit, plotzer, plotgrid)
%   tol = tolerance (default: tol = 1e-14)
%   maxit = maximal number of Newton steps (default: maxit = 20)
%   plotzer = 1: plot zeros zer
%   plotcrit = 1: plot critical curve.


% Set default values:
if ( ( nargin <= 5 ) || isempty(tol) )
    tol = 1e-14;
end

if ( ( nargin <= 6 ) || isempty(maxit) )
    maxit = 50;
end

if ( ( nargin <= 7 ) || isempty(plotzer) )
    plotzer = 0;
end

if ( ( nargin <= 8 ) || isempty(plotcrit) )
    plotcrit = 0;
end

if ( ( nargin <= 9 ) || isempty(tolsep) )
    tolsep = 1e-3;
end


% Run harmonic Newton on Zgrid:

[Zk,numiter] = harmonicNewton(f, dh, dg, Zgrid, maxit, tol);

% Determine closest zero:
zer = zer(:);
fzer = NaN(size(Zgrid));
for jj = 1:length(zer)
    fzer(abs(Zk - zer(jj)) <= tolsep) = jj;
end

fzer(isnan(fzer)) = maxit*numel(zer) + 1; % Mark `leftovers' as not converged.

% shading
fzer = fzer + numiter*numel(zer); 

% color for shading
colormap(hsv(numel(zer)));
new_color = zeros((maxit+1)*numel(zer)+1,3);
hsv_map = rgb2hsv(colormap);
for k = 0:maxit
     if k <= maxit/2
        hsv_map(:,2) = 0.2 + (2*k/maxit)*0.8;
    else
        hsv_map(:,3) = 1-(2*k-maxit)/maxit;
    end
%     hsv_map(:,3) = 1-k/maxit;
    new_color(k*numel(zer)+1:(k+1)*numel(zer),:) = hsv2rgb(hsv_map);
end
colormap(new_color)
% Plot
h = surf(real(Zgrid), imag(Zgrid), zeros(size(fzer)), fzer);
set(h, 'EdgeColor', 'none');
   
view(0,90), axis equal, axis off
caxis([1 (maxit+1)*length(zer)+1]); % numer of colors = number of colors + 1
  
    
hold on
if ( plotcrit == 1 )
    % Plot critical set:
    Jacval = abs(dh(Zgrid).^2) - abs(dg(Zgrid)).^2;
    contour(real(Zgrid), imag(Zgrid), Jacval, [0, 0], 'k-', ...
        'LineWidth', 2)
end
if ( plotzer == 1 )
    % Plot zeros:
    plotcp(zer, 'kx', 'LineWidth', 2)
end
hold off

end