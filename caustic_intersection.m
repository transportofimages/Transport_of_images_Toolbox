function crossings = caustic_intersection(fun, crit, caus, segment)
%CAUSTIC_INTERSECTION   Intersection of caustic and segment.
%   crossings = caustic_intersection(fun, crit, caus, [a, b]) lists the
%   intersections of the caustics with the segment [a, b].  Each row of
%   the matrix crossings corresponds to one intersection and has the form:
%       [z0, fun.f(z0), numsolchange, c, first_complex_dilatation]
%   where z0 is a point on the critical set, fun.f(z0) the corresponding
%   caustic point, and c the direction from [Sète & Zur, A Newton method
%   for harmonic mappings in the plane, Theorem 5.2]. numsolchange = +- 2
%   indicates how the number of solutions of fun.f(z) = w changes at the
%   caustic crossing when w travels from a to b.

% Endpoints of the segment:
a = segment(1);
b = segment(2);

%% Mapping to the real line
% T(z) = exp(-1i*angle(b-a)) * (z-a) maps [a, b] to [0, abs(b-a)].
T = @(z) exp(-1i*angle(b-a)) * (z-a);

% Map the caustics:
mapped_caus = caus;
for jj = 1:length(caus)
    mapped_caus{jj} = exp(-1i*angle(b-a)) * (caus{jj} - a);
end

%% Intersections:

crossings = [];

for jj = 1:length(mapped_caus)
    % Caustic intersects real axis if sign(imag(caus)) changes.
    signImCaus = sign(imag(mapped_caus{jj}));
    signImCaus(end) = [];     % Closed curve: last = first entry
    
    % Find intersections:
    shiftdown = circshift(signImCaus, 1);
    shiftup = circshift(signImCaus, -1);
    
    % Exact intersection:
    I0 = find(signImCaus == 0);
    % Make sure the sign changes:
    candidates = shiftup .* shiftdown;
    I0(candidates(I0) > 0) = [];    % Discard, as there is no crossing.
    if ( sum(candidates(I0) == 0) > 0 )
        warning("caustic_intersection:two_consec_intersects", ...
            "Two consecutive caustic points on the segment.")
    end
    
    % Keep only intersections of the mapped caustic with [0, abs(b-a)]:
    I0(real(mapped_caus{jj}(I0)) < 0) = [];
    I0(real(mapped_caus{jj}(I0)) > abs(b-a)) = [];
    
    % Change in number of solutions:
    numsolchange = -2 * shiftup(I0);
    
    crossings = [crossings; crit{jj}(I0), caus{jj}(I0), numsolchange];
    
    
    % Sign change: Intersection found, but not the exact point.
    Ic = find(signImCaus .* shiftup < 0);
    
    if ( ~isempty(Ic) )
        critshiftup = circshift(crit{jj}, -1); % can be improved
        
        [omega_num, omega_den] = dilatation2(fun);
        [domega_num, domega_den] = polyder(omega_num, omega_den);
        omega = @(z) polyval(omega_num, z) ./ polyval(omega_den, z);
        domega = @(z) polyval(domega_num, z) ./ polyval(domega_den, z);
        
        for kk = 1:length(Ic)
            z1 = crit{jj}(Ic(kk));
            z2 = critshiftup(Ic(kk));
            % Compute point closer to the intersection with biscetion:
            z0 = critcaus_local_refinement(z1, z2, @(z) T(fun.f(z)), omega, domega);
            fz0 = fun.f(z0);
            
            if ( ( real(T(fz0)) >= 0 ) && real(T(fz0)) <= abs(b-a) )
                % Keep this crossing since T(fz0) in [0, abs(b-a)].
                numsolchange = -2*shiftup(Ic(kk));
                crossings = [crossings; z0, fz0, numsolchange];
            end
        end
    end
end

%% Crossing direction c and (first) complex dilatation:

z0 = crossings(:,1);
complex_dilat = conj(fun.dg(z0))./fun.dh(z0);
c = -(fun.ddh(z0).*complex_dilat + conj(fun.ddg(z0).*complex_dilat))/2;
crossings = [crossings, c, complex_dilat];

%% Sort intersections from a to b:

[~, perm] = sort(abs(crossings(:, 2) - a));
crossings = crossings(perm, :);

end




function z0 = critcaus_local_refinement(z1, z2, f, omega, domega, tol)
%CRITCAUS_LOCAL_REFINEMENT   Critical and caustic points of intersection.
%   z0 = critcaus_local_refinement(z1, z2, f, omega, domega, tol)
%   critical point z0 such that f(z0) is on [a,b].

% f = @(z) T(fun.f(z));

if ( ( nargin < 6 ) || isempty(tol) )
    tol = 1e-12;
end

% omega(z_j) = e^{i t_j}:
t1 = angle(omega(z1));
t2 = angle(omega(z2));
if ( t2 < t1 )
    % From the parametrisation, t_1 < t_2 holds.  Otherwise, we have
    % crossed the disontinuity of the argument funcion.
    % Make it continuous:
    t2 = t2 + 2*pi;
end

f1 = f(z1);

% Bisection:
for jj = 1:10 % maxtstep
    t0 = (t1 + t2)/2;
    z0 = newton(@(z) omega(z) - exp(1i*t0), domega, z1);
    fz0 = f(z0);
    
    if ( abs(imag(fz0)) <= tol )
        % Crossing found.
        return;
    elseif ( imag(f1) * imag(fz0) < 0 )
        % Crossing in [t1, t0]:
        t2 = t0;
    else
        % Crossing in [t0, t2]:
        t1 = t0;
        z1 = z0;
        f1 = fz0;
    end
end
end