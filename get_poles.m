function polematrix = get_poles(fun)
%GET_POLES   Poles, orders, relevant coefficients of a rational function.
%   polematrix = GET_POLES(fun) gets the poles of the harmonic mapping fun.
%   Each row of polematrix has the form [pol, n, a, b, c], where pol is the
%   pole of order n, and c = -(a0 + conj(b0)).  a and b are the
%   coefficients of the leading terms of fun.f at pol.  If pol = inf:
%       f(z) = a*z^n + b*z^n + O(z^(n-1)) + conj(O(z^(n-1))),
%   and if pol is finite:
%       f(z) = a/(z-pol)^n + conj(b/(z-pol)^n)
%                       + O(1/(z-pol)^(n-1)) + conj(O(1/(z-pol)^(n-1))).


% Partial fraction decomposition:
[res_r, pol_r, poly_r] = residue(fun.rnum, fun.rden);
[res_s, pol_s, poly_s] = residue(fun.snum, fun.sden);

% Replace empty polynomials by zero polynomials:
if ( isempty(poly_r) )
    poly_r = 0;
end
if ( isempty(poly_s) )
    poly_s = 0;
end

% Point at inf:
if ( length(poly_r) > length(poly_s) )
    % inf is a pole of fun.f
    polematrix = [inf, length(poly_r)-1, poly_r(1), 0, ...
        -(poly_r(end) + poly_s(end))];
elseif ( length(poly_r) < length(poly_s) )
    % inf is a pole of fun.f
    polematrix = [inf, length(poly_s)-1, 0, poly_s(1), ...
        -(poly_r(end) + poly_s(end))];
elseif ( ( length(poly_r) == length(poly_s) ) && ( length(poly_r) > 1 ) )
    % inf is a pole of fun.f
    polematrix = [inf, length(poly_r)-1, poly_r(1), poly_s(1), ...
        -(poly_r(end) + poly_s(end))];
else
    % inf is not a pole of fun.f
    polematrix = [];
end

% Finite poles:
tol_poles = 1e-8;

% First, go through poles of r:
pos = 1;
done_s = [];

while ( pos <= length(pol_r) )
    % Current pole:
    z0 = pol_r(pos);
    
    % Occurances in r and s:
    Ir = find(abs(pol_r - z0) < tol_poles);
    Is = find(abs(pol_s - z0) < tol_poles);
    
    % Order and leading coefficients:
    if ( length(Ir) > length(Is) )
        order = length(Ir);
        a = res_r(Ir(end));
        b = 0;
    elseif ( length(Ir) < length(Is) )
        order = length(Is);
        a = 0;
        b = res_s(Is(end));
    else
        order = length(Ir);
        a = res_r(Ir(end));
        b = res_s(Is(end));
    end
    
    done_s = [done_s; Is];
    
    % Compute c = - (a0 + conj(b0)):
    res_r_remain = res_r;
    pol_r_remain = pol_r;
    res_s_remain = res_s;
    pol_s_remain = pol_s;
    % Remove current pole:
    res_r_remain(Ir) = [];
    pol_r_remain(Ir) = [];
    res_s_remain(Is) = [];
    pol_s_remain(Is) = [];
    % Evaluate remainder at the pole:
    [num_r_remain, den_r_remain] = residue(res_r_remain, pol_r_remain, poly_r);
    [num_s_remain, den_s_remain] = residue(res_s_remain, pol_s_remain, poly_s);
    
    a0 = polyval(num_r_remain, z0) / polyval(den_r_remain, z0);
    b0 = polyval(num_s_remain, z0) / polyval(den_s_remain, z0);
    c = -(a0 + conj(b0));
    
    % Add pole to the list of poles:
    polematrix = [polematrix; z0, order, a, b, c];
    
    pos = pos + length(Ir);
end

% Go through remaining poles of s:
remaining_poles_of_s = 1:length(pol_s);
remaining_poles_of_s(done_s) = [];

pos = 1;
while ( pos <= length(remaining_poles_of_s) )
    % Current pole:
    z0 = pol_s(remaining_poles_of_s(pos));
    
    % Check if pole is multiple:
    Is = find(abs(pol_s - pol_s(remaining_poles_of_s(pos))) < tol_poles);
    
    order = length(Is);
    a = 0;
    b = res_s(Is(end));
    
    % Compute c = - (a0 + conj(b0)):
    res_s_remain = res_s;
    pol_s_remain = pol_s;
    % Remove current pole:
    res_s_remain(Is) = [];
    pol_s_remain(Is) = [];
    % Evaluate remainder at the pole:
    [num_s_remain, den_s_remain] = residue(res_s_remain, pol_s_remain, poly_s);
    
    a0 = polyval(fun.rnum, z0) / polyval(fun.rden, z0);
    b0 = polyval(num_s_remain, z0) / polyval(den_s_remain, z0);
    c = -(a0 + conj(b0));
    
    % Add pole to the list of poles:
    polematrix = [polematrix; z0, order, a, b, c];
    
    pos = pos + length(Is);
end

end