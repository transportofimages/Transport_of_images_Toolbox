function zer = approx_zeros_at_pole(pol, n, a, b, c)
%APPROX_ZEROS_AT_POLE   Approximate zeros close to a pole.
%   zer = APPROX_ZEROS_AT_POLE(inf, n, a, b, c) returns approximations
%   to the zeros of f at infinity, where
%       f(z) = a*z^n + conj(b*z^n) + O(z^(n-1)) + conj(O(z^(n-1))),
%   and c = - (a0 + conj(b0)).
%
%   zer = APPROX_ZEROS_AT_POLE(zj, n, a, b, c) is the same but at a finite
%   pole zj, i.e.,
%       f(z) = a/(z-zj)^n + conj(b/(z-zj)^n) + O(1/(z-zj)^(n-1)) ...
%                                            + conj(O(1/(z-zj)^(n-1))).


% Solution of a*y + conj(b*y) = c, see [Sète & Zur, A Newton method for
% harmonic mappings in the plane, Lemma 4.2].
w = (conj(a)*c - conj(b*c)) / ((abs(a) + abs(b))*(abs(a) - abs(b)));

if ( isinf(pol) )
    % Solve z^n = w
    p = eye(1, n + 1);
    p(end) = - w;
    zer = roots(p);
    
elseif ( isnan(pol) )
    error("approx_zeros_at_pole:PolIsNaN", "Pole = NaN not allowed.")
    
else
    % Finite pole, solve (z-pol)^{-n} = w, i.e., (z-pol)^n = 1/w
    p = eye(1, n + 1);
    p(end) = -1/w;
    zer = roots(p);
    zer = zer + pol;
end

end