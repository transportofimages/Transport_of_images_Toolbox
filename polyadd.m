function p = polyadd(a, b)
%POLYADD   Add two polynomials.
%   p = polyadd(a, b) performs p = a + b.

if ( nargin < 2 )
    error("polyadd:input", "Not enough Input Argmuents.")
end

a=reshape(a, 1, []);    % Make sure inputs are polynomial row vectors.
b=b(:).';               % This makes a row as well.

na = length(a);
nb = length(b);

p = [zeros(1, nb - na) a] + [zeros(1, na - nb) b];  % Pad with zeros.

end