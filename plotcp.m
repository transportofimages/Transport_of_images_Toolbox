function varargout = plotcp(pts, varargin)
%PLOTCP   Plot points in the complex plane.
%   PLOTCP(PTS, ...) is short for PLOT(REAL(PTS), IMAG(PTS), ...).

h = plot(real(pts), imag(pts), varargin{:});

if ( nargout == 1 )
    varargout{1} = h;
end

end