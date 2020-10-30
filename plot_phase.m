function plot_phase(f, box, res)
%PLOT_PHASE   Phase plot of a complex function.
%   PLOT_PHASE(F, BOX) plots the phase of the function handle F in the box
%   defined by BOX (optional, by default BOX = 1.5); see below.
%
%   PLOT_PHASE(F, BOX, RES) specifies the resolution (RES = 256 by default).
%   If length(res) = 1, then xres = yres = res, if length(res) = 2, then
%   [xres, yres] = res.
%
%
% The input BOX is interpreted based on its length (let bj := box(j)):
%   length(box) == 0  square of width 3 around 0
%   length(box) == 1  square around 0 with width abs(b1)
%   length(box) == 2  square around (real(b1), imag(b1)), width abs(b2)
%   length(box) == 3  square around (b1,b2) with width abs(b3)
%   length(box) == 4  interpret as values [xmin, xmax, ymin, ymax]

% Creating adequate bounds for plot-domain.
if ( ( nargin < 2 ) || isempty(box) )
    box = 1.5;
end

box = interpret_plot_box(box);
xmin = box(1);
xmax = box(2);
ymin = box(3);
ymax = box(4);

if ( ( nargin == 3 ) && ~isempty(res) )
    if ( length(res) == 1 )
        xres = res;
        yres = res;
    elseif ( length(res) == 2 )
        xres = res(1);
        yres = res(2);
    else
        warning("plot_phase:resSizeWrong", 'Wrong size of res, switching to default resolution.')
        xres = 256; yres = 256;
    end
else
    xres = 256; yres = 256;
    % xres = 1024; yres = 1024;
end

% Creating meshgrid
x = linspace(xmin, xmax, xres);
y = linspace(ymin, ymax, yres);
[x, y] = meshgrid(x, y);
z = x + 1i*y;

% Evaluating function
fz = f(z);

% Phaseplot
% p = surf(real(z), imag(z), 0*fz, angle(-fz));
p = surf(real(z), imag(z), zeros(size(fz)), angle(-fz));
set(p, 'EdgeColor', 'none');
caxis([-pi, pi]);
colormap hsv(256)

view(0,90)
axis equal
% axis off

hold on
contour(real(z), imag(z), angle(-fz), 0, 'LineColor', 'w')
% This has no actual impact on the plot itself, but preserves the square
% appearance of the figure.
hold off

end

function box_out = interpret_plot_box(box_in)
% function box_out = interpret_plot_box(box_in)
%
% Generate a box for plot usage in plot_phase and related functions.
%
% The input box is interpreted based on its length (let bj := box_in(j) ):
%   length(box_in) == 0  square of width 3 around 0
%   length(box_in) == 1  square around 0 with width abs(b1)
%   length(box_in) == 2  square around (real(b1), imag(b1)), width abs(b2)
%   length(box_in) == 3  square around (b1,b2) with width abs(b3)
%   length(box_in) == 4  interpret as values [xmin, xmax, ymin, ymax]

if length(box_in) == 1
    box_out = [-abs(box_in),abs(box_in),-abs(box_in),abs(box_in)];
elseif length(box_in) == 2
    midx = real(box_in(1));
    midy = imag(box_in(1));
    width = abs(box_in(2));
    box_out = [midx - width, midx + width, midy - width, midy + width];
elseif length(box_in) == 3
    midx = box_in(1);
    midy = box_in(2);
    width = abs(box_in(3));
    box_out = [midx - width, midx + width, midy - width, midy + width];
elseif length(box_in) == 4
    box_out = box_in;
else
    warning("plot_phase:BoundFormat", 'Bad format for bounds, switching to default values.')
    box_out = [-1.5, 1.5, -1.5, 1.5];
end

if box_out(2) < box_out(1)
    warning("plot_phase:x_order", 'xmin > xmax. Switching xmin and xmax.')
    box_out([1 2]) = box_out([1 2]);
end

if box_out(4) < box_out(3)
    warning("plot_phase:y_order", 'ymin > ymax. Switching ymin and ymax.')
    box_out([3 4]) = box_out([4 3]);
end

end