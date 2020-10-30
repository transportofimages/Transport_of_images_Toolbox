% Reproducing Figure 8 from [Sète & Zur, The transport of images method:
... computing all zeros of harmonic mappings by continuation]
    
% Plot specifications
LW = 2;
fs = 30;
ms = 6;
mss = 9;

% Wilmshurst's polynomial, n = 3
n = 3;
max_zer = n^2;
fun = wilmshurstpoly(n);

% Computing critical curves and caustics
caus_pts = 2^8;
[crit, caus] = critical_curves(fun, caus_pts, "poly");

% Applying the of images method
theta = pi/50;
etatmp = 6*exp(1i*theta);
sol = approx_zeros_at_pole(inf, n, 1, 0, -fun.f(0) + etatmp);
zer = harmonicNewton(@(z) fun.f(z) - etatmp,fun.dh, fun.dg, sol);
[tpath, tpath_pts, crossings] = transport_path(fun, crit, caus, [etatmp,0]);
[~, sol_hist, prediction_hist, tpath_pts, numsol, total_steps, numiter_newton, ~] = transport_phase_hist(fun, zer, tpath, tpath_pts);

% Create Figure 8
subplot(3,4,[1,2]);

plotcpcell(caus,'k', 'LineWidth', 2);
hold on
axis equal
axis([0 7 -1.4 1.4]);

plotcp(tpath_pts, 'ko', 'MarkerFaceColor', 'r', 'Markersize', ms);
text(real(tpath_pts(1))-0.1,imag(tpath_pts(1))-0.27,'$\eta_1$','interpreter','latex');
text(real(tpath_pts(2))-0.1,imag(tpath_pts(2))-0.27,'$\eta_2$','interpreter','latex');
text(real(tpath_pts(3))-0.17,imag(tpath_pts(3))-0.27,'$\eta_3$','interpreter','latex');

subplot(3,4,3);
plotcpcell(caus,'k');
hold on
axis equal
axis([0 .12 -.06 .06]);

plotcp(tpath_pts, 'ko', 'MarkerFaceColor', 'r', 'Markersize', ms);

text(real(tpath_pts(4)),imag(tpath_pts(4))+0.01,'$\eta_4$','interpreter','latex');
text(real(tpath_pts(5)),imag(tpath_pts(5))-0.01,'$\eta_5$','interpreter','latex');
text(real(tpath_pts(6))-0.005,imag(tpath_pts(6))-0.01,'$\eta_6$','interpreter','latex');
text(real(tpath_pts(7))-0.01,imag(tpath_pts(7))+0.012,'$\eta_7$','interpreter','latex');
text(real(tpath_pts(8))-0.005,imag(tpath_pts(8))-0.01,'$\eta_8$','interpreter','latex');
text(real(tpath_pts(9)),imag(tpath_pts(9))-0.01,'$\eta_9$','interpreter','latex');

% Domain for Newton fractals
dom = [-1.5, 2.5, -2, 2];
[Xg, Yg] = ndgrid(dom(1):0.005:dom(2), dom(3):0.005:dom(4));
Zg = Xg + 1i*Yg;

% Create Newton fractals
zer_old = [];
for kk = 1:length(tpath_pts)
    subplot(3,4,3 + kk);
    
    % shifted functions for Newton fractals
    tmpfun = @(z)fun.f(z) - tpath_pts(kk);
    
    soltmp = zeros(max_zer,1);
    soltmp(1:length(sol_hist{kk})) = sol_hist{kk};
    newtonfractal(tmpfun, fun.dh, fun.dg, Zg, soltmp, [], 20);
    
    title("$f - \eta_" + num2str(kk) + "$",'interpreter','latex')
    hold on
    
    % Critical curves
    plotcpcell(crit,'k', 'LineWidth', 2);
    
    if kk > 1
        plotcp(prediction_hist{kk-1}, 'ks', 'MarkerFaceColor', 'w', 'Markersize', mss);
    else
        plotcp(sol, 'ks', 'MarkerFaceColor', 'w', 'Markersize', mss);
    end
    plotcp(sol_hist{kk}, 'ko', 'MarkerFaceColor', 'k', 'Markersize', ms);
    axis(dom);
    
    axis on
    box on
    hold off
end


