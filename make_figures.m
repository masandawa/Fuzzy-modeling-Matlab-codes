function make_figures(FD)
%MAKE_FIGURES  Draw the seven data figures of the manuscript and save them as PDF.
%   make_figures        loads results/figure_data.mat (written by the main script)
%   make_figures(FD)    uses the structure FD directly
%
%   Figures (manuscript numbering; Fig. 2, the model flow diagram, is drawn in LaTeX):
%     Fig. 1  fig2_membership.pdf        fuzzy intensity groups and intensity-dependent parameters
%     Fig. 3  fig6_validation.pdf        model evaluation
%     Fig. 4  fig_bifurcation.pdf        bifurcation diagram and global convergence
%     Fig. 5  fig3_R0.pdf                reproduction numbers versus intensity
%     Fig. 6  fig4_fuzzy_equilibria.pdf  fuzzy equilibria
%     Fig. 7  fig_coverage.pdf           critical deworming coverage
%     Fig. 8  fig8_scenarios.pdf         fuzzy scenario trajectories 2026-2030

here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
if nargin < 1
    S = load(fullfile(here, 'results', 'figure_data.mat')); FD = S.FD;
end
figdir = fullfile(here, 'figures'); if ~exist(figdir, 'dir'), mkdir(figdir); end

% ---------------- colours and helpers ----------------
BLUE = hex('2a78d6'); ORANGE = hex('eb6834'); AQUA = hex('1baf7a'); YELLOW = hex('eda100');
MAGENTA = hex('e87ba4'); GREEN = hex('008300'); RED = hex('e34948'); INK = hex('0b0b0b'); INK2 = hex('52514e');
SC = {BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN};          % scenario colours B,D1,D2,M,I1,I2
GCOL = {AQUA, YELLOW, RED};                                  % light/low, moderate, heavy/high
HG = [24 1000 3000; 1000 3000 5000; 3000 5000 8000];          % hookworm groups (epg)
MG = [1 2.5 3.5; 2.5 3.5 4.5; 3.5 4.5 5];                    % parasite-density groups (log10)
HC = [0 2000; 2000 4000; 4000 9000]; MC = [10 1e3; 1e3 1e4; 1e4 1e5];   % crisp WHO classes
HN = {'light', 'moderate', 'heavy'}; MN = {'low', 'moderate', 'high'};
set(0, 'DefaultAxesFontSize', 9, 'DefaultLineLineWidth', 1.6, 'DefaultAxesBox', 'off');
q = FD.q;

% =============== Fig. 1: membership and intensity-dependent parameters ===============
f = newfig('Fig. 1: fuzzy intensity groups', [100 100 900 650]);
e = linspace(0, 9000, 1801); r = logspace(0.7, 5, 2601);
subplot(2,2,1); hold on; h = zeros(1,3);
for g = 1:3
    m = trimf(e, HG(g,:));
    fill([e fliplr(e)], [m zeros(size(m))], light(GCOL{g}, 0.35), 'EdgeColor', 'none');
    h(g) = plot(e, m, 'Color', GCOL{g}, 'LineWidth', 1.4);
end
xlabel('Hookworm intensity \epsilon (eggs per gram)'); ylabel('Membership \Gamma(\epsilon)');
title('(a) Linguistic intensity groups (WHO cut-offs)'); legend(h, HN, 'Location', 'northeast', 'Box', 'off'); grid on;
subplot(2,2,2); hold on;
pih = clip01((e - 24)/(4000 - 24)); th = q.theta0 + (q.thetaM - q.theta0)*clip01(e/4000);
h1 = plot(e, pih, 'Color', BLUE); h2 = plot(e, th, 'Color', ORANGE);
plot([0 9000], [1 1], ':', 'Color', INK2, 'LineWidth', 0.8);
xlabel('Hookworm intensity \epsilon (epg)'); ylabel('Value'); title('(b) Worm-burden-dependent parameters');
legend([h1 h2], {'\pi(\epsilon)/\pi_{max} egg contamination', '\theta(\epsilon) malaria susceptibility'}, 'Location', 'east', 'Box', 'off'); grid on;
subplot(2,2,3); hold on; h = zeros(1,3);
for g = 1:3
    m = trimf(log10(r), MG(g,:));
    fill([r fliplr(r)], [m zeros(size(m))], light(GCOL{g}, 0.35), 'EdgeColor', 'none');
    h(g) = plot(r, m, 'Color', GCOL{g}, 'LineWidth', 1.4);
end
set(gca, 'XScale', 'log'); xlim([10 1e5]);
xlabel('Parasite density \rho (parasites/\muL)'); ylabel('Membership \Gamma(\rho)');
title('(c) Linguistic parasite-density groups'); legend(h, MN, 'Location', 'northwest', 'Box', 'off'); grid on;
subplot(2,2,4); hold on;
ch = clip01((r - 10)/(2000 - 10)); gm = 0.0035 + (1/30 - 0.0035)*clip01((r - 10)/(10000 - 10));
h1 = plot(r, ch, 'Color', BLUE); h2 = plot(r, gm/(1/30), 'Color', ORANGE);
set(gca, 'XScale', 'log'); xlim([5 6e4]);
xlabel('Parasite density \rho (parasites/\muL)'); ylabel('Normalised value'); title('(d) Density-dependent parameters');
legend([h1 h2], {'c(\rho)/c_{max} infectiousness', '\gamma_m(\rho)/\gamma_s clearance'}, 'Location', 'east', 'Box', 'off'); grid on;
savefig_(f, fullfile(figdir, 'fig2_membership.pdf'));

% =============== Fig. 3: model evaluation ===============
f = newfig('Fig. 3: model evaluation', [120 120 900 360]);
subplot(1,2,1); hold on;
h1 = plot(2015.5 + FD.tt/365, 100*FD.Pmh, 'Color', BLUE);
h2 = errorbar([2015.5 2021.5], [21.6 11.8], [0.3 0.4], [0.3 0.2], 'o', 'Color', INK, 'MarkerFaceColor', INK, 'LineWidth', 1);
h3 = errorbar(2023.5, 11.3, 0.3, 0.2, 's', 'Color', ORANGE, 'MarkerFaceColor', ORANGE, 'LineWidth', 1);
xlabel('Year'); ylabel('Malaria prevalence, SAC (%)'); title('(a) National school-age malaria');
legend([h1 h2 h3], {'model (fitted to 2015, 2021)', 'survey (fitting)', 'survey 2023 (validation)'}, 'Box', 'off'); grid on;
subplot(1,2,2); hold on;
patch([15.5 36.6 36.6 15.5], [-0.08 -0.08 1.45 1.45], light(ORANGE, 0.2), 'EdgeColor', 'none');
plot([25 25], [-0.08 1.45], 'Color', ORANGE, 'LineWidth', 1.2);
VB = FD.VB;
for i = 1:size(VB, 1)
    hA = plot(100*VB(i,2:3), [VB(i,1) VB(i,1)], 'Color', BLUE, 'LineWidth', 3);
    hB = plot(100*VB(i,4:5), VB(i,1) + [0.04 0.04], 'Color', AQUA, 'LineWidth', 3);
end
ylim([-0.08 1.45]);
xlabel('Hookworm reinfection by week 18 (%)'); ylabel('Membership degree \alpha'); title('(b) Reinfection, Pemba');
legend([hA hB], {'window 105 d', 'window 63 d (prepatency)'}, 'Location', 'northeast', 'Box', 'off'); grid on;
savefig_(f, fullfile(figdir, 'fig6_validation.pdf'));

% =============== Fig. 4: bifurcation diagram and global convergence ===============
f = newfig('Fig. 4: bifurcation and global convergence', [140 140 950 720]);
BD = FD.BD; thr = interp1(BD(:,2), BD(:,1), 1);
subplot(2,2,1); hold on;
patch([thr 1 1 thr], [0 0 30 30], light(ORANGE, 0.35), 'EdgeColor', 'none');
h1 = plot(BD(:,1), 100*BD(:,4), '--', 'Color', INK2, 'LineWidth', 1.5);
h2 = plot(BD(:,1), 100*BD(:,3), 'Color', BLUE, 'LineWidth', 2);
ylim([0 max(100*BD(:,3))*1.05]);
text(1.45, 3, {'hookworm-mediated', 'persistence (R_M<1<R_M^W)'}, 'Color', ORANGE, 'FontSize', 8);
xlabel('R_M (varying m_v)'); ylabel('Malaria prevalence (%)'); title('(a) Bifurcation diagram');
legend([h1 h2], {'without hookworm (E_M)', 'with hookworm (core)'}, 'Location', 'northwest', 'Box', 'off'); grid on;
subplot(2,2,2); hold on;
plot(BD(:,1), 1e3*BD(:,5), 'Color', BLUE); plot([min(BD(:,1)) max(BD(:,1))], [0 0], 'Color', INK2, 'LineWidth', 0.8);
xlabel('R_M'); ylabel('max Re(\lambda) (10^{-3} day^{-1})'); title('(b) Stability of the attractor'); grid on;
GC = FD.GC; GF = FD.GF;
ttl = {'(c) E_0', '(d) E_W', '(e) E_M', '(f) R_M<1<R_M^W', '(g) E^*'};
for rg = 1:5
    subplot(2,5,5 + rg); hold on;
    sub = GC(GC(:,1) == rg, :); col = SC{rg};
    for n = unique(sub(:,2))'
        t = sub(sub(:,2) == n, :);
        plot(100*t(:,5), 100*t(:,4), 'Color', light(col, 0.55), 'LineWidth', 0.7);
        plot(100*t(1,5), 100*t(1,4), '.', 'Color', light(col, 0.9), 'MarkerSize', 6);
    end
    g = GF(GF(:,1) == rg, :);
    plot(100*g(6), 100*g(5), 'o', 'MarkerFaceColor', INK, 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
    title(ttl{rg}, 'FontSize', 8); xlabel('Hookworm (%)', 'FontSize', 7);
    if rg == 1, ylabel('Malaria (%)', 'FontSize', 7.5); end
    set(gca, 'FontSize', 7); grid on;
end
savefig_(f, fullfile(figdir, 'fig_bifurcation.pdf'));

% =============== Fig. 5: reproduction numbers versus intensity ===============
f = newfig('Fig. 5: reproduction numbers', [160 160 950 330]);
FR = FD.FR;
subplot(1,3,1); hold on;
ymax = max(FD.RWg(FD.eg <= 9000))*1.05;
for g = 1:3, patch([HC(g,1) HC(g,2) HC(g,2) HC(g,1)], [0 0 ymax ymax], light(GCOL{g}, 0.15), 'EdgeColor', 'none'); end
plot(FD.eg, FD.RWg, 'Color', BLUE);
for i = 1:3, plot(FR(i,2), FR(i,4), 'o', 'MarkerFaceColor', GCOL{i}, 'MarkerEdgeColor', 'w', 'MarkerSize', 7); end
plot([0 9000], [1 1], ':', 'Color', INK2, 'LineWidth', 0.8);
xlim([0 9000]); ylim([0 ymax]);
xlabel('\epsilon (epg)'); ylabel('Reproduction number'); title('(a) Hookworm, R_W'); grid on;
subplot(1,3,2); hold on;
plot(FD.eg, FD.PhiG, 'Color', ORANGE);
yl = [min(FD.PhiG) max(FD.PhiG)*1.02];
plot([0 9000], [1 1], ':', 'Color', INK2, 'LineWidth', 0.8);
plot([FD.eps_RW1 FD.eps_RW1], yl, '--', 'Color', INK2, 'LineWidth', 0.8);
text(FD.eps_RW1 + 200, 0.98*max(FD.PhiG), sprintf('R_W=1 at\n%d epg', round(FD.eps_RW1)), 'Color', INK2, 'FontSize', 8, 'VerticalAlignment', 'top');
xlim([0 9000]); ylim(yl);
xlabel('\epsilon (epg)'); ylabel('Invasion factor \Phi(\epsilon)'); title('(b) Hookworm effect on malaria'); grid on;
subplot(1,3,3); hold on;
ok = FD.rg > 0; ymax = max(FD.RMg)*1.05;
for g = 1:3, patch([MC(g,1) MC(g,2) MC(g,2) MC(g,1)], [0 0 ymax ymax], light(GCOL{g}, 0.15), 'EdgeColor', 'none'); end
plot(FD.rg(ok), FD.RMg(ok), 'Color', BLUE);
for i = 1:3, plot(FR(3+i,2), FR(3+i,4), 'o', 'MarkerFaceColor', GCOL{i}, 'MarkerEdgeColor', 'w', 'MarkerSize', 7); end
plot([10 1e5], [1 1], ':', 'Color', INK2, 'LineWidth', 0.8);
set(gca, 'XScale', 'log'); xlim([10 1e5]); ylim([0 ymax]);
xlabel('\rho (parasites/\muL)'); ylabel('R_M(\rho)'); title('(c) Malaria'); grid on;
savefig_(f, fullfile(figdir, 'fig3_R0.pdf'));

% =============== Fig. 6: fuzzy equilibria ===============
f = newfig('Fig. 6: fuzzy equilibria', [180 180 900 620]);
Ee = FD.EQe; Er = FD.EQr; FE = FD.FE;
subplot(2,2,1); hold on;
h = [plot(Ee(:,1), 100*Ee(:,3), 'Color', BLUE), plot(Ee(:,1), 100*Ee(:,2), 'Color', ORANGE), plot(Ee(:,1), 100*Ee(:,4), 'Color', AQUA)];
xlabel('Hookworm intensity \epsilon (epg)'); ylabel('Equilibrium prevalence (%)'); title('(a) Equilibria vs worm burden');
legend(h, {'hookworm', 'malaria', 'co-infection'}, 'Location', 'east', 'Box', 'off'); grid on;
subplot(2,2,2); hold on; ok = Er(:,1) > 0;
plot(Er(ok,1), 100*Er(ok,3), 'Color', BLUE); plot(Er(ok,1), 100*Er(ok,2), 'Color', ORANGE); plot(Er(ok,1), 100*Er(ok,4), 'Color', AQUA);
set(gca, 'XScale', 'log'); xlim([10 1e5]);
xlabel('Parasite density \rho (parasites/\muL)'); ylabel('Equilibrium prevalence (%)'); title('(b) Equilibria vs parasite density'); grid on;
subplot(2,2,3); hold on; h = zeros(1,3);
for g = 1:3, h(g) = membplot(FE(FE(:,1) == 1 & FE(:,2) == g, :), 8, 9, GCOL{g}); end
xlabel('Co-infection prevalence (%)'); ylabel('Membership degree \alpha'); title('(c) Fuzzy co-infection, hookworm groups');
legend(h, HN, 'Box', 'off'); grid on;
subplot(2,2,4); hold on; h = zeros(1,3);
for g = 1:3, h(g) = membplot(FE(FE(:,1) == 2 & FE(:,2) == g, :), 4, 5, GCOL{g}); end
xlabel('Malaria prevalence (%)'); ylabel('Membership degree \alpha'); title('(d) Fuzzy malaria, density groups');
legend(h, MN, 'Box', 'off'); grid on;
savefig_(f, fullfile(figdir, 'fig4_fuzzy_equilibria.pdf'));

% =============== Fig. 7: critical deworming coverage ===============
f = newfig('Fig. 7: critical deworming coverage', [200 200 950 360]);
PC = FD.PC; PI = FD.PCI;
subplot(1,2,1); hold on;
for g = 1:3, patch([HC(g,1) HC(g,2) HC(g,2) HC(g,1)], [0 0 105 105], light(GCOL{g}, 0.15), 'EdgeColor', 'none'); end
h1 = plot(PI(:,1), 100*nanabove(PI(:,2)), 'Color', BLUE);
h2 = plot(PI(:,1), 100*nanabove(PI(:,3)), 'Color', ORANGE);
h3 = plot(PC(:,1), 100*nanabove(PC(:,2)), '--', 'Color', BLUE, 'LineWidth', 1);
plot(PC(:,1), 100*nanabove(PC(:,3)), '--', 'Color', ORANGE, 'LineWidth', 1);
plot([0 4500], [75 75], '--', 'Color', RED, 'LineWidth', 0.9); text(2600, 78, 'WHO 75% coverage', 'Color', RED, 'FontSize', 7);
plot([0 4500], [22.5 22.5], ':', 'Color', INK2, 'LineWidth', 0.8); text(2600, 25.5, 'school-only max (SAC 30%)', 'Color', INK2, 'FontSize', 7);
xlim([0 4500]); ylim([0 105]);
xlabel('Hookworm intensity \epsilon (epg)'); ylabel('Critical community coverage \phi_c (%)');
title('(a) Coverage needed to interrupt hookworm');
legend([h1 h2 h3], {'annual (periodic rounds)', 'semi-annual (periodic rounds)', 'averaged approximation'}, 'Location', 'east', 'Box', 'off'); grid on;
subplot(1,2,2); hold on;
al = linspace(0, 1, 401); g = HG(1,:);
lo = g(1) + al*(g(2) - g(1)); hi = g(3) - al*(g(3) - g(2));
cols = {BLUE, ORANGE}; h = zeros(1,2);
for k = 1:2
    fl = interp1(PI(:,1), PI(:,1+k), lo); fh = interp1(PI(:,1), PI(:,1+k), hi);
    x = [100*min(fl, 1.05), fliplr(100*min(fh, 1.05))]; y = [al, fliplr(al)];
    fill(x, y, light(cols{k}, 0.25), 'EdgeColor', 'none'); h(k) = plot(x, y, 'Color', cols{k});
end
plot([75 75], [0 1], '--', 'Color', RED, 'LineWidth', 0.9);
xlim([0 106]); ylim([0 1]);
xlabel('Critical coverage \phi_c (%)  [>100%: unattainable]'); ylabel('Membership degree \alpha');
title('(b) Fuzzy \phi_c, light-intensity group'); legend(h, {'annual', 'semi-annual'}, 'Location', 'northwest', 'Box', 'off'); grid on;
savefig_(f, fullfile(figdir, 'fig_coverage.pdf'));

% =============== Fig. 8: fuzzy scenario trajectories ===============
f = newfig('Fig. 8: control scenarios 2026-2030', [220 220 950 680]);
i0 = find(abs(FD.alv - 0) < 1e-9); i5 = find(abs(FD.alv - 0.5) < 1e-9); i1 = find(abs(FD.alv - 1) < 1e-9);
t = 2026 + FD.tgrid/365;
vars = {1, 100, 'Malaria prevalence (%)', '(a)'; 2, 100, 'Hookworm prevalence (%)', '(b)'; ...
        3, 100, 'Co-infection prevalence (%)', '(c)'; 4, 1, 'Malaria infections /1000 PY', '(d)'};
show = [1 4 6]; LG = {'B: baseline', 'M: vector control', 'I2: vector control + community deworming + WASH'};
for k = 1:4
    subplot(2,2,k); hold on; h = zeros(1,3); v = vars{k,1}; sc = vars{k,2};
    for j = 1:3
        s = show(j); c = SC{s};
        E0 = FD.ENV{i0, s}; E5 = FD.ENV{i5, s}; E1 = FD.ENV{i1, s};
        fill([t fliplr(t)], sc*[E0(v,:) fliplr(E0(5+v,:))], light(c, 0.12), 'EdgeColor', 'none');
        fill([t fliplr(t)], sc*[E5(v,:) fliplr(E5(5+v,:))], light(c, 0.28), 'EdgeColor', 'none');
        h(j) = plot(t, sc*E1(v,:), 'Color', c, 'LineWidth', 1.5);
    end
    xlabel('Year'); ylabel(vars{k,3}); title(vars{k,4}); grid on;
    if k == 1, legend(h, LG, 'Location', 'southwest', 'Box', 'off', 'FontSize', 7.5); end
end
savefig_(f, fullfile(figdir, 'fig8_scenarios.pdf'));
fprintf('Seven figures drawn and saved as PDF in %s\n', figdir);
end

% ======================= local helper functions =======================
function c = hex(s)
c = [hex2dec(s(1:2)) hex2dec(s(3:4)) hex2dec(s(5:6))]/255;
end
function c = light(c, f)
c = 1 - f*(1 - c);
end
function y = clip01(x)
y = min(max(x, 0), 1);
end
function m = trimf(x, abc)
m = clip01(min((x - abc(1))/(abc(2) - abc(1)), (abc(3) - x)/(abc(3) - abc(2))));
end
function y = nanabove(y)
y(y > 1) = NaN;
end
function f = newfig(name, pos)
f = figure('Name', name, 'Color', 'w', 'Position', pos);
end
function h = membplot(sub, clo, chi, col)
[~, o] = sort(sub(:,3)); sub = sub(o,:);
x = 100*[sub(:,clo); flipud(sub(:,chi))]'; y = [sub(:,3); flipud(sub(:,3))]';
fill(x, y, light(col, 0.3), 'EdgeColor', 'none'); h = plot(x, y, 'Color', col, 'LineWidth', 1.4);
end
function savefig_(f, fname)
drawnow;
try
    exportgraphics(f, fname, 'ContentType', 'vector');       % MATLAB R2020a and later
catch
    set(f, 'PaperPositionMode', 'auto'); print(f, fname, '-dpdf');   % older MATLAB and Octave
end
end
