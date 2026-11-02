%% MAIN_FUZZY_MALARIA_HOOKWORM
% Intensity-dependent fuzzy model of malaria-hookworm co-infection.
% Reproduces every numerical result of the manuscript: calibration, model evaluation,
% threshold quantities, fuzzy reproduction numbers, identifiability, bifurcation and
% global-convergence checks, critical deworming coverage, sensitivity analysis and the
% fuzzy (alpha-cut) evaluation of six control scenarios for 2026-2030.
% Key numbers are printed in the Command Window. The data needed for the figures are
% saved in results/figure_data.mat and drawn by make_figures.m (run_all.m does both).
% Optional CSV tables of all results: set WRITE_CSV = true in run_all.m.
%
% Model: SEIRS malaria x SIS hookworm with soil larvae and SEI mosquitoes (src/fmh_rhs.m).
% Requires MATLAB (Optimization Toolbox for fsolve) or GNU Octave >= 8.
clear; clc;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(fullfile(here, 'src'));
outdir = fullfile(here, 'results'); if ~exist(outdir, 'dir'), mkdir(outdir); end
rng(2026);

%% 1. Core calibration (Lake Zone, Tanzania)
% Targets: malaria prevalence 22.3% (Chacky et al. 2025, Lake Zone 2023),
% hookworm prevalence 18.7% (Ndokeji et al. 2016, Ilemela), OR = 1.32 (Kinung'hi et al. 2014).
p0 = fmh_params();
tgt_core = [0.223 0.187 1.32];
[q, ystar, res] = fmh_calibrate(p0, tgt_core);
xcore = [log(q.mv); log(q.pimax*fmh_shape(q, q.eps, q.rho)); q.thetaM];
d = fmh_derived(q);
[RM, RW] = fmh_R0(q, d);
[Phi, RMW] = fmh_invasion(q, d);
[hs, Ls, yM, yW] = fmh_boundary(q, d);
o = fmh_summaries(ystar, q, d);
lam_w = q.beta_w*ystar(9)/(1 + ystar(9));
lmax_core = fmh_maxeig(q, d, ystar);
[abif, bbif, mvc] = fmh_bifcoef(q);
fprintf('Calibration residuals: %g %g %g (max|rhs| %g); max Re(eig) %g\n', res, lmax_core);
fprintf('mv=%.4f pimax=%.5f thetaM=%.4f | RM=%.4f RW=%.4f Phi=%.4f RMW=%.4f | a=%.4g b=%.4g\n', ...
        q.mv, q.pimax, q.thetaM, RM, RW, Phi, RMW, abif, bbif);
% malaria burden attributable to hookworm co-endemicity (remove hookworm: pi = 0)
q0 = q; q0.pimax = 0; d0 = fmh_derived(q0);
y0w = fmh_eqrobust(q0, d0, ystar); o0 = fmh_summaries(y0w, q0, d0);
fmh_writecsv(fullfile(outdir,'core_summary.csv'), ...
  ['mv,pimax,thetaM,RM,RW,Phi,RMW,hs,PmM,s0,e0,i0,r0,s1,e1,i1,r1,L,sv,ev,iv,Pm,Pw,Pc,OR,inc,immune,' ...
   'theta_eps,cr,gm,gw,tau,lam_w,maxeig,a_bif,b_bif,mv_crit,Pm_noW,inc_noW'], ...
  [q.mv q.pimax q.thetaM RM RW Phi RMW hs yM(3) ystar' o.Pm o.Pw o.Pc o.OR o.inc o.R ...
   d.th d.cr d.gm d.gw d.tau lam_w lmax_core abif bbif mvc o0.Pm o0.inc]);

%% 2. Validation A: national school-age malaria trend 2015-2023 (Chacky et al. 2025)
tgt_nat = [0.216 0.187 1.32];                   % mid-2015 national prevalence 21.6%
[qn, yn] = fmh_calibrate(p0, tgt_nat, xcore);
% burn-in with impulsive deworming so that the 2015 state is periodic
[~, Yb] = fmh_simulate(qn, yn, [0 8*365], fmh_schedule('annual', qn.phi, 0, 8*365), 5);
yn0 = Yb(:, end);
f2021 = @(u) fmh_hindcast_at(qn, yn0, u, 6, 6*365) - 0.118;
umax = fzero(f2021, [0 0.6]);
[tt, Pmh] = fmh_hindcast(qn, yn0, umax, 8);
P2023 = interp1(tt, Pmh, 8*365);
fprintf('Hindcast: umax=%.4f, predicted 2023 prevalence %.4f (observed 0.113)\n', umax, P2023);
fmh_writecsv(fullfile(outdir,'hindcast.csv'), 'year,Pm', [2015.5 + tt'/365, Pmh']);
fmh_writecsv(fullfile(outdir,'hindcast_summary.csv'), 'umax,P2015,P2021,P2023,mv_nat', ...
  [umax, interp1(tt,Pmh,0), interp1(tt,Pmh,6*365), P2023, qn.mv]);

%% 3. Validation B: hookworm reinfection 18 weeks after treatment (Speich et al. 2016, Pemba)
% At any endemic equilibrium lambda_w = (gw+mu)*P/(1-P) (exact, see Section 4).
alphas = 0:0.25:1;
durTFN = [1 2 3]; CRTFN = [0.50 0.72 0.81]; PpTFN = [0.377 0.425 0.474];
VB = zeros(numel(alphas), 5);
for i = 1:numel(alphas)
    al = alphas(i);
    D = fmh_alphacut(durTFN, al); C = fmh_alphacut(CRTFN, al); P = fmh_alphacut(PpTFN, al);
    lamfun = @(dur, cr, pp) (1/(dur*365) - log(1 - q.phi*cr)/365 + q.mu)*pp/(1 - pp);
    lo = lamfun(D(2), C(1), P(1)); hi = lamfun(D(1), C(2), P(2));
    VB(i,:) = [al, 1 - exp(-lo*105), 1 - exp(-hi*105), 1 - exp(-lo*63), 1 - exp(-hi*63)];
end
fmh_writecsv(fullfile(outdir,'validation_reinfection.csv'), 'alpha,R105_lo,R105_hi,R63_lo,R63_hi', VB);
disp(VB)

%% 4. Fuzzy reproduction numbers and possibility of elimination by intensity group
% Linguistic groups (triangular membership functions on epg and parasites/uL)
% Overlapping groups: adjacent memberships cross 0.5 at the WHO cut-offs (2000, 4000 epg)
Hgrp = [24 1000 3000; 1000 3000 5000; 3000 5000 8000];           % light, moderate, heavy (epg)
% parasite-density groups on log10 scale, crossing 0.5 at 1000 and 10000 parasites/uL
Mgrp = [1 2.5 3.5; 2.5 3.5 4.5; 3.5 4.5 5];                     % low, moderate, high (log10)
eg = linspace(0, 10000, 5001); rg = unique([linspace(0, 2000, 2001), linspace(2000, 100000, 9801)]);
lrg = log10(max(rg, 1e-9));
RWg = zeros(size(eg)); PhiG = zeros(size(eg)); RMg = zeros(size(rg));
for k = 1:numel(eg)
    dk = fmh_derived(q, eg(k), q.rho); [~, RWg(k)] = fmh_R0(q, dk); PhiG(k) = fmh_invasion(q, dk);
end
for k = 1:numel(rg)
    dk = fmh_derived(q, q.eps, rg(k)); RMg(k) = fmh_R0(q, dk);
end
FR = zeros(6, 7);   % [group, modal value, R(modal), Rf (FEV), Poss(R<1), Nec(R<1), centroid]
for gI = 1:3
    G = fmh_trimf(eg, Hgrp(gI,:)); Rmax = max(RWg);
    Fv = fmh_fev(RWg/Rmax, G);
    FR(gI,:) = [gI, Hgrp(gI,2), interp1(eg,RWg,Hgrp(gI,2)), Rmax*Fv, ...
                max([0, G(RWg < 1)]), 1 - max([0, G(RWg >= 1)]), sum(G.*RWg)/sum(G)];
    G = fmh_trimf(lrg, Mgrp(gI,:)); Rmax = max(RMg);
    Fv = fmh_fev(RMg/Rmax, G);
    FR(3+gI,:) = [3+gI, 10^Mgrp(gI,2), interp1(rg,RMg,10^Mgrp(gI,2)), Rmax*Fv, ...
                max([0, G(RMg < 1)]), 1 - max([0, G(RMg >= 1)]), sum(G.*RMg)/sum(G)];
end
fmh_writecsv(fullfile(outdir,'fuzzy_R0_groups.csv'), 'group,modal,R_modal,R_fuzzy,PossElim,NecElim,R_centroid', FR);
idx = 1:10:numel(eg); idr = 1:numel(rg);
fmh_writecsv(fullfile(outdir,'R_curves_eps.csv'), 'eps,RW,Phi', [eg(idx)', RWg(idx)', PhiG(idx)']);
fmh_writecsv(fullfile(outdir,'R_curves_rho.csv'), 'rho,RM', [rg(idr)', RMg(idr)']);
% critical intensities
eps_th1 = q.eps_M*(1 - q.theta0)/(q.thetaM - q.theta0);           % theta(eps)=1
iPhi = find(PhiG >= 1, 1); eps_phi1 = eg(iPhi);                     % Phi(eps)=1
iRW = find(RWg >= 1, 1); eps_RW1 = eg(iRW);                         % R_W(eps)=1
fmh_writecsv(fullfile(outdir,'critical_intensities.csv'), 'eps_theta1,eps_Phi1,eps_RW1', [eps_th1 eps_phi1 eps_RW1]);
fprintf('critical intensities: theta=1 at %.1f, Phi=1 at %.1f, RW=1 at %.1f epg\n', eps_th1, eps_phi1, eps_RW1);

% fuzzy invasion reproduction number R_M^W and possibility/necessity of facilitation (Phi > 1)
RMWg = RM*sqrt(PhiG);
FI = zeros(3, 6);   % [group, RMW(modal), RMW fuzzy (FEV), Poss(Phi>1), Nec(Phi>1), Phi fuzzy]
for gI = 1:3
    G = fmh_trimf(eg, Hgrp(gI,:));
    FI(gI,:) = [gI, interp1(eg, RMWg, Hgrp(gI,2)), max(RMWg)*fmh_fev(RMWg/max(RMWg), G), ...
                max([0, G(PhiG > 1)]), 1 - max([0, G(PhiG <= 1)]), max(PhiG)*fmh_fev(PhiG/max(PhiG), G)];
end
fmh_writecsv(fullfile(outdir,'fuzzy_invasion_groups.csv'), 'group,RMW_modal,RMW_fuzzy,PossFacil,NecFacil,Phi_fuzzy', FI);

% critical deworming coverage (Corollary): crisp curve and fuzzy alpha-cuts by intensity group
ec = linspace(0, 8000, 801); PC = zeros(numel(ec), 3);
for k = 1:numel(ec)
    dk = fmh_derived(q, ec(k), q.rho);
    PC(k,:) = [ec(k), fmh_critcov(q, dk, 365), fmh_critcov(q, dk, 182.5)];
end
PCw = PC; PCw(~isfinite(PCw)) = 2;     % coded as 2 (> 1, unattainable)
fmh_writecsv(fullfile(outdir,'critical_coverage_curve.csv'), 'eps,phic_annual,phic_semi', PCw);
FC = [];
for gI = 1:3
    for al = alphas
        I = fmh_alphacut(Hgrp(gI,:), al);
        dl = fmh_derived(q, I(1), q.rho); dh = fmh_derived(q, I(2), q.rho);
        FC = [FC; gI, al, fmh_critcov(q, dl, 365), fmh_critcov(q, dh, 365), ...
                     fmh_critcov(q, dl, 182.5), fmh_critcov(q, dh, 182.5)]; %#ok<AGROW>
    end
end
FC(~isfinite(FC)) = 2;
fmh_writecsv(fullfile(outdir,'fuzzy_critical_coverage.csv'), 'group,alpha,ann_lo,ann_hi,semi_lo,semi_hi', FC);
% exact threshold of the impulsive (periodic-round) system via the Floquet multiplier
PCI = zeros(numel(ec), 3);
for k = 1:numel(ec)
    dk = fmh_derived(q, ec(k), q.rho);
    PCI(k,:) = [ec(k), fmh_critcov_imp(q, dk, 365), fmh_critcov_imp(q, dk, 182.5)];
end
PCI(~isfinite(PCI)) = 2;
fmh_writecsv(fullfile(outdir,'critical_coverage_impulsive.csv'), 'eps,phic_imp_annual,phic_imp_semi', PCI);
fprintf('impulsive critical coverage at core eps: annual %.3f, semi-annual %.3f\n', fmh_critcov_imp(q, d, 365), fmh_critcov_imp(q, d, 182.5));
fprintf('critical coverage at core eps: annual %.3f, semi-annual %.3f\n', fmh_critcov(q, d, 365), fmh_critcov(q, d, 182.5));

%% 5. Fuzzy co-existence equilibria (extension principle over intensity groups)
ee = linspace(0, 8000, 161); EQe = zeros(numel(ee), 5);
for k = 1:numel(ee)
    dk = fmh_derived(q, ee(k), q.rho);
    yk = fmh_eqrobust(q, dk, ystar);
    ok = fmh_summaries(yk, q, dk);
    EQe(k,:) = [ee(k), ok.Pm, ok.Pw, ok.Pc, ok.inc];
end
rr = [linspace(0, 2000, 81), linspace(2050, 100000, 160)]; EQr = zeros(numel(rr), 5);
for k = 1:numel(rr)
    dk = fmh_derived(q, q.eps, rr(k));
    yk = fmh_eqrobust(q, dk, ystar);
    ok = fmh_summaries(yk, q, dk);
    EQr(k,:) = [rr(k), ok.Pm, ok.Pw, ok.Pc, ok.inc];
end
fmh_writecsv(fullfile(outdir,'equilibria_eps.csv'), 'eps,Pm,Pw,Pc,inc', EQe);
fmh_writecsv(fullfile(outdir,'equilibria_rho.csv'), 'rho,Pm,Pw,Pc,inc', EQr);
% alpha-cuts of fuzzy equilibrium prevalences for each linguistic group
FE = [];
for gI = 1:3
    for al = alphas
        I = fmh_alphacut(Hgrp(gI,:), al); V = fmh_cutvals(EQe, I);
        FE = [FE; 1, gI, al, V]; %#ok<AGROW>
        I = 10.^fmh_alphacut(Mgrp(gI,:), al); V = fmh_cutvals(EQr, I);
        FE = [FE; 2, gI, al, V]; %#ok<AGROW>
    end
end
fmh_writecsv(fullfile(outdir,'fuzzy_equilibria_alphacuts.csv'), 'disease,group,alpha,Pm_lo,Pm_hi,Pw_lo,Pw_hi,Pc_lo,Pc_hi', FE);

% invasion factor over (theta, eta) at the core hookworm equilibrium
thg = linspace(0.6, 3, 61); etg = linspace(1, 1.5, 51); PhiTE = zeros(numel(etg), numel(thg));
for i = 1:numel(etg)
    for j = 1:numel(thg)
        qq = q; qq.eta = etg(i); dd = d; dd.th = thg(j); PhiTE(i,j) = fmh_invasion(qq, dd);
    end
end
fmh_writecsv(fullfile(outdir,'Phi_theta_eta.csv'), 'eta\theta', [[NaN thg]; etg' PhiTE]);

%% 5b. Bifurcation diagram: malaria prevalence against R_M with and without hookworm
RMgrid = linspace(0.4, 2.6, 89); BD = zeros(numel(RMgrid), 6); yprev = ystar; yprev0 = y0w;
for k = 1:numel(RMgrid)
    qk = q; qk.mv = q.mv*(RMgrid(k)/RM)^2; dk = fmh_derived(qk);
    yk = fmh_eqrobust(qk, dk, yprev); ok = fmh_summaries(yk, qk, dk);
    [~, RMWk] = fmh_invasion(qk, dk);
    qk0 = qk; qk0.pimax = 0; dk0 = fmh_derived(qk0);
    [~, ~, yMk] = fmh_boundary(qk0, dk0);
    BD(k,:) = [RMgrid(k), RMWk, ok.Pm, yMk(3), fmh_maxeig(qk, dk, yk), ok.Pc];
    if ok.Pm > 1e-6, yprev = yk; end
end
fmh_writecsv(fullfile(outdir,'bifurcation.csv'), 'RM,RMW,Pm_coinf,Pm_noW,maxeig,Pc', BD);

%% 5c. Global convergence: random initial conditions in five threshold regimes
reg = struct('name', {'R0<1','EW','EM','HS','coexist'}, 'RMt', {0.8, NaN, RM, NaN, RM}, ...
             'RWt', {0.8, RW, 0.8, RW, RW}, 'RMWt', {NaN, 0.85, NaN, 1.06, NaN});
GC = []; GF = [];
tg = 0:60:50*365;
for r = 1:numel(reg)
    qr = q;
    qr.pimax = q.pimax*(reg(r).RWt/RW)*(d.gw + q.mu)/(d.gw + q.mu);
    dr = fmh_derived(qr);
    if ~isnan(reg(r).RMt)
        qr.mv = q.mv*(reg(r).RMt/RM)^2;
    else
        [~, RMWr] = fmh_invasion(qr, dr); qr.mv = q.mv*(reg(r).RMWt/RMWr)^2;
    end
    dr = fmh_derived(qr);
    [RMr, RWr] = fmh_R0(qr, dr); [~, RMWr] = fmh_invasion(qr, dr);
    for n = 1:20
        y0 = fmh_ic(rand(12,1).*[ones(8,1); 2; ones(3,1)]);
        y0(9) = 2*rand;
        [~, Y] = ode45(@(t,yy) fmh_rhs(t,yy,qr,dr), tg, y0, odeset('RelTol',1e-5,'AbsTol',1e-8));
        Pm = Y(:,3) + Y(:,7); Pw = sum(Y(:,5:8), 2);
        GC = [GC; r*ones(numel(tg),1), n*ones(numel(tg),1), tg', Pm, Pw]; %#ok<AGROW>
    end
    ye = fmh_eqrobust(qr, dr, Y(end,:)'); oe = fmh_summaries(ye, qr, dr);
    GF = [GF; r, RMr, RWr, RMWr, oe.Pm, oe.Pw, fmh_maxeig(qr, dr, ye), max(abs(Pm(end) - oe.Pm), abs(Pw(end) - oe.Pw))]; %#ok<AGROW>
end
fmh_writecsv(fullfile(outdir,'global_convergence.csv'), 'regime,ic,t,Pm,Pw', GC);
fmh_writecsv(fullfile(outdir,'global_regimes.csv'), 'regime,RM,RW,RMW,Pm_eq,Pw_eq,maxeig,dist_last', GF);
disp(GF)

%% 5d. Identifiability of the interaction modifiers and shared-exposure bracket
% For the core targets, every (eta, zeta, theta0) combination is recalibrated exactly; the
% spread of derived quantities shows what the three targets cannot identify. The
% attributable malaria incidence is computed against the hookworm-free equilibrium E_M.
etas = linspace(1, 1.3, 7); zetas = linspace(0.7, 1, 7); th0s = [0.8 0.95 1.0];
IDT = [];
for t0 = th0s
    for et = etas
        for ze = zetas
            pp = p0; pp.eta = et; pp.zeta = ze; pp.theta0 = t0;
            [qs, ys, rs] = fmh_calibrate(pp, tgt_core);
            if any(isnan(rs)), continue; end
            ds = fmh_derived(qs); os = fmh_summaries(ys, qs, ds);
            [RMs, RWs] = fmh_R0(qs, ds); [Phis, RMWs] = fmh_invasion(qs, ds);
            q0s = qs; q0s.pimax = 0; d0s = fmh_derived(q0s); [~, ~, yMs] = fmh_boundary(q0s, d0s);
            o0s = fmh_summaries(yMs, q0s, d0s);
            IDT = [IDT; t0, et, ze, ds.th, qs.mv, RMs, Phis, RMWs, 100*(os.inc - o0s.inc)/os.inc, 1/sqrt(Phis)]; %#ok<AGROW>
        end
    end
end
fmh_writecsv(fullfile(outdir,'identifiability.csv'), 'theta0,eta,zeta,theta_eps,mv,RM,Phi,RMW,attr_inc,RM_threshold', IDT);
fprintf('identifiability: theta(eps) %.3f-%.3f, Phi %.3f-%.3f, attr %.2f-%.2f%%\n', min(IDT(:,4)), max(IDT(:,4)), min(IDT(:,7)), max(IDT(:,7)), min(IDT(:,9)), max(IDT(:,9)));

%% 6. Sensitivity analysis
% (a) normalised forward sensitivity indices (central finite differences on log scale)
names = {'a','b','cmax','mv','nuh','psi','nuv','muv','ga','gs','beta_w','pimax','muL','gw0','CR','phi','theta0','thetaM','eta','zeta'};
SI = zeros(numel(names), 3);
for k = 1:numel(names)
    h = 1e-4; qp = q; qm = q;
    qp.(names{k}) = q.(names{k})*(1 + h); qm.(names{k}) = q.(names{k})*(1 - h);
    dp = fmh_derived(qp); dm = fmh_derived(qm);
    [RMp, RWp] = fmh_R0(qp, dp); [RMm, RWm] = fmh_R0(qm, dm);
    Pp = fmh_invasion(qp, dp); Pmn = fmh_invasion(qm, dm);
    SI(k,:) = [(log(RMp) - log(RMm)), (log(RWp) - log(RWm)), (log(Pp) - log(Pmn))]/(2*h);
end
fmh_writecsv(fullfile(outdir,'sensitivity_indices.csv'), 'idx,S_RM,S_RW,S_Phi', [(1:numel(names))' SI]);
if isappdata(0,'fmh_write_csv') && getappdata(0,'fmh_write_csv'), fid = fopen(fullfile(outdir,'sensitivity_names.txt'),'w'); fprintf(fid, '%s\n', names{:}); fclose(fid); end

% (b) PRCC with Latin hypercube sampling (+/-25% around core values)
pr = {'a','mv','muv','nuh','psi','gs','ga','beta_w','pimax','muL','gw0','CR','phi','theta0','thetaM','eta','zeta'};
N = 400; k = numel(pr);
U = zeros(N, k);
for j = 1:k, U(:,j) = (randperm(N)' - rand(N,1))/N; end
X = zeros(N, k); Yo = nan(N, 3);
for n = 1:N
    qq = q;
    for j = 1:k
        v = q.(pr{j})*(0.75 + 0.5*U(n,j));
        if strcmp(pr{j}, 'theta0'), v = 0.8 + 0.2*U(n,j); end
        if strcmp(pr{j}, 'eta'),    v = 1.0 + 0.3*U(n,j); end
        if strcmp(pr{j}, 'zeta'),   v = 0.7 + 0.3*U(n,j); end
        if strcmp(pr{j}, 'CR'),     v = 0.5 + 0.31*U(n,j); end
        qq.(pr{j}) = v; X(n,j) = v;
    end
    dd = fmh_derived(qq);
    yy = fmh_eqrobust(qq, dd, ystar);
    oo = fmh_summaries(yy, qq, dd);
    Yo(n,:) = [oo.Pc, oo.Pm, oo.inc];
end
ok = all(isfinite(Yo), 2);
PR = [fmh_prcc(X(ok,:), Yo(ok,1)); fmh_prcc(X(ok,:), Yo(ok,2)); fmh_prcc(X(ok,:), Yo(ok,3))]';
fmh_writecsv(fullfile(outdir,'prcc.csv'), 'idx,PRCC_Pc,PRCC_Pm,PRCC_inc', [(1:k)' PR]);
if isappdata(0,'fmh_write_csv') && getappdata(0,'fmh_write_csv'), fid = fopen(fullfile(outdir,'prcc_names.txt'),'w'); fprintf(fid, '%s\n', pr{:}); fclose(fid); end

%% 7. Integrated control scenarios under fuzzy uncertainty (alpha-cuts, 2026-2030)
% Malaria intensification (u = umax) repeats the reduction in effective biting that
% reproduced the 2015-2021 national decline (Section 2).
rng(7);
S = struct('name', {'B','D1','D2','M','I1','I2'}, ...
           'kind', {'annual','semi','semi','annual','semi','semi'}, ...
           'phi',  {0.225, 0.225, 0.75, 0.225, 0.75, 0.75}, ...
           'u',    {0, 0, 0, umax, umax, umax}, ...
           'wash', {0, 0, 0, 0, 0, 0.30});
FIN = struct('theta0',[0.80 0.95 1.00],'eta',[1.00 1.11 1.30],'zeta',[0.70 0.85 1.00], ...
             'CR',[0.50 0.72 0.81],'dur',[1 2 3],'eps',[200 600 1500], ...
             'Pm',[0.216 0.223 0.229],'Pw',[0.062 0.187 0.380],'OR',[0.98 1.32 1.78]);
fn = fieldnames(FIN); nf = numel(fn);
alv = [0 0.25 0.5 0.75 1]; ncorner = 40; nint = 12;
tgrid = 0:5:5*365; nt = numel(tgrid);
ENV = cell(numel(alv), numel(S));     % envelopes [lo; hi] for Pm, Pw, Pc, inc, MHI
OUT = [];                              % rows: alpha, scen, sample, outcomes
CAL = [];                              % calibrated parameters per sample
for ia = 1:numel(alv)
    al = alv(ia);
    if al == 1
        Z = 0.5*ones(1, nf);                       % core
    else
        Zc = double(rand(ncorner, nf) > 0.5);      % random corners of the alpha-cut box
        Zi = zeros(nint, nf);
        for j = 1:nf, Zi(:,j) = (randperm(nint)' - rand(nint,1))/nint; end
        Z = [Zc; Zi; 0.5*ones(1, nf)];
    end
    for sI = 1:numel(S), ENV{ia,sI} = [inf(5,nt); -inf(5,nt)]; end
    for n = 1:size(Z,1)
        pp = p0; tg = zeros(1,3);
        for j = 1:nf
            I = fmh_alphacut(FIN.(fn{j}), al);
            if al == 1, v = FIN.(fn{j})(2); else, v = I(1) + Z(n,j)*(I(2) - I(1)); end
            switch fn{j}
                case 'dur', pp.gw0 = 1/(v*365);
                case 'Pm',  tg(1) = v;
                case 'Pw',  tg(2) = v;
                case 'OR',  tg(3) = v;
                otherwise,  pp.(fn{j}) = v;
            end
        end
        [qs, ys, rs] = fmh_calibrate(pp, tg, xcore);
        if any(isnan(rs)) || max(abs(rs(1:3))) > 1e-6, fprintf('alpha %.2f sample %d: targets infeasible\n', al, n); continue; end
        ds = fmh_derived(qs); lmax = fmh_maxeig(qs, ds, ys);
        [RMs, RWs] = fmh_R0(qs, ds); [Phis, RMWs] = fmh_invasion(qs, ds); as = fmh_bifcoef(qs);
        CAL = [CAL; al, n, qs.mv, qs.pimax, qs.thetaM, tg, lmax, ds.th, RMs, RWs, Phis, RMWs, as]; %#ok<AGROW>
        [~, Yb] = fmh_simulate(qs, ys, [0 8*365], fmh_schedule('annual', qs.phi, 0, 8*365), 5);
        y0 = Yb(:, end);
        for sI = 1:numel(S)
            [t1, Y1] = fmh_simulate(qs, y0, [0 365], fmh_schedule('annual', qs.phi, 0, 365), 5);
            q2 = qs; q2.u = S(sI).u; q2.wash = S(sI).wash;
            [t2, Y2] = fmh_simulate(q2, Y1(:,end), [365 5*365], fmh_schedule(S(sI).kind, S(sI).phi, 365, 5*365), 5);
            Y = [Y1, Y2(:,2:end)];
            d1 = fmh_derived(qs, [], [], false); o1 = fmh_summaries(Y1, qs, d1);
            d2 = fmh_derived(q2, [], [], false); o2 = fmh_summaries(Y2(:,2:end), q2, d2);
            Pm = [o1.Pm o2.Pm]; Pw = [o1.Pw o2.Pw]; Pc = [o1.Pc o2.Pc]; inc = [o1.inc o2.inc];
            MHI = fmh_mhi(Pw, qs, qs.eps, tg(2));
            V = [Pm; Pw; Pc; inc; MHI];
            ENV{ia,sI} = [min(ENV{ia,sI}(1:5,:), V); max(ENV{ia,sI}(6:10,:), V)];
            y26 = tgrid < 365; y30 = tgrid >= 4*365;
            red = 100*(1 - mean(inc(y30))/mean(inc(y26)));
            OUT = [OUT; al, sI, n, red, mean(Pm(y30)), mean(Pw(y30)), mean(Pc(y30)), ...
                   mean(MHI(y30)), mean(inc(y30)), mean(inc(y26)), mean(MHI(y26))]; %#ok<AGROW>
        end
    end
    fprintf('alpha %.2f done (%d samples)\n', al, size(Z,1));
end
fmh_writecsv(fullfile(outdir,'scenario_samples.csv'), ...
  'alpha,scen,sample,inc_red,Pm30,Pw30,Pc30,MHI30,inc30,inc26,MHI26', OUT);
fmh_writecsv(fullfile(outdir,'calibration_samples.csv'), 'alpha,sample,mv,pimax,thetaM,Pm_t,Pw_t,OR_t,maxeig,theta_eps,RM,RW,Phi,RMW,a_bif', CAL);
for ia = 1:numel(alv)
    for sI = 1:numel(S)
        fmh_writecsv(fullfile(outdir, sprintf('env_a%d_s%d.csv', round(100*alv(ia)), sI)), ...
          't,Pm_lo,Pw_lo,Pc_lo,inc_lo,MHI_lo,Pm_hi,Pw_hi,Pc_hi,inc_hi,MHI_hi', [tgrid', ENV{ia,sI}']);
    end
end

%% 8. Data for the figures (drawn by make_figures.m)
FD = struct('q', q, 'umax', umax, 'tt', tt, 'Pmh', Pmh, 'VB', VB, ...
            'eg', eg, 'RWg', RWg, 'PhiG', PhiG, 'rg', rg, 'RMg', RMg, 'FR', FR, 'eps_RW1', eps_RW1, ...
            'EQe', EQe, 'EQr', EQr, 'FE', FE, 'PC', PCw, 'PCI', PCI, 'BD', BD, 'GC', GC, 'GF', GF, ...
            'alv', alv, 'tgrid', tgrid);
FD.ENV = ENV;
save(fullfile(outdir, 'figure_data.mat'), 'FD', '-v7');
fprintf('Results saved to %s\n', fullfile(outdir, 'figure_data.mat'));
