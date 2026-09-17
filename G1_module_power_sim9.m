%% G1 Module Power Bin Simulation vs Target
% JX / T1 Energy
% MP Solar (procured) + G2 cells, split by monthly MW plan (separate populations).
% MP cells age from Nov-26 (storage decay). Single baseline scenario.
% Module power is Gaussian by CLT, so distributions are computed analytically
% (normal PDF/CDF) — no Monte Carlo. Output: monthly bin distribution vs target.

clear; clc; close all;

%% ------------------------------------------------------------------
%  1. Cell geometry -> power conversion
%  ------------------------------------------------------------------
% G12R cell: 210 x 182 mm rectangle, cut from 272 mm pseudo-square ingot.
% 4 chamfered corners, each ~4mm x 4mm right triangle.
w_mm = 210; h_mm = 182;
corner_tri = 0.5*4*4;                 % mm^2 per corner
area_mm2   = w_mm*h_mm - 4*corner_tri;
area_cm2   = area_mm2/100;            % = 381.88 cm^2

% STC: 1000 W/m^2 = 0.1 W/cm^2. Cell power = eff * area_cm2 * 0.1
Pcell = @(eff) eff .* area_cm2 .* 0.1;   % W per cell

N_CELLS = 66;                         % cells per module
CTM     = 0.958;                      % cell-to-module ratio, from real production
                                      % (2-mo cells -> module bins; floor-edge binning)
SIGMA_MOD = 2.65;                     % module-level power std (W), fit to real production
PREMIUM_BOM = false;                  % (unused; BOM cases removed. Re-add rows to CASES to enable)
BOM_GAIN_W  = 3;                      % +3W from thicker wire / high-Tx glass (if BOM re-enabled)

%% ------------------------------------------------------------------
%  2. Efficiency distributions (discrete PMFs)
%  ------------------------------------------------------------------
% MP Solar procured cells (fixed distribution)
mp_eff  = [0.253 0.254 0.255 0.256 0.257 0.258];
mp_pmf  = [0.0063 0.0276 0.2383 0.4526 0.1626 0.0988];
mp_pmf  = mp_pmf/sum(mp_pmf);

% Screened MP distributions: procurement accepts only cells above a threshold.
% Keep the qualifying bins and renormalize.
mp_pmf_hi6 = mp_pmf .* (mp_eff >= 0.256 - 1e-9);   % >= 25.6%
mp_pmf_hi6 = mp_pmf_hi6 / sum(mp_pmf_hi6);
mp_pmf_hi8 = mp_pmf .* (mp_eff >= 0.258 - 1e-9);   % >= 25.8% (top bin only)
mp_pmf_hi8 = mp_pmf_hi8 / sum(mp_pmf_hi8);

% G2 cells: same distribution SHAPE as MP, but shifted so the MEDIAN
% tracks the monthly ramp target. We shift the discrete support by the
% delta between the target median and the MP baseline median.
mp_median = weighted_median(mp_eff, mp_pmf);   % ~0.256

%% ------------------------------------------------------------------
%  3. Monthly plan and G2 availability (MW)
%  ------------------------------------------------------------------
months     = {'Nov26','Dec26','Jan','Feb','Mar','Apr','May','Jun', ...
              'Jul','Aug','Sep','Oct','Nov','Dec'};
% G1 module plan (MW). Nov/Dec 2026 use ~460 MW baseline (production start).
g1_plan_MW = [460 460 462 414 478 462 478 462 478 478 462 478 446 466];
% G2 cells (MW). None available in 2026.
g2_MW      = [  0   0   0   0   0 139 204 170 172 175 178 225 182 195];
nM = numel(months);

% Months since cell production (Nov-26). Nov-26 = 0 (fresh), Dec-26 = 1, etc.
% This is the storage age for cells consumed each month.
month_age  = 0:(nM-1);

%% ------------------------------------------------------------------
%  4. G2 efficiency ramp (median)
%  ------------------------------------------------------------------
% G2 median ramp (single ramp for all cases): 25.3% Apr -> 25.6% Dec.
g2_med_ramp = ramp_median(months, {'Apr',0.253}, {'Dec',0.256});

%% ------------------------------------------------------------------
%  5. Storage decay on MP cells (absolute efficiency drop, month-anchored)
%  ------------------------------------------------------------------
% Decay is an ABSOLUTE cell-efficiency drop (Δη) tied to production MONTH.
% Expert best & worst estimates; mid-point (their average) drives the sim.
%   Best  -0.10 %pt @Apr, -0.15 %pt @Dec   (Δη 0.0010 -> 0.0015)
%   Worst -0.20 %pt @Apr, -0.25 %pt @Dec   (Δη 0.0020 -> 0.0025)
%   Mid   -0.15 %pt @Apr, -0.20 %pt @Dec   (Δη 0.0015 -> 0.0020)
% Δη = 0 at cell production start (Nov-26), ramping to the Apr anchor (kink at
% Apr), then extrapolated to Dec. Real data informs Jan->Apr; Apr->Dec is a guess.
% DECAY_SPLIT_MONTH is the array index of April (real/extrapolated boundary).
% With Nov26,Dec26,Jan..Dec, April is the 6th element.
DECAY_SPLIT_MONTH = 6;                  % Apr: boundary of real vs extrapolated
dEff_best  = decay_anchored(month_age, 0.0010, 0.0015);
dEff_worst = decay_anchored(month_age, 0.0020, 0.0025);
dEff_mid   = decay_anchored(month_age, 0.0015, 0.0020);   % drives the simulation

%% ------------------------------------------------------------------
%  6. Target module-power bin distribution (client)
%  ------------------------------------------------------------------
% bins (W) and fractions from client data. 605 added as an explicit
% target bin at 0% (client did not request 605, but sim can land there).
tgt_bin  = [600 605 610 615 620];
tgt_cnt  = [613 0 475501 1957657 578847];        % 600,605,610,615,620 counts
tgt_frac = tgt_cnt/sum(tgt_cnt);

%% ------------------------------------------------------------------
%  7. Analytic module-power distributions per month — 2 cases
%  ------------------------------------------------------------------
% Each module = sum of 66 i.i.d. cell powers + a module-level Gaussian.
% By the CLT the module power is Gaussian, so no Monte Carlo is needed:
%   mean = 66 * mean(cell power) * decay * CTM   (+BOM)
%   var  = 66 * var(cell power) * (decay*CTM)^2  + SIGMA_MOD^2
% We keep MP and G2 as separate populations (weights = MW share) so the MP
% distribution is independent of G2 efficiency.

% case definitions: name, MP-screen ('none'|'hi6'|'hi8')
% All cases share the single G2 ramp and single mid-point storage decay.
CASES = { ...
    'base',  'none'; ...   % MP as-is (full distribution)
    'mp256', 'hi6' };      % MP screened >= 25.6%
nCase = size(CASES,1);
results = struct();

% floor bins
binLbl   = [600 605 610 615 620];      % bin labels (W)
binEdges = [binLbl 625];               % floor edges: [600,605),...,[620,625)


for cidx = 1:nCase
    scen    = CASES{cidx,1};
    screen  = CASES{cidx,2};
    g2_med  = g2_med_ramp;                          % shared G2 ramp
    mp_dEff = dEff_mid;                             % shared mid-point decay
    bom_add = 0;
    switch screen
        case 'hi6', mpPMF = mp_pmf_hi6;
        case 'hi8', mpPMF = mp_pmf_hi8;
        otherwise,  mpPMF = mp_pmf;
    end

    muMP = zeros(nM,1); sdMP = zeros(nM,1);   % pure-MP module Gaussian
    muG2 = zeros(nM,1); sdG2 = zeros(nM,1);   % pure-G2 module Gaussian
    wG2  = zeros(nM,1);                        % G2 population weight (MW share)
    monthAvg = zeros(nM,1); monthStd = zeros(nM,1);
    monthDistMP = zeros(nM,numel(binLbl));
    monthDistG2 = zeros(nM,numel(binLbl));
    monthDist   = zeros(nM,numel(binLbl));

    for m = 1:nM
        g2_share = min(max(g2_MW(m)/g1_plan_MW(m),0),1);
        wG2(m) = g2_share;

        % --- MP module Gaussian: storage decay drops each cell efficiency by
        %     dEff (absolute), so recompute MP cell-power moments at reduced eff.
        %     Uses screened MP distribution when screenMP is set.
        pc_mpD  = Pcell(mp_eff - mp_dEff(m));
        muC_mpD = sum(pc_mpD .* mpPMF);
        varC_mpD= sum((pc_mpD - muC_mpD).^2 .* mpPMF);
        muMP(m) = N_CELLS*muC_mpD*CTM + bom_add;
        sdMP(m) = sqrt(N_CELLS*varC_mpD*(CTM)^2 + SIGMA_MOD^2);

        % --- G2 module Gaussian (shifted efficiency, no decay) ---
        % G2 always uses the full (unscreened) shape; only MP is screened.
        pc_g2   = Pcell(mp_eff + (g2_med(m) - mp_median));
        muC_g2  = sum(pc_g2 .* mp_pmf);
        varC_g2 = sum((pc_g2 - muC_g2).^2 .* mp_pmf);
        muG2(m) = N_CELLS*muC_g2*CTM + bom_add;
        sdG2(m) = sqrt(N_CELLS*varC_g2*(CTM)^2 + SIGMA_MOD^2);

        % --- floor-bin fractions via normal CDF, weighted MP/G2 ---
        cdfMP = normcdf(binEdges, muMP(m), sdMP(m));
        cdfG2 = normcdf(binEdges, muG2(m), sdG2(m));
        fMP = diff(cdfMP);  fG2 = diff(cdfG2);
        monthDistMP(m,:) = (1-g2_share)*fMP;
        monthDistG2(m,:) =    g2_share *fG2;
        monthDist(m,:)   = monthDistMP(m,:) + monthDistG2(m,:);

        % --- mixture mean / std (for the row labels) ---
        w = [1-g2_share, g2_share];
        mus = [muMP(m), muG2(m)];  sds = [sdMP(m), sdG2(m)];
        mix_mu = sum(w.*mus);
        mix_var= sum(w.*(sds.^2 + mus.^2)) - mix_mu^2;
        monthAvg(m) = mix_mu;
        monthStd(m) = sqrt(mix_var);
    end

    results.(scen).muMP = muMP; results.(scen).sdMP = sdMP;
    results.(scen).muG2 = muG2; results.(scen).sdG2 = sdG2;
    results.(scen).wG2  = wG2;
    results.(scen).dist   = monthDist;
    results.(scen).distMP = monthDistMP;
    results.(scen).distG2 = monthDistG2;
    results.(scen).avg  = monthAvg;
    results.(scen).std  = monthStd;
    results.(scen).lbl  = binLbl;
end

%% ------------------------------------------------------------------
%  8. Plots
%  ------------------------------------------------------------------
% Target fractions mapped onto the same 5 bin labels [600 605 610 615 620]
tgt_map = zeros(1,numel(binLbl));
for i = 1:numel(tgt_bin)
    j = find(binLbl==tgt_bin(i));
    if ~isempty(j), tgt_map(j) = tgt_map(j) + tgt_frac(i); end
end

% Shared plot constants (used by the bin plot below).
caseList = fieldnames(results);
ROW = 0.9;                                   % vertical spacing between months
SCL = 2.0;                                   % height scale (unused by bin plot)
cMP = [0 0 0]; cG2 = [0.871 0.996 0.424];   % MP=black (procured), G2=#DEFE6C green

% Module Power Bin vs Month — RIDGELINE over 5W bins, one axes per case.

% Each month is a stacked step-area over the 5W bins [600 605 610 615 620],
% MP (blue) + G2 (orange), offset vertically. One croppable graphic per case.
% 5W bins concentrate mass into single tall bars, so use a smaller vertical
% scale than (a) to keep the 615 bar within its row.
binX = binLbl;                               % 600 605 610 615 620
SCL_B = 0.75;                                % height scale for 5W-bin bars
mStart = 2;                                  % plot from Dec-26 (index 2); Nov-26 stays in model
nPlot = nM - mStart + 1;                     % rows shown
for ci = 1:numel(caseList)
  scen = caseList{ci};
  for showG2 = [true false]                  % full (MP+G2), then MP-only
    if showG2, tag=''; else, tag=' (MP only)'; end
    figure('Name',['G1 5W-bin ridgeline - ' scen tag],'Color','w', ...
           'Position',[780+40*ci 60 720 900]);
    ax = axes; hold(ax,'on'); set(ax,'Color','w');
    for m = nM:-1:mStart
        base = (nM - m) * ROW;
        dmp = results.(scen).distMP(m,:);    % MP fraction per bin (prod-weighted)
        dg2 = results.(scen).distG2(m,:);    % G2 fraction per bin
        if ~showG2
            % MP-only: show pure MP distribution (renormalized to sum 1), no G2
            s = sum(dmp); if s>0, dmp = dmp/s; end
            dg2 = zeros(size(dg2));
        end
        bw = 4.2;                            % bar width (of 5W bin)
        for b = 1:numel(binX)
            x0 = binX(b)-bw/2; x1 = binX(b)+bw/2;
            yb = base;
            yMP = base + SCL_B*dmp(b);
            yTP = base + SCL_B*(dmp(b)+dg2(b));
            if dmp(b)>0
                fill([x0 x1 x1 x0],[yb yb yMP yMP],cMP,'EdgeColor','none');
            end
            if dg2(b)>0
                fill([x0 x1 x1 x0],[yMP yMP yTP yTP],cG2,'EdgeColor','k','LineWidth',0.5);
            end
        end
        text(596.5, base+0.12, sprintf('Month %d', m-mStart+1), ...
             'HorizontalAlignment','right','FontSize',9,'VerticalAlignment','bottom', ...
             'Color','k');
    end
    xlim([596 623]); ylim([-0.1 (nM-mStart)*ROW + SCL_B*1.0 + 0.2]);
    set(ax,'XTick',binLbl,'XTickLabel',["600" "605" "610" "615" "620"], ...
        'YTick',[],'XColor','k','YColor','none','FontSize',9, ...
        'TickLabelInterpreter','tex');
    xlabel('Power bin (W)','FontSize',10,'Color','k'); box off;
    if showG2
        title(sprintf('G1 Module Power Bin by Cell Source vs Month — %s', ...
              caseLabel(scen)),'FontSize',11,'Color','k');
    else
        title(sprintf('G1 Module Power Bin (MP cells only) vs Month — %s', ...
              caseLabel(scen)),'FontSize',11,'Color','k');
    end
    % legend swatches
    yl = ylim; xL = 619; yb2 = yl(2)-0.35;
    fill(xL+[0 1.0 1.0 0], yb2+[0 0 0.12 0.12], cMP,'EdgeColor','none');
    text(xL+1.3, yb2+0.06,'MP','FontSize',8,'VerticalAlignment','middle','Color','k');
    if showG2
        fill(xL+[0 1.0 1.0 0], yb2-0.18+[0 0 0.12 0.12], cG2,'EdgeColor','k','LineWidth',0.5);
        text(xL+1.3, yb2-0.12,'G2','FontSize',8,'VerticalAlignment','middle','Color','k');
    end
  end
end

%% ------------------------------------------------------------------
%  8c. Input-assumption plots: G2 ramp & MP storage decay
%  ------------------------------------------------------------------
xm = 1:nM;

% (c) G2 efficiency ramp (G2 production starts Apr)
figure('Name','G2 efficiency ramp','Color','w', ...
       'Position',[300 120 640 380]);
ax = axes; hold(ax,'on'); set(ax,'Color','w');
g2start = find(g2_MW>0, 1);                   % first month with G2 (Apr-27)
idx = g2start:nM;
cG2dark = [0.42 0.55 0.12];                    % darker shade of G2 green for line-on-white
plot(xm(idx), 100*g2_med_ramp(idx), '-o','Color',cG2dark,'MarkerFaceColor',cG2dark, ...
     'LineWidth',1.6,'MarkerSize',5);
% mark the ramp anchors (Apr 25.3%, then climbing toward Dec 25.6%)
yline(25.3,'--','Color',[0.6 0.6 0.6]);
xlim([g2start-0.5 nM+0.5]); ylim([25.2 25.7]);
set(ax,'XTick',idx,'XTickLabel',months(idx),'XColor','k','YColor','k','FontSize',9);
xlabel('Month (2027)','FontSize',10,'Color','k');
ylabel('G2 median cell efficiency (%)','FontSize',10,'Color','k');
title('G2 Efficiency Ramp (25.3% Apr \rightarrow 25.6% Dec)', ...
      'FontSize',11,'Color','k');
grid on; box off;

% (d) MP cell storage decay -> cell efficiency drop (% absolute).
% Best, worst, and mid-point curves. Mid drives the sim (bold).
% Jan-Apr informed by data (solid); Apr-Dec extrapolated (dotted). Kink at Apr.
figure('Name','MP storage decay','Color','w', ...
       'Position',[360 120 680 410]);
ax = axes; hold(ax,'on'); set(ax,'Color','w');
s = DECAY_SPLIT_MONTH;                         % Apr = real/extrapolated boundary
xr = 1:s;  xe = s:nM;
cWorst=[0.75 0.2 0.2]; cBest=[0.42 0.55 0.12]; cMid=[0 0 0];
dW=-100*dEff_worst; dB=-100*dEff_best; dM=-100*dEff_mid;
% worst & best: thin reference curves
plot(xr,dW(xr),'-','Color',cWorst,'LineWidth',1.0); plot(xe,dW(xe),':','Color',cWorst,'LineWidth',1.0);
plot(xr,dB(xr),'-','Color',cBest,'LineWidth',1.0);  plot(xe,dB(xe),':','Color',cBest,'LineWidth',1.0);
% mid: bold, drives the simulation
plot(xr,dM(xr),'-o','Color',cMid,'MarkerFaceColor',cMid,'LineWidth',1.8,'MarkerSize',5);
plot(xe,dM(xe),':o','Color',cMid,'MarkerFaceColor',cMid,'LineWidth',1.8,'MarkerSize',5);
yline(0,'-','Color',[0.75 0.75 0.75]);
xline(s,'--','Color',[0.7 0.7 0.7]);
text(s+0.1, 0.02, ' extrapolated \rightarrow','Color',[0.45 0.45 0.45], ...
     'FontSize',8,'VerticalAlignment','top');
% legend (proxy handles)
hW=plot(nan,nan,'-','Color',cWorst,'LineWidth',1.0);
hM=plot(nan,nan,'-o','Color',cMid,'LineWidth',1.8,'MarkerFaceColor',cMid);
hB=plot(nan,nan,'-','Color',cBest,'LineWidth',1.0);
lgD = legend([hW hM hB],{'Worst','Mid (used)','Best'},'Location','southwest','FontSize',8,'Box','off');
set(lgD,'TextColor','k');
xlim([0.5 nM+0.5]); ylim([min(dW)*1.2 0.06]);
set(ax,'XTick',xm,'XTickLabel',months,'XColor','k','YColor','k','FontSize',9);
xlabel('Month (2027)','FontSize',10,'Color','k');
ylabel('MP cell efficiency drop (% absolute)','FontSize',10,'Color','k');
title('MP Cell Efficiency Drop from Storage (solid = data, dotted = extrapolated)', ...
      'FontSize',10.5,'Color','k');
grid on; box off;

% (e) Monthly MP (procured) cell consumption in MW = G1 plan - G2.
% MP fills the full plan until G2 kicks in (Apr), then MP = G1 - G2.
figure('Name','MP cell consumption','Color','w', ...
       'Position',[420 120 720 380]);
ax = axes; hold(ax,'on'); set(ax,'Color','w');
mp_MW = max(g1_plan_MW - g2_MW, 0);           % MP fills the remainder
hbar = bar(xm, [mp_MW(:) g2_MW(:)], 'stacked', 'BarWidth',0.7);
hbar(1).FaceColor = [0 0 0];                   % MP = black
hbar(1).EdgeColor = 'none';
hbar(2).FaceColor = cG2;                       % G2 = green
hbar(2).EdgeColor = 'k'; hbar(2).LineWidth = 0.4;
xlim([0.5 nM+0.5]); ylim([0 max(g1_plan_MW)*1.15]);
set(ax,'XTick',xm,'XTickLabel',months,'XColor','k','YColor','k','FontSize',9);
xtickangle(ax,45);
xlabel('Month','FontSize',10,'Color','k');
ylabel('Cell consumption (MW)','FontSize',10,'Color','k');
lgC = legend({'MP (procured)','G2'},'Location','northwest','FontSize',8,'Box','off');
set(lgC,'TextColor','k');
title('Monthly Cell Consumption — MP fills plan until G2 ramps (Apr-27)', ...
      'FontSize',10.5,'Color','k');
grid on; box off;

%% ------------------------------------------------------------------
%  9. Summary table — average module power, all cases
%  ------------------------------------------------------------------
fprintf('\n=== Monthly avg module power (W), sigma_mod=%.1fW, mid-point decay ===\n',SIGMA_MOD);
fprintf('MP screen keeps: >=25.6%% -> %.0f%% of MP cells.\n', ...
        100*sum(mp_pmf(mp_eff>=0.256-1e-9)));
fprintf('%-6s %9s %9s %7s\n','Month','MP all','MP>=25.6','G2sh');
for m=1:nM
    fprintf('%-6s %9.1f %9.1f %6.0f%%\n', months{m}, ...
        results.base.avg(m), results.mp256.avg(m), ...
        100*g2_MW(m)/g1_plan_MW(m));
end
fprintf('\nTarget avg ~615W (bins: 16%% @610, 65%% @615, 19%% @620; 605 target 0%%).\n');

%% ==================================================================
%  10. CELL CONSUMPTION STRATEGY (alpha: low-first <-> uniform)
%  ------------------------------------------------------------------
%  Two inventory scenarios: (a) MP all bins, with <25.6% grouped as one low
%  tier, then 25.6/25.7/25.8; (b) MP screened >=25.6% (3 tiers). Total = MP
%  demand Dec26->Dec. Each month's consumption blends two strategies by ALPHA
%  (applied to ALL months):
%      ALPHA=0  -> low-first  (consume lowest available tier, hoard high bins)
%      ALPHA=1  -> uniform    (consume all tiers proportional to remaining stock)
%      0<ALPHA<1 -> blend.  used = (1-ALPHA)*lowfirst + ALPHA*uniform, then
%      clamped to inventory with any shortfall topped up low-first.
%  Low ALPHA protects the year-end (banks high bins vs worst decay) but sags
%  mid-year; high ALPHA lifts early months but lets year-end decline.
%  Consumed cells carry storage decay, then lock in (laminate).

ALPHA_ALL = 0.8;                              % all-bins scenario (0=low-first, 1=uniform)
ALPHA_HI6 = 0.9;                              % >=25.6% scenario  (0=low-first, 1=uniform)
ALPHA_LOW = 0.9;                              % low-only (<25.6%) scenario

mp_demand = max(g1_plan_MW - g2_MW, 0);       % MP MW per month
total_mp  = sum(mp_demand(2:nM));             % Dec26 onward (exclude Nov26)

% --- Define the two inventory scenarios (tiers low->high) ---
% All-bins: group 25.3/25.4/25.5 into one low tier, then 25.6, 25.7, 25.8.
lowMask = mp_eff < 0.256;
lowFrac = sum(mp_pmf(lowMask));
lowEff  = sum(mp_eff(lowMask).*mp_pmf(lowMask))/lowFrac;
S = struct();
S.all.eff  = [lowEff 0.256 0.257 0.258];
S.all.frac = [lowFrac mp_pmf(mp_eff==0.256) mp_pmf(mp_eff==0.257) mp_pmf(mp_eff==0.258)];
S.all.frac = S.all.frac/sum(S.all.frac);
S.all.name = 'MP all bins';
S.all.leg  = {'<25.6% (grp)','25.6%','25.7%','25.8%'};
S.all.alpha = ALPHA_ALL;
S.all.splitStart = 1;                          % all-bins: alpha from month 1 (all months)
hi = mp_pmf(mp_eff>=0.256); hi = hi/sum(hi);
S.hi6.eff  = [0.256 0.257 0.258];
S.hi6.frac = hi;
S.hi6.name = 'MP \geq25.6%';
S.hi6.leg  = {'25.6%','25.7%','25.8%'};
S.hi6.alpha = ALPHA_HI6;
S.hi6.splitStart = 4;                          % >=25.6%: pure low-first months 1-3, alpha from Mar

% Low-only: procure only cells below 25.6% (25.3/25.4/25.5), 3 separate tiers.
% Inventory = full MP demand distributed by the low bins' relative proportions.
loMask = mp_eff < 0.256;
S.low.eff  = mp_eff(loMask);                    % [0.253 0.254 0.255]
S.low.frac = mp_pmf(loMask)/sum(mp_pmf(loMask));
S.low.name = 'MP low-only (<25.6%)';
S.low.leg  = {'25.3%','25.4%','25.5%'};
S.low.alpha = ALPHA_LOW;
S.low.splitStart = 1;                          % alpha from month 1

stratScen = {'all','hi6','low'};
STR = struct();                                % strategy results per scenario

for si = 1:numel(stratScen)
    sc = stratScen{si};
    tierEff = S.(sc).eff; tierFrac = S.(sc).frac; nT = numel(tierEff);
    ALPHA = S.(sc).alpha;                       % this scenario's alpha
    splitStart = S.(sc).splitStart;             % plot-month where alpha blend begins
    invR = tierFrac * total_mp;                 % MW per tier (low->high)
    mpP = nan(nM,1); blend = nan(nM,1);
    used_all = zeros(nM,nT);
    dist = zeros(nM,numel(binLbl));
    for m = 2:nM
        demand = mp_demand(m); if demand<=0, continue; end
        de = dEff_mid(m);
        % low-first candidate
        L = zeros(1,nT); need = demand;
        for t = 1:nT, x=min(invR(t),need); L(t)=x; need=need-x; end
        % uniform candidate (proportional to remaining stock)
        tot = sum(invR);
        if tot>0, U = min(demand*invR/tot, invR); else, U = zeros(1,nT); end
        pm = m - 1;                             % plot-month (Dec26 = 1)
        if pm < splitStart
            used = L;                           % pure low-first (early months)
        else
            used = min((1-ALPHA)*L + ALPHA*U, invR);   % alpha blend
        end
        % top up shortfall low-first
        short = demand - sum(used); rem = invR - used;
        for t = 1:nT
            if short<=0, break; end
            x=min(rem(t),short); used(t)=used(t)+x; rem(t)=rem(t)-x; short=short-x;
        end
        invR = invR - used;
        used_all(m,:) = used; um = sum(used);
        if um>0
            mpP(m) = N_CELLS * sum((tierEff-de).*used)/um * area_cm2 * 0.1 * CTM;
            wg = g2_MW(m)/g1_plan_MW(m);
            g2P = N_CELLS * g2_med_ramp(m) * area_cm2 * 0.1 * CTM;
            blend(m) = (1-wg)*mpP(m) + wg*g2P;
            for t = 1:nT
                w = used(t)/um;
                if w>0
                    mu = N_CELLS*(tierEff(t)-de)*area_cm2*0.1*CTM;
                    dist(m,:) = dist(m,:) + w*diff(normcdf(binEdges, mu, SIGMA_MOD));
                end
            end
        end
    end
    STR.(sc).mpP = mpP; STR.(sc).blend = blend; STR.(sc).used = used_all;
    STR.(sc).dist = dist; STR.(sc).invR = invR; STR.(sc).inv0 = tierFrac*total_mp;
    STR.(sc).eff = tierEff; STR.(sc).leg = S.(sc).leg; STR.(sc).name = S.(sc).name;
    STR.(sc).alpha = ALPHA;
end

mm = 2:nM;
mLabels = arrayfun(@(k) sprintf('Month %d',k), mm-1, 'UniformOutput', false);
posX = [380 460 540];

% ---- Plot: strategy monthly module power (both scenarios) ----
for si = 1:numel(stratScen)
    sc = stratScen{si};
    figure('Name',['Strategy power - ' sc],'Color','w', ...
           'Position',[posX(si) 90 760 420]);
    ax = axes; hold(ax,'on'); set(ax,'Color','w');
    plot(mm, STR.(sc).mpP(mm), '-o','Color',[0 0 0],'MarkerFaceColor',[0 0 0], ...
         'LineWidth',1.7,'MarkerSize',5);
    plot(mm, STR.(sc).blend(mm), '-s','Color',cG2dark,'MarkerFaceColor',cG2dark, ...
         'LineWidth',1.7,'MarkerSize',5);
    yline(615,'--','Color',[0.5 0.5 0.5]);
    text(nM, 615, ' 615W target','Color',[0.4 0.4 0.4],'FontSize',8, ...
         'VerticalAlignment','bottom','HorizontalAlignment','right');
    xlim([1.5 nM+0.5]); ylim([610 620]);
    set(ax,'XTick',mm,'XTickLabel',mLabels,'XColor','k','YColor','k','FontSize',9);
    xtickangle(ax,45);
    xlabel('Month','FontSize',10,'Color','k');
    ylabel('Module power (W)','FontSize',10,'Color','k');
    lg = legend({'MP only','MP + G2 blended'},'Location','southwest','FontSize',8,'Box','off');
    set(lg,'TextColor','k');
    title(sprintf('Consumption Strategy (\\alpha=%.2f) — %s', STR.(sc).alpha, STR.(sc).name), ...
          'FontSize',10,'Color','k');
    grid on; box off;
end

% ---- Plot: bin consumption per month (stacked MW), both scenarios ----
cTier4 = [0.55 0.55 0.55; 0.2 0.2 0.2; 0.55 0.65 0.30; 0.871 0.996 0.424];  % grp,25.6,25.7,25.8
for si = 1:numel(stratScen)
    sc = stratScen{si};
    nT = numel(STR.(sc).eff);
    if nT==4, cols = cTier4; else, cols = cTier4(2:4,:); end
    figure('Name',['Strategy bin consumption - ' sc],'Color','w', ...
           'Position',[posX(si)+40 90 760 400]);
    ax = axes; hold(ax,'on'); set(ax,'Color','w');
    hb = bar(mm, STR.(sc).used(mm,:), 'stacked', 'BarWidth',0.7);
    for k=1:nT, hb(k).FaceColor = cols(k,:); hb(k).EdgeColor='k'; hb(k).LineWidth=0.3; end
    xlim([1.5 nM+0.5]);
    set(ax,'XTick',mm,'XTickLabel',mLabels,'XColor','k','YColor','k','FontSize',9);
    xtickangle(ax,45);
    xlabel('Month','FontSize',10,'Color','k');
    ylabel('MP cells consumed (MW)','FontSize',10,'Color','k');
    lgB = legend(STR.(sc).leg,'Location','northeast','FontSize',8,'Box','off');
    set(lgB,'TextColor','k');
    title(sprintf('Strategy Cell Consumption by Bin — %s', STR.(sc).name), ...
          'FontSize',10.5,'Color','k');
    grid on; box off;
end

% ---- Plot: bin consumption per month as % of month (100% stacked) ----
for si = 1:numel(stratScen)
    sc = stratScen{si};
    nT = numel(STR.(sc).eff);
    if nT==4, cols = cTier4; else, cols = cTier4(2:4,:); end
    % normalize each month's consumption to % of that month's total
    U = STR.(sc).used(mm,:);
    rowsum = sum(U,2); rowsum(rowsum==0) = 1;
    Upct = 100 * U ./ rowsum;
    figure('Name',['Strategy bin share %% - ' sc],'Color','w', ...
           'Position',[posX(si)+120 90 760 400]);
    ax = axes; hold(ax,'on'); set(ax,'Color','w');
    hb = bar(mm, Upct, 'stacked', 'BarWidth',0.7);
    for k=1:nT, hb(k).FaceColor = cols(k,:); hb(k).EdgeColor='k'; hb(k).LineWidth=0.3; end
    % percentage labels centered in each segment (skip <4%)
    for r = 1:numel(mm)
        cum = 0;
        for k = 1:nT
            v = Upct(r,k);
            if v >= 4
                txtCol = 'k'; if k<=1 && nT==4, txtCol='w'; end   % grp tier is mid-grey
                if k==2 && nT==4, txtCol='w'; end                 % 25.6 dark -> white
                if k==1 && nT==3, txtCol='w'; end                 % 25.6 dark -> white
                text(mm(r), cum+v/2, sprintf('%.0f',v), 'HorizontalAlignment','center', ...
                     'VerticalAlignment','middle','FontSize',6.5,'Color',txtCol);
            end
            cum = cum + v;
        end
    end
    xlim([1.5 nM+0.5]); ylim([0 100]);
    set(ax,'XTick',mm,'XTickLabel',mLabels,'XColor','k','YColor','k','FontSize',9);
    xtickangle(ax,45);
    xlabel('Month','FontSize',10,'Color','k');
    ylabel('Share of month''s cells (%)','FontSize',10,'Color','k');
    lgP = legend(STR.(sc).leg,'Location','eastoutside','FontSize',8,'Box','off');
    set(lgP,'TextColor','k');
    title(sprintf('Strategy Bin Share per Month — %s', STR.(sc).name), ...
          'FontSize',10.5,'Color','k');
    box off;
end

% ---- Plot: strategy month-by-month power bin ridgeline, both scenarios ----
for si = 1:numel(stratScen)
    sc = stratScen{si};
    figure('Name',['Strategy bin ridgeline - ' sc],'Color','w', ...
           'Position',[posX(si)+80 60 720 900]);
    ax = axes; hold(ax,'on'); set(ax,'Color','w');
    for m = nM:-1:2
        base = (nM - m) * ROW; d = STR.(sc).dist(m,:); bw = 4.2;
        for b = 1:numel(binLbl)
            x0 = binLbl(b)-bw/2; x1 = binLbl(b)+bw/2; yTP = base + SCL_B*d(b);
            if d(b)>0
                fill([x0 x1 x1 x0],[base base yTP yTP],cMP,'EdgeColor','none');
                pct = 100*d(b);
                if pct >= 3
                    text(binLbl(b), yTP+0.02, sprintf('%.0f%%',pct), ...
                         'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                         'FontSize',6.5,'Color','k');
                end
            end
        end
        text(596.5, base+0.12, sprintf('Month %d', m-1), ...
             'HorizontalAlignment','right','FontSize',9,'VerticalAlignment','bottom','Color','k');
    end
    xline(615,'--','Color',[0.5 0.5 0.5]);
    xlim([596 623]); ylim([-0.1 (nM-2)*ROW + SCL_B*1.0 + 0.2]);
    set(ax,'XTick',binLbl,'XTickLabel',["600" "605" "610" "615" "620"], ...
        'YTick',[],'XColor','k','YColor','none','FontSize',9,'TickLabelInterpreter','tex');
    xlabel('Power bin (W)','FontSize',10,'Color','k'); box off;
    title(sprintf('Strategy Monthly Module Power Bin — %s', STR.(sc).name), ...
          'FontSize',10.5,'Color','k');
end

% ---- Console summary (both scenarios) ----
for si = 1:numel(stratScen)
    sc = stratScen{si};
    fprintf('\n=== Consumption strategy (%s, alpha=%.2f) ===\n', STR.(sc).name, STR.(sc).alpha);
    fprintf('%-6s %6s %8s %8s\n','Month','MP-MW','MP-W','blend-W');
    for m=2:nM
        fprintf('%-6s %6.0f %8.1f %8.1f\n', months{m}, mp_demand(m), ...
                STR.(sc).mpP(m), STR.(sc).blend(m));
    end
    % overall bin composition across all simulated months
    tierTot = sum(STR.(sc).used(2:nM,:), 1);    % total MW consumed per tier
    tierPct = 100 * tierTot / sum(tierTot);
    fprintf('Overall bin composition over all months (%s):\n', STR.(sc).name);
    labels = STR.(sc).leg;
    for t = 1:numel(labels)
        fprintf('   %-14s %6.0f MW  (%5.1f%%)\n', labels{t}, tierTot(t), tierPct(t));
    end
    fprintf('   %-14s %6.0f MW\n', 'TOTAL', sum(tierTot));
    % overall MODULE POWER BIN distribution across all months,
    % weighted by each month's MP production volume.
    wts = mp_demand(2:nM); wts = wts(:);
    dmat = STR.(sc).dist(2:nM,:);               % months x bins
    binOverall = sum(dmat .* wts, 1) / sum(wts);
    binOverall = 100 * binOverall / sum(binOverall);
    fprintf('Overall module power-bin distribution (%s):\n', STR.(sc).name);
    for b = 1:numel(binLbl)
        fprintf('   %dW bin: %5.1f%%\n', binLbl(b), binOverall(b));
    end
    % volume-weighted average module power over all months
    mpAvg    = sum(STR.(sc).mpP(2:nM)   .* wts) / sum(wts);   % MP-only
    blendAvg = sum(STR.(sc).blend(2:nM) .* wts) / sum(wts);   % MP + G2 blended
    fprintf('Avg module power (%s): MP-only = %.1f W, blended = %.1f W\n', ...
            STR.(sc).name, mpAvg, blendAvg);
end


function s = caseLabel(scen)
    switch scen
        case 'base',  s = 'MP all bins';
        case 'mp256', s = 'MP screened \geq25.6%';
        otherwise,    s = upper(scen);
    end
end

function y = normpdf(x, mu, s)
    y = exp(-0.5*((x-mu)./s).^2) ./ (s*sqrt(2*pi));
end

function p = normcdf(x, mu, s)
    p = 0.5*(1 + erf((x-mu)./(s*sqrt(2))));
end

function med = weighted_median(vals, pmf)
    c = cumsum(pmf);
    k = find(c>=0.5,1);
    med = vals(k);
end

function medvec = ramp_median(months, startPt, endPt)
    % Interpolate against ABSOLUTE month index so endpoints (e.g. Dec-27)
    % outside the plotted month list still work.
    i0 = mabs(startPt{1}); v0 = startPt{2};
    i1 = mabs(endPt{1});   v1 = endPt{2};
    nM = numel(months);
    medvec = nan(1,nM);
    for m=1:nM
        mi = mabs(months{m});
        if mi < i0
            medvec(m) = v0;                       % before start: hold
        elseif mi <= i1
            medvec(m) = v0 + (v1-v0)*(mi-i0)/(i1-i0);
        else
            medvec(m) = v1;                       % after end: stable
        end
    end
end

function idx = mabs(mm)
    % Resolve a month label to an absolute index with Jan-27 = 1.
    % Accepts 'Mon-YY', bare 'Mon' (defaults to yr 27), or 'MonYY' (e.g. Nov26).
    s = string(mm);
    mon = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
    if contains(s,'-')
        parts = split(s,'-');
        mi = find(mon==parts(1)); yy = str2double(parts(2));
    else
        pre = extractBefore(s,4);            % first 3 letters = month
        mi  = find(mon==pre);
        yytxt = extractAfter(s,3);           % trailing digits = year, if any
        if strlength(yytxt) >= 1
            yy = str2double(yytxt);
        else
            yy = 27;                         % bare month defaults to 2027
        end
    end
    idx = mi + 12*(yy-27);
end

function d = decay_anchored(month_age, dApr, dDec)
    % Absolute efficiency drop (Δη) vs storage age in months since cell
    % production (Nov-26 = age 0). Δη = 0 at Nov-26, ramping linearly to the
    % April anchor (kink at Apr), then extrapolated to the December anchor.
    %   Nov-26 = 0, Dec-26 = 1, Jan-27 = 2, Apr-27 = 5, Dec-27 = 13.
    aprA = 5; decA = 13;                    % Apr, Dec as age-in-months (from Nov-26)
    d = zeros(size(month_age));
    for k = 1:numel(month_age)
        a = month_age(k);
        if a <= aprA
            d(k) = dApr * a/aprA;           % 0 at Nov-26 -> dApr at Apr
        else
            d(k) = dApr + (dDec-dApr)*(a-aprA)/(decA-aprA);  % Apr -> Dec
        end
    end
end
