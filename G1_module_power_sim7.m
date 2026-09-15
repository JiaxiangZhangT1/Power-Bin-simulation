%% G1 Module Power Bin Simulation vs Target
% JX / T1 Energy
% MP Solar (procured) + G2 cells, split by monthly MW plan (separate populations).
% MP cells age from Nov-26 (storage decay). Best & worst case scenarios.
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
CTM     = 0.956;                      % baseline cell-to-module ratio
SIGMA_MOD = 1.5;                      % module-level power std (W): real cell->module, From Shundi's data
PREMIUM_BOM = false;                  % (unused; BOM cases removed. Re-add rows to CASES to enable)
BOM_GAIN_W  = 3;                      % +3W from thicker wire / high-Tx glass (if BOM re-enabled)

%% ------------------------------------------------------------------
%  2. Efficiency distributions (discrete PMFs)
%  ------------------------------------------------------------------
% MP Solar procured cells (fixed distribution)
mp_eff  = [0.253 0.254 0.255 0.256 0.257 0.258];
mp_pmf  = [0.0063 0.0276 0.2383 0.4526 0.1626 0.0988];
mp_pmf  = mp_pmf/sum(mp_pmf);

% Screened MP distribution: procurement accepts only cells >= 25.6%.
% Keep the >=0.256 bins and renormalize (drops ~28% of cells, lifts mean).
MP_SCREEN = 0.256;
mp_pmf_hi = mp_pmf .* (mp_eff >= MP_SCREEN - 1e-9);
mp_pmf_hi = mp_pmf_hi / sum(mp_pmf_hi);

% G2 cells: same distribution SHAPE as MP, but shifted so the MEDIAN
% tracks the monthly ramp target. We shift the discrete support by the
% delta between the target median and the MP baseline median.
mp_median = weighted_median(mp_eff, mp_pmf);   % ~0.256

%% ------------------------------------------------------------------
%  3. Monthly plan and G2 availability (MW)
%  ------------------------------------------------------------------
months     = {'Jan','Feb','Mar','Apr','May','Jun', ...
              'Jul','Aug','Sep','Oct','Nov','Dec'};
g1_plan_MW = [462 414 478 462 478 462 478 478 462 478 446 466];  % Dec est. from trend
g2_MW      = [  0   0   0 139 204 170 172 175 178 225 182 195];  % Dec est. from trend
nM = numel(months);

% Month index relative to Jan-27 (1..11); storage clock starts Nov-26
month_num  = 1:nM;

%% ------------------------------------------------------------------
%  4. G2 efficiency ramp (median), best & worst case
%  ------------------------------------------------------------------
% Best:  25.3% Apr -> 25.6% Sep, then stable
% Worst: 25.2% Apr -> 25.5% Dec (Dec is past Nov window; sets the ramp slope)
g2_med_best  = ramp_median(months, {'Apr',0.253}, {'Sep',0.256});
g2_med_worst = ramp_median(months, {'Apr',0.252}, {'Dec',0.255});

%% ------------------------------------------------------------------
%  5. Storage decay on MP cells (aging from Nov-26)
%  ------------------------------------------------------------------
% Decay is an ABSOLUTE cell-efficiency drop (Δη), not a power multiplier.
% e.g. Δη = 0.001 means 25.4% -> 25.3% (≈2.4W/module), Δη = 0.002 ≈ 4.8W.
% Best:  Δη 0.001 @3mo, 0.002 @6mo, then STABLE
% Worst: Δη 0.001 @3mo, 0.002 @6mo, then LINEAR (0.002 per additional 6mo)
% MP cells age from Nov-26. A cell consumed in month m (Jan-27 = 1) has aged
% (m + 2) months, since Nov-26 -> Jan-27 is 2 months.
storage_age = month_num + 2;          % Jan-27 = 3 months aged (from Nov-26)
dEff_best  = arrayfun(@(a) storage_dEff(a,'best'),  storage_age);  % abs eff drop
dEff_worst = arrayfun(@(a) storage_dEff(a,'worst'), storage_age);

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

% case definitions: name, ramp source, BOM on/off, MP screen >=25.6%
CASES = { ...
    'best',      'best',  false, false; ...
    'worst',     'worst', false, false; ...
    'best_hi',   'best',  false, true ; ...   % MP screened >=25.6%
    'worst_hi',  'worst', false, true };
nCase = size(CASES,1);
results = struct();

% floor bins
binLbl   = [600 605 610 615 620];      % bin labels (W)
binEdges = [binLbl 625];               % floor edges: [600,605),...,[620,625)


for cidx = 1:nCase
    scen    = CASES{cidx,1};
    basecase= CASES{cidx,2};
    useBOM  = CASES{cidx,3};
    screenMP= CASES{cidx,4};
    if strcmp(basecase,'best')
        g2_med = g2_med_best;  mp_dEff = dEff_best;
    else
        g2_med = g2_med_worst; mp_dEff = dEff_worst;
    end
    bom_add = useBOM * BOM_GAIN_W;
    if screenMP, mpPMF = mp_pmf_hi; else, mpPMF = mp_pmf; end   % MP dist for this case

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
    results.(scen).bom  = useBOM;
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

% (a) Power distribution vs target — RIDGELINE in a single axes.
% Curves are analytic Gaussians (no Monte Carlo): each month's MP and G2
% module distributions are normal, weighted by MW share and stacked.
% MP (blue) and G2 (orange) are separate populations; MP is independent of G2.
xg = 600:0.25:625;                           % smooth power grid for PDFs
caseList = fieldnames(results);
ROW = 0.9;                                   % vertical spacing between months
SCL = 2.0;                                   % height scale for each ridge
cMP = [0 0 0]; cG2 = [0.871 0.996 0.424];   % MP=black (procured), G2=#DEFE6C green
for ci = 1:numel(caseList)
    scen = caseList{ci};
    figure('Name',['G1 power ridgeline - ' scen],'Color','w', ...
           'Position',[100+40*ci 60 720 900]);
    ax = axes; hold(ax,'on'); set(ax,'Color','w');   % white plot background
    for m = nM:-1:1                          % draw bottom(Nov) first so top overlaps
        base = (nM - m) * ROW;               % Jan highest, Nov lowest
        wg = results.(scen).wG2(m);
        % analytic densities (per 1W, to match earlier probability scaling),
        % weighted by population share
        dMP = (1-wg) * normpdf(xg, results.(scen).muMP(m), results.(scen).sdMP(m));
        dG2 =    wg  * normpdf(xg, results.(scen).muG2(m), results.(scen).sdG2(m));
        yMP  = base + SCL*dMP;               % blue on baseline
        yTOP = base + SCL*(dMP + dG2);       % orange stacked above blue
        fill([xg fliplr(xg)], [yTOP fliplr(yMP)], cG2, 'EdgeColor','k','LineWidth',0.5);
        fill([xg fliplr(xg)], [yMP fliplr(base*ones(size(xg)))], cMP, 'EdgeColor','none');
        plot([xg(1) xg(end)], [base base], 'Color',[0.6 0.6 0.6], 'LineWidth',0.3);
        % month / avg-power label at left
        text(598.5, base+0.12, months{m}, ...
             'HorizontalAlignment','right','FontSize',9,'VerticalAlignment','bottom', ...
             'Color','k');
    end
    xline(615,'--','Color',[0.5 0.5 0.5]);
    xlim([600 625]); ylim([-0.1 (nM-1)*ROW + SCL*0.4 + 0.2]);
    set(ax,'YTick',[],'XColor','k','YColor','none','FontSize',9);
    xlabel('Module power (W)','FontSize',10,'Color','k');
    box off;
    title(sprintf('G1 Module Power Bin by Cell Source vs Month — %s (\\sigma=%.1fW)', ...
          caseLabel(scen),SIGMA_MOD),'FontSize',11,'Color','k');
    % manual legend swatches
    yl = ylim; xL = 621.5; yb = yl(2)-0.35;
    fill(xL+[0 1.2 1.2 0], yb+[0 0 0.12 0.12], cMP,'EdgeColor','none');
    text(xL+1.5, yb+0.06,'MP (procured)','FontSize',8,'VerticalAlignment','middle','Color','k');
    fill(xL+[0 1.2 1.2 0], yb-0.18+[0 0 0.12 0.12], cG2,'EdgeColor','k','LineWidth',0.5);
    text(xL+1.5, yb-0.12,'G2','FontSize',8,'VerticalAlignment','middle','Color','k');
end

% (b) Floor-binned (5W bins) vs target — RIDGELINE, same style as (a).
% Each month is a stacked step-area over the 5W bins [600 605 610 615 620],
% MP (blue) + G2 (orange), offset vertically. One croppable graphic per case.
% 5W bins concentrate mass into single tall bars, so use a smaller vertical
% scale than (a) to keep the 615 bar within its row.
binX = binLbl;                               % 600 605 610 615 620
SCL_B = 0.75;                                % height scale for 5W-bin bars
for ci = 1:numel(caseList)
    scen = caseList{ci};
    figure('Name',['G1 5W-bin ridgeline - ' scen],'Color','w', ...
           'Position',[780+40*ci 60 720 900]);
    ax = axes; hold(ax,'on'); set(ax,'Color','w');
    for m = nM:-1:1
        base = (nM - m) * ROW;
        dmp = results.(scen).distMP(m,:);    % MP fraction per bin
        dg2 = results.(scen).distG2(m,:);    % G2 fraction per bin
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
        text(596.5, base+0.12, months{m}, ...
             'HorizontalAlignment','right','FontSize',9,'VerticalAlignment','bottom', ...
             'Color','k');
    end
    xlim([596 623]); ylim([-0.1 (nM-1)*ROW + SCL_B*1.0 + 0.2]);
    set(ax,'XTick',binLbl,'XTickLabel',["600" "605" "610" "615" "620"], ...
        'YTick',[],'XColor','k','YColor','none','FontSize',9, ...
        'TickLabelInterpreter','tex');
    xlabel('Power bin (W)','FontSize',10,'Color','k'); box off;
    title(sprintf('G1 Module Power Bin (5W) by Cell Source vs Month — %s', ...
          caseLabel(scen)),'FontSize',11,'Color','k');
    % legend swatches
    yl = ylim; xL = 619; yb2 = yl(2)-0.35;
    fill(xL+[0 1.0 1.0 0], yb2+[0 0 0.12 0.12], cMP,'EdgeColor','none');
    text(xL+1.3, yb2+0.06,'MP','FontSize',8,'VerticalAlignment','middle','Color','k');
    fill(xL+[0 1.0 1.0 0], yb2-0.18+[0 0 0.12 0.12], cG2,'EdgeColor','k','LineWidth',0.5);
    text(xL+1.3, yb2-0.12,'G2','FontSize',8,'VerticalAlignment','middle','Color','k');
end

%% ------------------------------------------------------------------
%  8c. Input-assumption plots (worst case): G2 ramp & MP storage decay
%  ------------------------------------------------------------------
xm = 1:nM;

% (c) Worst-case G2 efficiency ramp (G2 production starts Apr)
figure('Name','Worst-case G2 efficiency ramp','Color','w', ...
       'Position',[300 120 640 380]);
ax = axes; hold(ax,'on'); set(ax,'Color','w');
g2start = 4;                                  % Apr = first month with G2
idx = g2start:nM;
cG2dark = [0.42 0.55 0.12];                    % darker shade of G2 green for line-on-white
plot(xm(idx), 100*g2_med_worst(idx), '-o','Color',cG2dark,'MarkerFaceColor',cG2dark, ...
     'LineWidth',1.6,'MarkerSize',5);
% mark the ramp anchors (Apr 25.2%, then climbing toward Dec 25.5%)
yline(25.2,'--','Color',[0.6 0.6 0.6]);
xlim([g2start-0.5 nM+0.5]); ylim([25.1 25.6]);
set(ax,'XTick',idx,'XTickLabel',months(idx),'XColor','k','YColor','k','FontSize',9);
xlabel('Month (2027)','FontSize',10,'Color','k');
ylabel('G2 median cell efficiency (%)','FontSize',10,'Color','k');
title('Worst-Case G2 Efficiency Ramp (25.2% Apr \rightarrow 25.5% Dec)', ...
      'FontSize',11,'Color','k');
grid on; box off;

% (d) Worst-case MP cell storage decay -> cell efficiency drop (% absolute)
figure('Name','Worst-case MP storage decay','Color','w', ...
       'Position',[360 120 640 380]);
ax = axes; hold(ax,'on'); set(ax,'Color','w');
% Absolute cell-efficiency change (percentage points, negative = drop):
dEta = -100 * dEff_worst;                    % e.g. dEff 0.002 -> -0.20 %pt
plot(xm, dEta, '-o','Color',cMP,'MarkerFaceColor',cMP, ...
     'LineWidth',1.6,'MarkerSize',5);
yline(0,'-','Color',[0.7 0.7 0.7]);
% reference: 6-month level (dEff = 0.002 = -0.20 %pt)
yline(-0.20,'--','Color',[0.6 0.6 0.6]);
text(nM, -0.20, '  6-mo level (-0.20 %pt)','Color',[0.4 0.4 0.4],'FontSize',8, ...
     'VerticalAlignment','top','HorizontalAlignment','right');
xlim([0.5 nM+0.5]); ylim([min(dEta)*1.15 0.05]);
set(ax,'XTick',xm,'XTickLabel',months,'XColor','k','YColor','k','FontSize',9);
xlabel('Month (2027)','FontSize',10,'Color','k');
ylabel('MP cell efficiency drop (% absolute)','FontSize',10,'Color','k');
title('Worst-Case MP Cell Efficiency Drop from Storage (aging from Nov-26)', ...
      'FontSize',11,'Color','k');
grid on; box off;

%% ------------------------------------------------------------------
%  9. Summary table — average module power, all cases
%  ------------------------------------------------------------------
fprintf('\n=== Monthly avg module power (W), sigma_mod=%.1fW ===\n',SIGMA_MOD);
fprintf('MP screen >=%.1f%% keeps %.0f%% of MP cells.\n', ...
        100*MP_SCREEN, 100*sum(mp_pmf(mp_eff>=MP_SCREEN-1e-9)));
fprintf('%-6s %8s %8s %9s %10s %7s\n','Month', ...
        'best','worst','best_hi','worst_hi','G2sh');
for m=1:nM
    fprintf('%-6s %8.1f %8.1f %9.1f %10.1f %6.0f%%\n', months{m}, ...
        results.best.avg(m), results.worst.avg(m), ...
        results.best_hi.avg(m), results.worst_hi.avg(m), ...
        100*g2_MW(m)/g1_plan_MW(m));
end
fprintf('\nTarget avg ~615W (bins: 16%% @610, 65%% @615, 19%% @620; 605 target 0%%).\n');
fprintf('_hi cases: MP procurement screened to >=25.6%% cells (storage decay still applies).\n');

%% ================= helper functions =================
function s = caseLabel(scen)
    switch scen
        case 'best',     s = 'BEST case';
        case 'worst',    s = 'WORST case';
        case 'best_hi',  s = 'BEST case, MP screened \geq25.6%';
        case 'worst_hi', s = 'WORST case, MP screened \geq25.6%';
        otherwise,       s = upper(scen);
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
    % 'Mon-YY' or bare 'Mon' -> absolute month index with Jan(-27) = 1.
    % Bare month names default to year 27.
    parts = split(string(mm),'-');
    mon = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
    mi  = find(mon==parts(1));
    if numel(parts) >= 2
        yy = str2double(parts(2));
    else
        yy = 27;                 % default plan year
    end
    idx = mi + 12*(yy-27);
end

function dEff = storage_dEff(ageMonths, scen)
    % returns ABSOLUTE cell-efficiency drop (Δη), e.g. 0.001 = 25.4->25.3%
    % Best:  0.001 @3mo, 0.0015 @6mo, 0.002 @12mo, then hold.
    % Worst: 0.001 @3mo, 0.002  @6mo, 0.003 @12mo, then continue.
    a = ageMonths;
    if strcmp(scen,'best')
        if a <= 3
            dEff = 0.001*a/3;                       % 0 -> 0.001 at 3mo
        elseif a <= 6
            dEff = 0.001 + 0.0005*(a-3)/3;          % 0.001 -> 0.0015 at 6mo
        elseif a <= 12
            dEff = 0.0015 + 0.0005*(a-6)/6;         % 0.0015 -> 0.002 at 12mo
        else
            dEff = 0.002;                           % stabilize
        end
    else % worst
        if a <= 3
            dEff = 0.001*a/3;                       % 0 -> 0.001 at 3mo
        elseif a <= 6
            dEff = 0.001 + 0.001*(a-3)/3;           % 0.001 -> 0.002 at 6mo
        elseif a <= 12
            dEff = 0.002 + 0.001*(a-6)/6;           % 0.002 -> 0.003 at 12mo
        else
            dEff = 0.003 + 0.001*(a-12)/6;          % keep decaying past 12mo
        end
    end
end
