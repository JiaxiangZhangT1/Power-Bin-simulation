%% MP Cell Efficiency Distribution — Shift-Higher Candidates
% JX / T1 Energy
% Generates candidate MP efficiency distributions shifted higher by a single
% KNOB (mix fraction toward a higher-shifted version), while keeping shape and
% spread similar. Three methods are compared:
%   1. Empirical shift   - blend real PMF with a one-bin-up version
%   2. Skew-normal blend - smooth left-skewed fit, blend with a shifted copy
%   3. Beta blend        - bounded-range fit, blend with a shifted copy
% Tune KNOB (0 = original, 1 = fully shifted) and re-run to compare.

clear; clc; close all;

%% ------------------------------------------------------------------
%  Inputs: measured MP cell efficiency distribution
%  ------------------------------------------------------------------
eff    = [0.253 0.254 0.255 0.256 0.257 0.258];   % bin centers
pmf    = [0.0063 0.0276 0.2383 0.4526 0.1626 0.0988];
pmf    = pmf/sum(pmf);
effpct = eff*100;                                  % in %

SHIFT_PT = 0.1;                 % how far the "fully shifted" version moves (%pt)
KNOBS    = [0 0.3 0.6 1.0];     % knob values to overlay (mix fraction)

%% ------------------------------------------------------------------
%  Generate and plot all three methods
%  ------------------------------------------------------------------
methods = {'Empirical shift','Skew-normal blend','Beta blend'};
cols = lines(numel(KNOBS));

figure('Color','w','Position',[150 40 1100 1000]);
nK = numel(KNOBS);
for ki = 1:nK                                  % row = knob value
    for mi = 1:3                               % col = method
        subplot(nK, 3, (ki-1)*3 + mi); hold on; set(gca,'Color','w');
        fr = KNOBS(ki);
        switch mi
            case 1, p = empirical_shift(pmf, fr);
            case 2, p = skew_blend(effpct, pmf, fr, SHIFT_PT);
            case 3, p = beta_blend(effpct, fr, SHIFT_PT);
        end
        bar(effpct, p, 'BarWidth', 0.75, 'FaceColor', cols(ki,:), ...
            'EdgeColor','k','LineWidth',0.3);
        xticks(effpct); ylim([0 0.5]);
        set(gca,'XColor','k','YColor','k','FontSize',8); grid on; box off;
        if ki==1, title(methods{mi},'Color','k','FontSize',10); end
        if mi==1, ylabel(sprintf('knob=%.1f\nProbability',fr),'Color','k','FontSize',8); end
        if ki==nK, xlabel('Cell efficiency (%)','Color','k','FontSize',8); end
    end
end
sgtitle('MP Cell Efficiency — Shift-Higher Candidates (rows: knob, cols: method)','Color','k');

%% ------------------------------------------------------------------
%  Console: mean vs knob for each method
%  ------------------------------------------------------------------
fprintf('Mean cell efficiency (%%) vs knob:\n');
fprintf('%-6s %12s %12s %12s\n','knob','empirical','skew-normal','beta');
for ki = 1:numel(KNOBS)
    fr = KNOBS(ki);
    em = sum(effpct.*empirical_shift(pmf,fr));
    sk = sum(effpct.*skew_blend(effpct,pmf,fr,SHIFT_PT));
    bt = sum(effpct.*beta_blend(effpct,fr,SHIFT_PT));
    fprintf('%-6.1f %12.3f %12.3f %12.3f\n', fr, em, sk, bt);
end

%% ================= distribution methods ==========================
function p = empirical_shift(pmf, frac)
    % Blend real PMF with a version where mass moves up one bin.
    up = zeros(size(pmf));
    up(2:end) = pmf(1:end-1);       % move each bin's mass up one bin
    up(end)   = up(end) + pmf(end); % top bin accumulates
    p = (1-frac)*pmf + frac*up;
    p = p/sum(p);
end

function p = skewnorm_pmf(effpct, a_sk, loc, scale, locShift)
    edges = [effpct-0.05, effpct(end)+0.05];       % bin edges around centers
    c = skewcdf(edges, a_sk, loc+locShift, scale);
    p = diff(c); p = max(p, 0); p = p/sum(p);      % clamp tiny negatives
end

function p = skew_blend(effpct, pmf, frac, shiftPt)
    % Fit a left-skewed normal to data moments; blend with shifted copy.
    mu  = sum(effpct.*pmf);
    v   = sum((effpct-mu).^2 .* pmf);
    a_sk = -3;                                      % left skew (long low tail)
    d = a_sk/sqrt(1+a_sk^2);
    scale = sqrt(v/(1 - 2*d^2/pi));
    loc = mu - scale*d*sqrt(2/pi);
    base = skewnorm_pmf(effpct, a_sk, loc, scale, 0);
    hi   = skewnorm_pmf(effpct, a_sk, loc, scale, shiftPt);
    p = (1-frac)*base + frac*hi; p = p/sum(p);
end

function p = beta_pmf_local(effpct, aB, bB, loB, hiB, shiftPt)
    edges = [effpct-0.05, effpct(end)+0.05];
    z = (edges - loB)/(hiB - loB);
    z = min(max(z + shiftPt/(hiB-loB), 0), 1);
    c = betacdf_local(z, aB, bB);
    p = diff(c); if sum(p)>0, p = p/sum(p); end
end

function p = beta_blend(effpct, frac, shiftPt)
    aB = 5; bB = 2.5; loB = 25.2; hiB = 25.9;       % left-skewed, peak ~25.6
    base = beta_pmf_local(effpct, aB, bB, loB, hiB, 0);
    hi   = beta_pmf_local(effpct, aB, bB, loB, hiB, shiftPt);
    p = (1-frac)*base + frac*hi; if sum(p)>0, p = p/sum(p); end
end

%% ================= stat helpers (no toolbox needed) ==============
function c = skewcdf(x, a, loc, scale)
    % Skew-normal CDF: Phi(z) - 2*T(z,a), z=(x-loc)/scale
    z = (x - loc)/scale;
    c = 0.5*(1+erf(z/sqrt(2))) - 2*owenT(z, a);
end

function T = owenT(h, a)
    % Owen's T function by numerical integration (vectorized over h).
    T = zeros(size(h));
    ng = 60; t = linspace(0,1,ng);
    for i = 1:numel(h)
        f = exp(-0.5*h(i)^2*(1+t.^2))./(1+t.^2);
        T(i) = a/(2*pi) * trapz(t, f);
    end
end

function c = betacdf_local(z, a, b)
    % Regularized incomplete beta via numerical integration (vectorized).
    c = zeros(size(z));
    ng = 200; v = linspace(0,1,ng);
    den = trapz(v, v.^(a-1).*(1-v).^(b-1));
    for i = 1:numel(z)
        if z(i)<=0, c(i)=0;
        elseif z(i)>=1, c(i)=1;
        else
            u = linspace(0, z(i), ng);
            c(i) = trapz(u, u.^(a-1).*(1-u).^(b-1))/den;
        end
    end
end
