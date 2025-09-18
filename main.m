function main
% Stirling Engine Analysis & Flywheel Design — (Scaffold + Volumes)
% One-file submission; all helpers are local functions below.
% After this step: volumes are validated and a diagnostic figure is saved.

    %% HOUSE KEEPING
    format compact
    rng(0);                              
    thisFile = mfilename('fullpath'); 
    thisDir  = fileparts(thisFile);
    OUTDIR = fullfile(thisDir, 'outputs');
    if ~exist(OUTDIR, 'dir')
        mkdir(OUTDIR);
    end
    fprintf('Output Folder Read: %s\n',OUTDIR);

    %% CONFIG/SPECS (Value are GPT generated)
    spec.RPM  = 650;          % Ω (rpm)            
    spec.Cs   = 0.003;        % Cf (fraction)     
    spec.Pmin = 5.00e5;       % Pa at BDC          
    spec.CR   = 1.20;         % Compression ratio  
    
    % Temperatures [K] (Table B.1)
    T.Th = 900;               % hot-side           << edit here
    T.Tc = 300;               % cold-side          << edit here
    T.Tr = 0.5*(T.Th + T.Tc); % regenerator (lumped assumption)
    
    % Derived kinematics
    const.f     = spec.RPM/60;      % Hz
    const.omega = 2*pi*const.f;     % rad/s
    
    % ---- Working fluid (Table B.1) ----
    fluid = "air";                    % "air" | "helium" | "hydrogen"  << edit here
    switch lower(fluid)
      case "air",      gas.R = 287;
      case "helium",   gas.R = 2077;
      case "hydrogen", gas.R = 4124;
      otherwise, error('Unknown fluid "%s" — add its R.', fluid);
    end
    gas.m = NaN;   % will be calibrated from Pmin at BDC (Step 3)
    
    % ---- Geometry (Table B.1) ----
    % Power piston (crank-slider)
    geom.r   = 0.025;          % power-piston crank length [m]        << edit here
    geom.L   = 0.075;          % power-piston conrod length [m]       << edit here
    geom.bore_p = 0.050;       % cylinder bore [m] (φ)                << edit here
    geom.A_p = pi*(geom.bore_p/2)^2;  % derived piston area [m^2]
    
    % Displacer (use volume-driven shuttle from Table B)
    Vdisp      = 3.0e-5;       % displacer volume [m^3]               << edit here
    geom.Vs_amp = Vdisp/2;     % sinusoid amplitude ⇒ zero-mean shuttle
    
    geom.phi = deg2rad(90);    % phase shift ψ [deg → rad]            << edit here
    geom.Vr  = 2.00e-5;        % regenerator dead volume [m^3]        << edit here
    
    % Guard
    assert(geom.L >= 2*geom.r, 'Geometry error: need L >= 2*r for slider-crank.');
    
    % ---- Angle grid (keep high resolution) ----
    theta = linspace(0,2*pi,2001).';   % rad
    deg   = theta*180/pi;
    
    % ---- Use CR to set clearances so model exactly hits CR ----
    % Build power-piston displacement/volume once
    xp = geom.r*cos(theta) + sqrt(max(geom.L^2 - (geom.r*sin(theta)).^2, 0));
    xp = xp - min(xp);
    Vp = geom.A_p * xp;
    Vp_max = max(Vp);
    
    % For beta with conservative displacer (shuttle only):
    % Vtot(θ) = Vconst + Vp(θ), where Vconst = Vh0 + Vc0 + Vr
    % CR = (Vconst + Vp_max)/Vconst  ⇒  Vconst = Vp_max / (CR - 1)
    Vconst = Vp_max / (spec.CR - 1);
    sum_clearances = Vconst - geom.Vr;   % = Vh0 + Vc0
    
    % Split clearances (tune alpha if spec suggests otherwise)
    alpha = 0.5;                 % 50/50 split hot:cold (change if needed)
    geom.Vh0 = alpha     * sum_clearances;
    geom.Vc0 = (1 - alpha) * sum_clearances;
    
    % Print small summary of derived quantities
    fprintf('Derived from CR: Vconst=%.3e m^3, Vp_max=%.3e m^3, Vh0=%.3e, Vc0=%.3e\n', ...
            Vconst, Vp_max, geom.Vh0, geom.Vc0);
    
    % Assumptions (kept explicit for your report)
    assume.idealGas     = true;
    assume.isothermal   = true;
    assume.perfectRegen = true;
    assume.uniformP     = true;
    assume.noFriction   = true;
    
    % Banner (optional)
    banner(spec, assume, T, gas, geom, OUTDIR);

    %% Volume (θ)

    [Vh, Vc, Vtot, Vp, Vs] = volumes(theta, geom);

    % Check Point
    assert(all(Vh  > 0), 'Negative hot volume detected. Increase Vh0 or adjust geometry.');
    assert(all(Vc  > 0), 'Negative cold volume detected. Increase Vc0 or adjust geometry.');
    assert(all(Vtot> 0), 'Negative total volume detected. Check clearances/geometry.');
    
    Vc_eff = Vc+ Vp;     % cold space includes the piston chamber
    Vh_eff = Vh;          % hot space unchanged
    Vtot   = Vh_eff + Vc_eff + geom.Vr;  % keep total consistent

    fprintf('Volumes (m^3): Vtot_min=%.3e, Vtot_max=%.3e, stroke ΔV=%.3e\n', ...
          min(Vtot), max(Vtot), max(Vtot)-min(Vtot));
    
    % Diagnostic plot
    plot_volumes(theta*180/pi, Vh, Vc, Vp, Vtot, OUTDIR);

    %% Masss Calibration 
    
  
    % Find BDC index: where total volume is maximum
    [~, iBDC] = max(Vtot);
    Vh_bdc = Vh_eff(iBDC);
    Vc_bdc = Vc_eff(iBDC);
    Vr_bdc = geom.Vr;
    
    % Back-solve m from Schmidt relation at BDC:
    % p = mR / (Vh/Th + Vc/Tc + Vr/Tr)  =>  m = p/R * (Vh/Th + Vc/Tc + Vr/Tr)
    gas.m = calibrate_mass_from_Pmin(spec.Pmin, Vh_bdc, Vc_bdc, Vr_bdc, T, gas.R);
    
    fprintf('Calibrated gas mass from Pmin: m = %.6e kg (BDC)\n', gas.m);

    %% Compute Pressure p(θ) Schmidt - style
    p = pressure_schmidt(Vh_eff, Vc_eff, geom.Vr, T, gas.m, gas.R);    % Pa

    
    %====================%
    % DEBUG CHECKS       %
    %====================%
    
    % 1) Make sure we are using consistent volumes
    Vc_eff = Vc + Vp;          % cold space includes piston chamber (again, ensure defined here)
    Vh_eff = Vh;
    Vtot   = Vh_eff + Vc_eff + geom.Vr;
    
    % 2) Sanity on sizes & finiteness
    v_spec = Vtot ./ gas.m;    % specific volume for PV x-axis
    assert(isvector(v_spec) && isvector(p) && numel(v_spec)==numel(p), 'PV arrays size mismatch.');
    assert(all(isfinite(v_spec)) && all(isfinite(p)), 'NaN/Inf in p or v_spec.');
    
    % 3) Check the BDC anchor really matches spec.Pmin (within tiny tolerance)
    [~, iBDC] = max(Vtot);
    p_bdc = p(iBDC);
    fprintf('DEBUG: Pmin spec=%.2f kPa, p(theta_BDC)=%.2f kPa (Δ=%.2f%%)\n', ...
      spec.Pmin/1e3, p_bdc/1e3, 100*(p_bdc-spec.Pmin)/spec.Pmin);
    
    % 4) Quick time-series plots (theta vs v, theta vs p) to ensure each is smooth
    debug_plot_series(theta, v_spec, p, OUTDIR);
    
    % 5) Plot engine PV only (no overlay), closing the loop explicitly with markers
    plot_pv_engine_only(v_spec, p, OUTDIR);


    % sanity checks
    assert(all(p > 0), 'Non-physical pressure detected. Check mass calibration or temps.');
    fprintf('p(θ):  min=%.2f kPa, max=%.2f kPa\n', min(p)/1e3, max(p)/1e3);

    [pmx, imx] = max(p);  [pmn, imn] = min(p);
    [vmx, jmx] = max(v_spec); [vmn, jmn] = min(v_spec);
    fprintf('DEBUG extremes:\n');
    fprintf('  p_min=%.2f kPa at θ=%.1f deg | p_max=%.2f kPa at θ=%.1f deg\n', ...
            pmn/1e3, (imn-1)/(numel(p)-1)*360, pmx/1e3, (imx-1)/(numel(p)-1)*360);
    fprintf('  v_min=%.3e at θ=%.1f deg | v_max=%.3e at θ=%.1f deg\n', ...
            vmn, (jmn-1)/(numel(v_spec)-1)*360, vmx, (jmx-1)/(numel(v_spec)-1)*360);

    % PV PLOT — Engine vs Ideal Stirling (Fig A)
    v_spec = Vtot ./ gas.m;                         % specific volume (m^3/kg)
    plot_pv_cycle(v_spec, p, gas, T, Vtot, OUTDIR,'physical'); % saves FigA_PV.png
end



%% ------------------ Local Funtion ---------------------

% Summary Info Banner
function banner(spec, assume, T, gas, geom, OUTDIR)
    fprintf('\n--- Stirling Project: Step 1 Scaffold ---\n');
    fprintf('Output Folder Created: %s\n\n', OUTDIR);
    fprintf('Specs: RPM=%g, Cs=%.3f | f=%.2f Hz, omega=%.2f rad/s\n', ...
       spec.RPM, spec.Cs, spec.RPM/60, 2*pi*(spec.RPM/60));
    fprintf('Temps [K]: Th=%g, Tc=%g, Tr=%g\n', T.Th, T.Tc, T.Tr);
    fprintf('Gas: R=%.1f J/(kg·K), m=%.3e kg\n', gas.R, gas.m);
    fprintf('Geom: r=%.3f, L=%.3f, A_p=%.2e, Vh0=%.2e, Vc0=%.2e, Vr=%.2e\n', ...
       geom.r, geom.L, geom.A_p, geom.Vh0, geom.Vc0, geom.Vr);
    fprintf('Phase lead phi = %.1f deg\n', geom.phi*180/pi);
    fprintf('Assumptions: idealGas=%d, isothermal=%d, perfectRegen=%d, uniformP=%d, noFriction=%d\n', ...
       assume.idealGas, assume.isothermal, assume.perfectRegen, assume.uniformP, assume.noFriction);
    fprintf('-----------------------------------------\n');
end

% Beta-type Stirling volumes from kinematics.
function [Vh, Vc, Vtot, Vp, Vs] = volumes(theta, g)
% Outputs: hot Vh(θ), cold Vc(θ), total Vtot(θ), power-piston Vp(θ), shuttle Vs(θ).

    theta = ensure_col(theta);
    
    % --- Power piston (slider-crank exact) ---
    % Reference so min(xp)=0 ⇒ Vp >= 0 and clearances remain in Vh0/Vc0
    xp = g.r*cos(theta) + sqrt(max(g.L^2 - (g.r*sin(theta)).^2, 0));
    xp = xp - min(xp);
    Vp = g.A_p * xp;

    % --- Displacer (sinusoid with lead φ) ---
    Vs = g.Vs_amp * cos(theta + g.phi);  % signed, zero mean
    Vh = g.Vh0 + Vs;                     % add to hot
    Vc = g.Vc0 - Vs;                     % subtract from cold

    % --- Total volume includes regenerator dead volume ---
    Vtot = Vh + Vc + Vp + g.Vr;
end

function plot_volumes(deg, Vh, Vc, Vp, Vtot, OUTDIR)
% Save a diagnostic figure to confirm kinematics/volumes are sensible.
    f = figure('Color','w'); hold on; grid on; box on;
    plot(deg, Vh,  'LineWidth',1.5, 'DisplayName','V_h(\theta)');
    plot(deg, Vc,  'LineWidth',1.5, 'DisplayName','V_c(\theta)');
    plot(deg, Vp,  'LineWidth',1.5, 'DisplayName','V_p(\theta)');
    plot(deg, Vtot,'k','LineWidth',2.0, 'DisplayName','V_{tot}(\theta)');
    xlabel('\theta (deg)'); ylabel('Volume (m^3)');
    title('Diagnostic: Volumes vs. Crank Angle');
    xlim([0 360]); legend('Location','best');
    saveas(f, fullfile(OUTDIR, 'Fig0_Volumes.png'));
    fprintf('Saved as Fig0_Volumes.png\n');
end

function x = ensure_col(x)
    if isrow(x), x = x.'; end
end


% Schmidt relation (Who helped me with this? ChatGPT did!)
function m = calibrate_mass_from_Pmin(Pmin, Vh_bdc, Vc_bdc, Vr, T, R)
  denom = Vh_bdc/T.Th + Vc_bdc/T.Tc + Vr/T.Tr;
  m = (Pmin / R) * denom;
end

function p = pressure_schmidt(Vh, Vc, Vr, T, m, R)
% Lumped uniform-pressure model:
  p = (m*R) ./ ( Vh./T.Th + Vc./T.Tc + Vr./T.Tr );
end

function plot_pv_cycle(v_spec, p, gas, T, Vtot, OUTDIR, mode)
% Fig A: PV diagram showing (1) engine loop and (2) ideal Stirling rectangle.
% mode: 'normalized' (default) scales ideal to engine p-range; 'physical' uses mRT/v.

  if nargin < 7 || isempty(mode), mode = 'normalized'; end

  % Ensure columns & close engine loop
  v_spec = v_spec(:); p = p(:);
  vv = [v_spec; v_spec(1)];
  pp = [p;      p(1)];

  f = figure('Color','w'); hold on; grid on; box on;
  plot(vv, pp/1e3, 'k-', 'LineWidth', 1.8, 'DisplayName','Engine loop');

  % Engine ranges
  p_min = min(p); p_max = max(p);
  Vmin  = min(Vtot); Vmax = max(Vtot);
  vmin  = Vmin/gas.m; vmax = Vmax/gas.m;

  switch lower(string(mode))
    case "physical"
      % True ideal isotherms (may be tiny if m is small and v is large)
      mR  = gas.m * gas.R;
      p_tl = (mR*T.Th)/vmin;  p_tr = (mR*T.Th)/vmax;
      p_bl = (mR*T.Tc)/vmin;  p_br = (mR*T.Tc)/vmax;

    otherwise  % "normalized" — scale ideal to engine pressure band
      % Put top isotherm at engine max pressure, bottom at (Tc/Th) times that.
      p_tl = p_max;                  p_tr = p_max;                    % Th line (flat)
      p_bl = p_max * (T.Tc/T.Th);    p_br = p_max * (T.Tc/T.Th);      % Tc line (flat)
  end

  % Draw the 4 sides explicitly (no extra vertical lines)
  plot([vmin vmax], [p_tl p_tr]/1e3, '--', 'LineWidth',1.3, 'DisplayName','Ideal (Th)');
  plot([vmax vmax], [p_tr p_br]/1e3,  ':', 'LineWidth',1.2, 'DisplayName','Ideal iso');
  plot([vmax vmin], [p_br p_bl]/1e3, '--', 'LineWidth',1.3, 'DisplayName','Ideal (Tc)');
  plot([vmin vmin], [p_bl p_tl]/1e3,  ':', 'LineWidth',1.2, 'DisplayName','Ideal iso');

  xlabel('Specific volume v = V_{tot}/m (m^3/kg)','Interpreter','tex');
  ylabel('Pressure p (kPa)','Interpreter','tex');
  title('Fig A — PV Diagram (Engine vs. Ideal Stirling)');
  legend('Location','best');

  % Nice axes padding
  ypad = 0.05*(p_max - p_min + eps);
  ylim([(p_min-ypad)/1e3, (p_max+ypad)/1e3]);

  saveas(f, fullfile(OUTDIR,'FigA_PV.png'));
  fprintf('Saved FigA_PV.png\n');
end



function debug_plot_series(theta, v_spec, p, OUTDIR)
    % v(theta)
    f1 = figure('Color','w'); grid on; box on;
    plot(theta*180/pi, v_spec, 'LineWidth', 1.4);
    xlabel('\theta (deg)'); ylabel('v = V_{tot}/m (m^3/kg)');
    title('DEBUG: Specific volume vs. crank angle');
    saveas(f1, fullfile(OUTDIR,'DBG_v_vs_theta.png'));
    
    % p(theta)
    f2 = figure('Color','w'); grid on; box on;
    plot(theta*180/pi, p/1e3, 'LineWidth', 1.4);
    xlabel('\theta (deg)'); ylabel('p (kPa)');
    title('DEBUG: Pressure vs. crank angle');
    saveas(f2, fullfile(OUTDIR,'DBG_p_vs_theta.png'));
    
    fprintf('Saved DBG_v_vs_theta.png and DBG_p_vs_theta.png\n');
end

function plot_pv_engine_only(v_spec, p, OUTDIR)
    % Plot engine loop only, closed explicitly, with sparse markers to see the path
    v_spec = v_spec(:); p = p(:);
    
    f = figure('Color','w'); hold on; grid on; box on;
    % close the loop explicitly by appending first point
    vv = [v_spec; v_spec(1)];
    pp = [p;      p(1)];
    
    plot(vv, pp/1e3, 'k-', 'LineWidth', 1.8, 'DisplayName','Engine loop');
    
    % put markers every ~100 points to visualize traversal order
    step = max(1, floor(numel(v_spec)/36));    % ~10° increments
    idx  = 1:step:numel(v_spec);
    plot(v_spec(idx), p(idx)/1e3, 'o', 'MarkerSize', 3, 'DisplayName','samples');
    
    xlabel('v = V_{tot}/m (m^3/kg)');
    ylabel('p (kPa)');
    title('DEBUG: PV Engine Loop (closed, with markers)');
    legend('Location','best');
    saveas(f, fullfile(OUTDIR,'DBG_PV_engine_only.png'));
    fprintf('Saved DBG_PV_engine_only.png\n');
end
