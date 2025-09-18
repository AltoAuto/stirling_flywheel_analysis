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
    spec.RPM = 600;          % Ω placeholder 
    spec.Cs  = 0.02;         % Cf allowable speed ripple (fraction)
    const.f     = spec.RPM/60;         % [Hz]
    const.omega = 2*pi*const.f;        % [rad/s]
    
    % Assumptions (kept explicit for the report)
    assume.idealGas     = true;
    assume.isothermal   = true;
    assume.perfectRegen = true;
    assume.uniformP     = true;
    assume.noFriction   = true;
    
    % Temps [K]
    T.Th = 700;
    T.Tc = 300;
    T.Tr = (T.Th+T.Tc)/2;     % regenerator temp

    % Gas (air baseline)
    gas.R = 287;      % J/(kg·K)
    gas.m = NaN;   % kg (This is only temporary, will eventually calculated pV=mRT)
    
    % Geometry (beta-type baseline)
    geom.r    = 0.010;        % crank radius [m]
    geom.L    = 0.060;        % conrod length [m]
    geom.A_p  = 1.20e-4;      % power piston area [m^2]
    geom.A_d  = 1.20e-4;      % displacer proj. area [m^2]
    geom.Vh0  = 1.50e-5;      % hot clearance [m^3]
    geom.Vc0  = 1.50e-5;      % cold clearance [m^3]
    geom.Vr   = 1.00e-5;      % regenerator vol [m^3]
    geom.phi  = deg2rad(90);  % displacer leads power piston [rad]
    
    assert(geom.L >= 2*geom.r, 'Geometry domain error: need L >= 2r.');

    %% Angle Grid
    theta = linspace(0,2*pi,2001).';  % one cycle, high resolution
    deg   = theta*180/pi;             % for plotting labe   

    %% Run Banner
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
    
    % Appendix B gives Pmin at BDC (max total volume). Plug real value later:
    spec.Pmin = 1.20e5;   % [Pa] placeholder; replace with Table B.1 value
    
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
    plot_pv_cycle(v_spec, p, gas, T, Vtot, OUTDIR); % saves FigA_PV.png
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
    fprintf('Geom: r=%.3f, L=%.3f, A_p=%.2e, A_d=%.2e, Vh0=%.2e, Vc0=%.2e, Vr=%.2e\n', ...
       geom.r, geom.L, geom.A_p, geom.A_d, geom.Vh0, geom.Vc0, geom.Vr);
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
    xd = g.r * cos(theta + g.phi);  % signed displacement
    Vs = (g.A_d * xd)/2;                 % signed shuttle volume
    Vh = g.Vh0 + Vs;        % add to hot when Vs>0
    Vc = g.Vc0 - Vs;       % add to cold when Vs<0

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

function plot_pv_cycle(v_spec, p, gas, T, Vtot, OUTDIR)
% Required Fig A: plot PV loop of the engine (p vs specific volume),
% and overlay an ideal Stirling rectangle for visual reference.

% Ideal overlay construction (simple & standard):
% - Use Vmin = min(Vtot), Vmax = max(Vtot).
% - Isothermal at Th from Vmin→Vmax, and at Tc from Vmax→VmNaNin.
% - Isochores at Vmin and Vmax connect the two isotherms (pressure ratio ~ Th/Tc).

    % ENGINE LOOP
    f = figure('Color','w'); hold on; grid on; box on;
    plot(v_spec, p/1e3, 'k-', 'LineWidth', 1.8, 'DisplayName','Engine loop');
    
    % IDEAL STIRLING OVERLAY (very simple construction)
    Vmin = min(Vtot); Vmax = max(Vtot);
    mR   = gas.m * gas.R;
    
    % Isothermal at Th (expansion) and Tc (compression)
    vmin = Vmin/gas.m; vmax = Vmax/gas.m;
    vTh  = linspace(vmin, vmax, 200);
    vTc  = linspace(vmax, vmin, 200);
    pTh  = (mR*T.Th) ./ vTh;
    pTc  = (mR*T.Tc) ./ vTc;
    
    % Isochores: vertical lines at vmin and vmax linking pressures
    pIso_min = [ (mR*T.Tc)/vmin, (mR*T.Th)/vmin ];
    pIso_max = [ (mR*T.Th)/vmax, (mR*T.Tc)/vmax ];
    
    plot(vTh, pTh/1e3,  '--', 'LineWidth',1.4, 'DisplayName','Ideal isotherm (Th)');
    plot(vTc, pTc/1e3,  '--', 'LineWidth',1.4, 'DisplayName','Ideal isotherm (Tc)');
    plot([vmin vmin], pIso_min/1e3, ':', 'LineWidth',1.2, 'DisplayName','Isochore');
    plot([vmax vmax], pIso_max/1e3, ':', 'LineWidth',1.2, 'DisplayName','Isochore');
    
    xlabel('Specific volume v = V_{tot}/m (m^3/kg)','Interpreter','tex');
    ylabel('Pressure p (kPa)','Interpreter','tex');
    
    title('Fig A — PV Diagram (Engine vs. Ideal Stirling)');
    legend('Location','best');
    
    fn = fullfile(OUTDIR, 'FigA_PV.png');
    saveas(f, fn);
    fprintf('Saved as FigA_PV.png');
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
