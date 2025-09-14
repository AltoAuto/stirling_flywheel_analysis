function main
% Stirling Engine Analysis & Flywheel Design — (Scaffold + Volumes)
% One-file submission; all helpers are local functions below.
% After this step: volumes are validated and a diagnostic figure is saved.

    % HOUSE KEEPING
    format compact
    rng(0);                              
    thisFile = mfilename('fullpath'); 
    thisDir  = fileparts(thisFile);
    OUTDIR = fullfile(thisDir, 'outputs');
    if ~exist(OUTDIR, 'dir')
        mkdir(OUTDIR);
    end
    
    % Scope & Guardrails  %
    spec = spec_from_tableB_placeholder(); % TODO: replace with real Table B.1 later
    const.f     = spec.RPM/60;
    const.omega = 2*pi*const.f;
    
    % Assumptions (boolean flags kept for clarity; not all used yet)
    assume.idealGas     = true;
    assume.isothermal   = true;
    assume.perfectRegen = true;
    assume.uniformP     = true;
    assume.noFriction   = true;
    
    % Temperatures (K) — baseline, adjust as needed
    T.Th = 700;      % hot space
    T.Tc = 300;      % cold space
    T.Tr = 500;      % regenerator
    
    % Gas (air as baseline)
    gas.R = 287;     % J/(kg·K)
    gas.m = 1.0e-2;  % kg (tunable later once we set a target pressure)
    
    % Geometry (beta-type baseline; ensures L >= 2r for slider-crank domain)
    geom = defaultGeometry();
    geom.R  = gas.R;
    geom.m  = gas.m;
    geom.Vr = 1.00e-5;    % regenerator volume (m^3)
    
    assert(geom.L >= 2*geom.r, 'Geometry domain error: require L >= 2r for slider-crank.');
    
    % Angle grid (high resolution to avoid numerical noise downstream)
    theta = linspace(0,2*pi,2001).';   % Nx1 radians
    deg   = theta * 180/pi;            % for plotting convenience
    
    banner(spec, assume, T, gas, geom, OUTDIR);
    
    %---------------------------------------%
    % Phase 1 — Kinematics → Volumes(θ)     %
    %---------------------------------------%
    [Vh, Vc, Vtot, Vp, Vs] = volumes(theta, geom);
    
    % Basic hygiene checks
    assert(all(Vh  > 0), 'Negative hot volume detected. Increase Vh0 or adjust geometry.');
    assert(all(Vc  > 0), 'Negative cold volume detected. Increase Vc0 or adjust geometry.');
    assert(all(Vtot> 0), 'Negative total volume detected. Check clearances and piston geometry.');
    
    % Quick summary
    fprintf('Volumes (m^3): Vtot_min=%.3e, Vtot_max=%.3e, stroke ΔV=%.3e\n', ...
          min(Vtot), max(Vtot), max(Vtot)-min(Vtot));
    
    % Diagnostic plot (helps validate shapes before we do thermodynamics)
    plot_volumes(deg, Vh, Vc, Vp, Vtot, OUTDIR);
    
    % Next steps (in following messages): pressure, PV loop, torque/power, energy/J, speed ripple, phase sweep.
    end
    
    % ======================= Local functions =======================
    
    function spec = spec_from_tableB_placeholder()
    % Placeholder for course specs (Table B.1) — update later with the real values.
    spec.RPM = 600;      % target operating speed
    spec.Cs  = 0.02;     % allowable speed fluctuation (fractional), e.g., 2%
    % Add target power or other limits here when Table B.1 become avaliable
    end
    
    function geom = defaultGeometry()
    % Baseline beta-type geometry and clearances (SI units)
    geom.r    = 0.010;      % crank radius [m]
    geom.L    = 0.060;      % conrod length [m]  (>= 2r)
    geom.A_p  = 1.20e-4;    % power piston area [m^2]
    geom.A_d  = 1.20e-4;    % displacer projected area [m^2]
    geom.Vh0  = 1.50e-5;    % hot clearance volume [m^3]
    geom.Vc0  = 1.50e-5;    % cold clearance volume [m^3]
    geom.phi  = deg2rad(90);% displacer leads power piston [rad]
    end
    
    function [Vh, Vc, Vtot, Vp, Vs] = volumes(theta, g)
    % Compute hot/cold/total volumes for a beta-type Stirling with:
    % - Exact slider-crank for power piston displacement (x_p)
    % - Sinusoidal displacer (x_d) that shuttles volume between hot and cold
    %
    % Outputs are Nx1 arrays matching theta.
    
    theta = ensure_col(theta);
    
    % --- Power piston slider-crank position (exact) ---
    % Geometry in-line with standard slider-crank: x increases with cosine
    % Choose reference so minimum x_p = 0 (so Vp is purely additional volume)
    xp = g.r*cos(theta) + sqrt(max(g.L^2 - (g.r*sin(theta)).^2, 0));
    xp = xp - min(xp);          % shift so min(xp) = 0
    Vp = g.A_p * xp;            % power-piston contributed volume (always >= 0)
    
    % --- Displacer sinusoid (shuttle volume) ---
    % Positive Vs adds to hot; negative adds to cold.
    xd = g.r * cos(theta + g.phi);  % signed displacement
    Vs = g.A_d * xd;                % signed "shuttle" volume
    
    % Split shuttle between hot/cold spaces with clearances
    Vh = g.Vh0 + max(Vs, 0);
    Vc = g.Vc0 + max(-Vs,0);
    
    % Regenerator is constant volume; total is sum of compartments + piston
    Vtot = Vh + Vc + Vp + g.Vr;
    end
    
    function plot_volumes(deg, Vh, Vc, Vp, Vtot, OUTDIR)
    % Save a diagnostic figure to validate kinematics/volumes before thermo.
    f = figure('Color','w'); hold on; grid on;
    plot(deg, Vh,  'LineWidth',1.5, 'DisplayName','V_h(\theta)');
    plot(deg, Vc,  'LineWidth',1.5, 'DisplayName','V_c(\theta)');
    plot(deg, Vp,  'LineWidth',1.5, 'DisplayName','V_p(\theta)');
    plot(deg, Vtot,'k','LineWidth',2.0, 'DisplayName','V_{tot}(\theta)');
    xlabel('\theta (deg)'); ylabel('Volume (m^3)');
    title('Diagnostic: Volumes vs. Crank Angle');
    xlim([0 360]);
    legend('Location','best'); box on;
    fn = fullfile(OUTDIR, 'Fig0_Volumes.png');
    saveas(f, fn);
    fprintf('Saved diagnostic figure: %s\n', fn);
    end
    
    function x = ensure_col(x)
    if isrow(x), x = x.'; end
    end
    
    function banner(spec, assume, T, gas, geom, OUTDIR)
    fprintf('--- Stirling Engine Project (Step 1: Scaffold + Volumes) ---\n');
    fprintf('OUTDIR: %s\n', OUTDIR);
    fprintf('Specs: RPM=%g, Cs(target)=%.3f\n', spec.RPM, spec.Cs);
    fprintf('Temps [K]: Th=%g, Tc=%g, Tr=%g\n', T.Th, T.Tc, T.Tr);
    fprintf('Gas: R=%.1f J/(kg·K), m=%.3e kg\n', gas.R, gas.m);
    fprintf('Geom: r=%.3f m, L=%.3f m, A_p=%.2e m^2, A_d=%.2e m^2\n', ...
          geom.r, geom.L, geom.A_p, geom.A_d);
    fprintf('Clearances: Vh0=%.2e m^3, Vc0=%.2e m^3, Vr=%.2e m^3\n', ...
          geom.Vh0, geom.Vc0, geom.Vr);
    fprintf('Phase lead (phi): %.1f deg\n', geom.phi*180/pi);
    fprintf('Assumptions: idealGas=%d, isothermal=%d, perfectRegen=%d, uniformP=%d, noFriction=%d\n', ...
          assume.idealGas, assume.isothermal, assume.perfectRegen, assume.uniformP, assume.noFriction);
    fprintf('-----------------------------------------------------------\n');
    end
