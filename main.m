% Stirling Engine Analysis & Flywheel Design 

function main
close all; clc;

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

    %% CONFIG/SPECS
    spec.RPM  = 650;          % Ω (rpm)            
    spec.Cf   = 0.003;        % Cf (fraction)     
    spec.Pmin = 5.00e5;       % Pa at BDC
    spec.CR   = 1.20;         % Compression ratio
    
    % Temperatures [K] (Table B.1)
    T.Th = 900;               % hot-side
    T.Tc = 300;               % cold-side
    T.Tr = 0.5*(T.Th + T.Tc); % regenerator (lumped assumption)
    
    % Derived kinematics
    const.f     = spec.RPM/60;      % Hz
    const.omega = 2*pi*const.f;     % rad/s
    
    % Working fluid
    gas.R = 287;   % Air
    gas.m = NaN;   % will be calibrated from Pmin at BDC
    
    %% Geometry (Table B.1)
    % Power piston (crank-slider)
    geom.r   = 0.025;          % power-piston crank length [m]
    geom.L   = 0.075;          % power-piston conrod length [m]
    geom.bore_p = 0.050;       % cylinder bore [m] (φ)
    geom.A_p = pi*(geom.bore_p/2)^2;  % derived piston area [m^2]
    
    % Displacer (use volume-driven shuttle from Table B)
    Vdisp      = 4.0e-5;       % displacer volume [m^3]
    geom.Vs_amp = Vdisp/2;     % sinusoid amplitude ⇒ zero-mean shuttle
    
    % Phase angle
    phase_angle = 90;           % phase angle [deg]
    geom.phi = deg2rad(phase_angle);    % phase shift ψ [deg → rad]
    geom.Vr  = 2.00e-5;        % regenerator dead volume [m^3]

    % Angle grid (keep high resolution)
    theta = linspace(0,2*pi,4000);   % rad
    deg   = theta*180/pi;
    
    % Use CR to set clearances so model exactly hits CR
    % Build power-piston displacement/volume once
    xp = geom.r*cos(theta) + sqrt(max(geom.L^2 - (geom.r*sin(theta)).^2, 0));   % slider crank position
    Vp_max = max(geom.A_p * (xp - min(xp))); 
    
    % For beta with conservative displacer (shuttle only):
    % Vtot(θ) = Vconst + Vp(θ), where Vconst = Vh0 + Vc0 + Vr
    % CR = (Vconst + Vp_max)/Vconst 
    Vconst = Vp_max / (spec.CR - 1);
    sum_clearances = Vconst - geom.Vr;   % = Vh0 + Vc0
    
    % Split clearances (alpha from spec)
    alpha = 0.5;                 % 50/50 split hot:cold (change if needed)
    geom.Vh0 = alpha     * sum_clearances;
    geom.Vc0 = (1 - alpha) * sum_clearances;


    % Assumptions
    assume.idealGas     = true;
    assume.isothermal   = true;
    assume.perfectRegen = true;
    assume.uniformP     = true;
    assume.noFriction   = true;
    
    % Banner
    banner(spec, assume, T, gas, geom, OUTDIR);
    % Print small summary of derived quantities
    fprintf('Volume (from CR): Vconst=%.3e m^3, Vp_max=%.3e m^3, Vh0=%.3e, Vc0=%.3e\n', ...
            Vconst, Vp_max, geom.Vh0, geom.Vc0);

    %% Volume (θ)
    [Vh, Vc, Vtot, Vp, Vs] = volumes(theta, geom);

    % Check Point
    assert(all(Vh  > 0), 'Negative hot volume detected. Increase Vh0 or adjust geometry.');
    assert(all(Vc  > 0), 'Negative cold volume detected. Increase Vc0 or adjust geometry.');
    assert(all(Vtot> 0), 'Negative total volume detected. Check clearances/geometry.');
    
    Vc_eff = Vc+ Vp ;     % cold space includes the piston chamber
    Vh_eff = Vh ;          % hot space + displacer
    Vtot   = Vh_eff + Vc_eff + geom.Vr;  % keep total consistent

    fprintf('Volumes (m^3): Vtot_min=%.3e, Vtot_max=%.3e, stroke ΔV=%.3e\n', ...
          min(Vtot), max(Vtot), max(Vtot)-min(Vtot));
    
    % Volume plot
    plot_volumes(theta*180/pi, Vh, Vc, Vp, Vtot, OUTDIR);

    %% Masss Calibration 
      
    % Find BDC index: where total volume is maximum
    [~, iBDC] = max(Vtot);
    Vh_bdc = Vh_eff(iBDC);
    Vc_bdc = Vc_eff(iBDC);
    Vr_bdc = geom.Vr;
    
    % Back-solve m from Schmidt relation at BDC:
    gas.m = calibrate_mass_from_Pmin(spec.Pmin, Vh_bdc, Vc_bdc, Vr_bdc, T, gas.R);
    
    fprintf('Calibrated gas mass from Pmin: m = %.6e kg (BDC)\n', gas.m);

    %% Compute Pressure p(θ) Schmidt - style
    p = pressure_schmidt(Vh_eff, Vc_eff, geom.Vr, T, gas.m, gas.R);    % Pa

    % specific volume
    v_spec = Vtot ./ gas.m;  
    assert(isvector(v_spec) && isvector(p) && numel(v_spec)==numel(p), 'PV arrays size mismatch.');
    assert(all(isfinite(v_spec)) && all(isfinite(p)), 'NaN/Inf in p or v_spec.');

    %% Pressure and Volume (FigA)
    plot_pv_cycle(v_spec, p, gas, T, Vtot, OUTDIR);

    %%  Torque and Power (Fig B)
    [Trq, Wcyc, P1, P2, Tmean] = torque_power(theta, p, Vtot, spec.RPM, const.omega, OUTDIR, true);
         
    % Specific engine work from your model (J/kg·cycle)
    W_engine_spec = Wcyc / gas.m;
    
    % v-range for ideal cycle reference (specific volume)
    vmin = min(v_spec); 
    vmax = max(v_spec);
    lnCRv = log(vmax/vmin);                      % "volumetric CR" on specific-basis
    
    % Ideal Stirling references (perfect regeneration, isothermal)
    W_ideal_spec = gas.R * (T.Th - T.Tc) * lnCRv;   % J/kg·cycle
    Qin_ideal_spec = gas.R * T.Th * lnCRv;          % J/kg·cycle
    
    % 1) “Cycle utilization” vs ideal Stirling work (how much of the ideal rectangle you realize)
    eta_util = W_engine_spec / max(W_ideal_spec, 1e-12);
    
    % 2) Indicated thermal efficiency vs the ideal heat input (with perfect regen)
    %    (This is your modeled W divided by the ideal Stirling Q_in.)
    eta_i = W_engine_spec / max(Qin_ideal_spec, 1e-12);
    
    % 3) Theoretical (Carnot/Stirling-ideal) thermal efficiency, for reference
    eta_ideal = 1 - T.Tc / T.Th;
    
    fprintf('\nEFFICIENCY CHECKS (per cycle)\n');
    fprintf('  W_engine_spec   = %.0f J/kg·cycle\n', W_engine_spec);
    fprintf('  W_ideal_spec    = %.0f J/kg·cycle\n', W_ideal_spec);
    fprintf('  Qin_ideal_spec  = %.0f J/kg·cycle\n', Qin_ideal_spec);
    fprintf('  eta_util        = %.3f  (W_engine / W_ideal)\n', eta_util);
    fprintf('  eta_i           = %.3f  (W_engine / Qin_ideal)\n', eta_i);
    fprintf('  eta_ideal_ref   = %.3f  (1 - Tc/Th)\n\n', eta_ideal);

    % Print summary for the report
    relErr = abs(P1 - P2) / max(abs(P2), 1e-12);
    fprintf('Torque/Power summary:\n');
    fprintf('  Work per cycle  Wcyc   = %.6g J/cycle\n', Wcyc);
    fprintf('  Mean torque     Tmean  = %.6g N·m\n', Tmean);
    fprintf('  Power (P1= W*f)       = %.6g W\n', P1);
    fprintf('  Power (P2= Tmean*ω)   = %.6g W\n', P2);
    fprintf('  |P1-P2|/P2            = %.3f %%\n\n', 100*relErr);

    %%  Energy deviation and required J (Fig C) 
    [Ed, Epp, Jreq] = energy_and_inertia(theta, Trq, Tmean, const.omega, spec.Cf, OUTDIR,true);
    
    [omega_th, Cf_sim] = simulateSpeed(theta, Trq, Tmean, Jreq, const.omega, OUTDIR);
    
    fprintf('Energy buffer & flywheel sizing summary:\n');
    fprintf('  Epp  = %.6g J\n', Epp);
    fprintf('  Jreq = %.6g kg·m^2\n', Jreq);
    fprintf('  Cf_sim = %.3f %%  (target %.3f %%)\n\n', 100*Cf_sim, 100*spec.Cf);


    %% Phase sweep (Fig D)
    phis_deg = 0:1:360;
    resultsD = phase_sweep(phis_deg, theta, geom, T, gas, const, spec, OUTDIR);
    
    [~, koptW] = max(resultsD.Wcyc);
    [~, koptJ] = min(resultsD.Jreq);
    fprintf('Phase sweep summary:\n');
    fprintf('  Max work at phi = %g deg: Wcyc = %.6g J/cycle\n', phis_deg(koptW), resultsD.Wcyc(koptW));
    fprintf('  Min Jreq at phi = %g deg: Jreq = %.6g kg·m^2\n', phis_deg(koptJ), resultsD.Jreq(koptJ));


    %% Flywheel geometry
    fly_in.rho      = 7850;   % steel
    fly_in.w        = 0.025;  % width [m]
    fly_in.t        = 0.050;  % rim thickness [m] (ri = ro - t)
    fly_in.yield    = 250e6;  % Pa yield
    fly_in.vtip_max = 120;    % m/s
    
    fly = sizeFlywheel_minD(Jreq, fly_in, const.omega);
    print_flywheel(fly);

end



%%%%%%%%%%%%%%%%%%%%% Local Funtion %%%%%%%%%%%%%%%%%%%%%

%% Summary Info Banner
function banner(spec, assume, T, gas, geom, OUTDIR)
    fprintf('\n--- Stirling Project: Step 1 Scaffold ---\n');
    fprintf('Output Folder Created: %s\n\n', OUTDIR);
    fprintf('Specs: RPM=%g, Cf=%.3f | f=%.2f Hz, omega=%.2f rad/s\n', ...
       spec.RPM, spec.Cf, spec.RPM/60, 2*pi*(spec.RPM/60));
    fprintf('Temps [K]: Th=%g, Tc=%g, Tr=%g\n', T.Th, T.Tc, T.Tr);
    fprintf('Gas: R=%.1f J/(kg·K), m=%.3e kg\n', gas.R, gas.m);
    fprintf('Geom: r=%.3f, L=%.3f, A_p=%.2e, Vh0=%.2e, Vc0=%.2e, Vr=%.2e\n', ...
       geom.r, geom.L, geom.A_p, geom.Vh0, geom.Vc0, geom.Vr);
    fprintf('Phase lead phi = %.1f deg\n', geom.phi*180/pi);
    fprintf('Assumptions: idealGas=%d, isothermal=%d, perfectRegen=%d, uniformP=%d, noFriction=%d\n', ...
       assume.idealGas, assume.isothermal, assume.perfectRegen, assume.uniformP, assume.noFriction);
    fprintf('-----------------------------------------\n');
end

%% Beta-type Stirling volumes from kinematics.
function [Vh, Vc, Vtot, Vp, Vs] = volumes(theta, g)
% Outputs: hot Vh(θ), cold Vc(θ), total Vtot(θ), power-piston Vp(θ), shuttle Vs(θ).

    theta = ensure_col(theta);
    
    % --- Power piston (slider-crank exact) ---
    % Reference so min(xp)=0 ⇒ Vp >= 0 and clearances remain in Vh0/Vc0
    xp = g.r*cos(theta) + sqrt(max(g.L^2 - (g.r*sin(theta)).^2, 0));
    xp = xp - min(xp);
    Vp = g.A_p * xp;

    % --- Displacer (sinusoid with lead φ) ---
    Vs = g.Vs_amp * cos(theta + g.phi);  % displacer volume (signed, positive-more in hot, negative-more on cold)
    Vh = g.Vh0 + Vs;                     % add to hot
    Vc = g.Vc0 - Vs;                     % subtract from cold

    % --- Total volume ---
    Vtot = Vh + Vc + Vp + g.Vr;
end
function plot_volumes(deg, Vh, Vc, Vp, Vtot, OUTDIR)
% Save a diagnostic figure to confirm kinematics/volumes are sensible.
    if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

    f = figure('Color','w','Visible','off'); 
    hold on; 
    plot(deg, Vh,  'LineWidth',1.5, 'DisplayName','V_h(\theta)');
    plot(deg, Vc,  'LineWidth',1.5, 'DisplayName','V_c(\theta)');
    plot(deg, Vp,  'LineWidth',1.5, 'DisplayName','V_p(\theta)');
    plot(deg, Vtot,'k','LineWidth',2.0, 'DisplayName','V_{tot}(\theta)');
    xlabel('\theta (deg)'); ylabel('Volume (m^3)');
    title('Volumes vs. Crank Angle');
    xlim([0 360]); legend('Location','best');
    grid on; box on;

    % Set figure size (inches)
    w = 5.5; h = 2.99;
    set(f,'Units','inches','Position',[1 1 w h]);
    set(f,'PaperUnits','inches','PaperSize',[w h], ...
           'PaperPosition',[0 0 w h],'PaperPositionMode','manual');

    % Save high-res PNG (600 dpi)
    fn_png = fullfile(OUTDIR,'Fig0_Volumes.png');
    print(f, fn_png, '-dpng','-r600');

    % Optional: also save a vector PDF
    fn_pdf = fullfile(OUTDIR,'Fig0_Volumes.pdf');
    print(f, fn_pdf, '-dpdf');
end


function x = ensure_col(x)
    if isrow(x), x = x.'; end
end


%% THE CORE WHERE MAGIC HAPPENS
function m = calibrate_mass_from_Pmin(Pmin, Vh_bdc, Vc_bdc, Vr, T, R)
  m = (Pmin / R) * (Vh_bdc/T.Th + Vc_bdc/T.Tc + Vr/T.Tr);
end

function p = pressure_schmidt(Vh, Vc, Vr, T, m, R)
% Lumped uniform-pressure model:
  p = (m*R) ./ ( Vh./T.Th + Vc./T.Tc + Vr./T.Tr );
end


%% Fig A: PV diagram showing (1) engine loop and (2) ideal Stirling rectangle.
function plot_pv_cycle(v_spec, p, gas, T, Vtot, OUTDIR)

  % Ensure columns & close engine loop
  v_spec = v_spec(:); p = p(:);
  vv = [v_spec; v_spec(1)];
  pp = [p;      p(1)];

  f = figure('Color','w'); hold on; grid off; box off;
  plot(vv, pp/1e3, 'k-', 'LineWidth', 1.8, 'DisplayName','Engine loop');

  % Engine ranges
  p_min = min(p); 
  p_max = max(p);
  vmin = min(v_spec);
  vmax = max(v_spec);

  % Build ideal Stirling isotherms (hyperbolas) between vmin and vmax
  Rspec = gas.R;

  % Base (physical) curves
  v_iso_inc = linspace(vmin, vmax, 400);       % increasing v
  v_iso_dec = linspace(vmax, vmin, 400);       % decreasing v
  p_hot_phys_inc  = (Rspec*T.Th) ./ v_iso_inc;    % Th, expand: vmin->vmax
  p_cold_phys_dec = (Rspec*T.Tc) ./ v_iso_dec;    % Tc, compress: vmax->vmin


  p_hot  = p_hot_phys_inc;
  p_cold = p_cold_phys_dec;

  % Corner pressures (shared exactly by segments to avoid gaps)
  p_hot_vmax  = p_hot(end);                    % Th at v = vmax
  p_cold_vmax = p_cold(1);                     % Tc at v = vmax
  p_cold_vmin = p_cold(end);                   % Tc at v = vmin
  p_hot_vmin  = p_hot(1);                      % Th at v = vmin

  % Plot hot isotherm (Th): vmin -> vmax
  plot(v_iso_inc, p_hot/1e3, '--', 'LineWidth',1.3, 'DisplayName','Ideal isotherm Th');

  % Isochoric cooling at v = vmax: down from Th(vmax) to Tc(vmax)
  plot([vmax vmax], [p_hot_vmax p_cold_vmax]/1e3, ':', 'LineWidth',1.2, 'DisplayName','Isochoric (cool)');

  % Plot cold isotherm (Tc): vmax -> vmin (cycle direction)
  plot(v_iso_dec, p_cold/1e3, '--', 'LineWidth',1.3, 'DisplayName','Ideal isotherm Tc');

  % Isochoric heating at v = vmin: up from Tc(vmin) to Th(vmin)
  plot([vmin vmin], [p_cold_vmin p_hot_vmin]/1e3, ':', 'LineWidth',1.2, 'DisplayName','Isochoric (heat)');


  % ---- Fig A — PV Diagram (Engine vs. Ideal Stirling) ----
    xlabel('Specific volume v = V_{tot}/m (m^3/kg)','Interpreter','tex');
    ylabel('Pressure p (kPa)','Interpreter','tex');
    title('PV Diagram (Engine vs. Ideal Stirling)');
    box off; 
    grid off;
    legend('Location','best');

  % Save Plot
    w = 5.5;   % width in inches
    h = 2.99;  % height in inches
    set(f, 'Units','inches', 'Position',[1 1 w h]);     % on-screen size
    set(f, 'PaperUnits','inches', 'PaperSize',[w h]);   % paper size
    set(f, 'PaperPosition',[0 0 w h], 'PaperPositionMode','manual');
    print(f, fullfile(OUTDIR,'FigA_PV.png'), '-dpng', '-r600');  % 600 dpi
end

%% TORQUE_POWER  Compute instantaneous torque from virtual work, and power by 2 methods.
function [Trq, Wcyc, P1, P2, Tmean] = torque_power(theta, p, Vtot, RPM, omega, OUTDIR, do_plot)
%   Inputs:
%     theta [Nx1 rad] : crank angle grid
%     p     [Nx1 Pa]  : pressure vs angle
%     Vtot  [Nx1 m^3] : total volume vs angle
%     RPM             : rotational speed (rev/min)
%     omega [rad/s]   : angular speed
%     OUTDIR          : output directory for figures
%   Outputs:
%     Trq  [Nx1 N·m] : torque vs angle
%     Wcyc  [J/cycle] : work over one cycle
%     P1    [W]       : power = Wcyc * f
%     P2    [W]       : power = Tmean * omega
%     Tmean [N·m]     : mean torque over cycle

    % dV/dθ using robust gradient on nonuniform grid
    dVdth = gradient(Vtot, theta);
    
    % T = p * dV/dθ
    Trq = p .* dVdth;
    
    % Calculate work per cycle using trapezoidal rule (0 - 2pi)
    Wcyc = trapz(theta, Trq);        % J per cycle
    
    % Calcualte mean torque
    Tmean = Wcyc / (2*pi);            % N·m
    
    % Two independent power calculations
    f  = RPM/60;                      % Hz
    P1 = Wcyc * f;                    % W  (work per cycle × cycles/sec)
    P2 = Tmean * omega;               % W  (avg torque × angular speed)
    

    % ---- Fig B — Torque vs. Crank Angle ----
    if do_plot
        if nargin < 6 || isempty(OUTDIR), OUTDIR = pwd; end
        if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end
    
        fB = figure('Color','w','Visible','on');  
        plot(theta*180/pi, Trq, 'LineWidth', 1.8);
        xlabel('\theta (deg)'); ylabel('Torque T(\theta) (N·m)');
        title('Torque vs. Crank Angle');
        xlim([0 360]); grid off; box off;
    
        % Set figure size in inches
        w = 5.5; h = 2.99;
        set(fB,'Units','inches','Position',[1 1 w h]);
        set(fB,'PaperUnits','inches','PaperSize',[w h], ...
               'PaperPosition',[0 0 w h],'PaperPositionMode','manual');
    
        % Save high-resolution PNG (600 dpi)
        fn = fullfile(OUTDIR, 'FigB_Torque.png');
        print(fB, fn, '-dpng','-r600');
    end
end

%% ENERGY_AND_INERTIA
function [Ed, Epp, Jreq] = energy_and_inertia(theta, Trq, Tmean, omega, Cf_target, OUTDIR, do_plot)
%   Compute cumulative energy deviation Ed(θ) from load torque and
%   size the flywheel inertia J to meet Cf_target (small-ripple formula).
%
%   dE/dθ = T(θ) - Tload, with Tload = Tmean for steady operation.
%   Epp   = max(Ed) - min(Ed)
%   Jreq  = Epp / (omega^2 * Cf_target)

    theta = ensure_col(theta);
    Trq  = ensure_col(Trq);
    
    Ed    = cumtrapz(theta, Trq - Tmean);       % J (per rad)
    Epp   = max(Ed) - min(Ed);                   % J (peak-to-peak over cycle)
    
    Jreq  = Epp / (omega^2 * Cf_target);         % kg·m^2 (small-ripple sizing)
    
    % plot
    if do_plot && ~isempty(OUTDIR)
        fE = figure('Color','w','Visible','off'); 
        plot(theta*180/pi, Ed, 'LineWidth', 1.8);
        xlabel('\theta (deg)'); ylabel('\DeltaE(\theta) (J)');
        title('Energy Deviation over Cycle');
        xlim([0 360]); grid off; box off;
    
        % Set figure size (inches)
        w = 5.5; h = 2.99;
        set(fE,'Units','inches','Position',[1 1 w h]);
        set(fE,'PaperUnits','inches','PaperSize',[w h], ...
               'PaperPosition',[0 0 w h],'PaperPositionMode','manual');
    
        % Save high-res PNG (600 dpi)
        fn_png = fullfile(OUTDIR,'FigE_EnergyDeviation.png');
        print(fE, fn_png, '-dpng','-r600');
    
        % Optional: save vector PDF for perfectly crisp lines
        fn_pdf = fullfile(OUTDIR,'FigE_EnergyDeviation.pdf');
        print(fE, fn_pdf, '-dpdf');
    
    end
end

%% SIMULATESPEED
function [omega_th, Cf_sim] = simulateSpeed(theta, Trq, Tload, J, omega_mean, OUTDIR)
%   Integrate energy balance to get ω(θ).
%   d/dθ (½ J ω^2) = T(θ) - Tload  ⇒  ω^2(θ) = ω_mean^2 + (2/J)∫(T - Tload) dθ
%
%   Returns ω(θ) and simulated ripple Cf_sim = (ωmax - ωmin) / ωmean,
%   and saves Fig C — Speed vs angle.

    theta    = ensure_col(theta);
    Trq     = ensure_col(Trq);
    Tload    = ensure_col(Tload);
    E        = cumtrapz(theta, Trq - Tload); % J (per rad)
    omega2   = omega_mean^2 + (2/J) * E;      % rad^2/s^2
    omega2   = max(omega2, 0);            % guard
    omega_th = sqrt(omega2);
    
    om_mean  = mean(omega_th);

    % DO NOT USE min(om_mean, 1e-12)
    Cf_sim   = (max(omega_th) - min(omega_th)) / max(om_mean, 1e-12);
    
    % ---- Fig C: Speed vs crank angle ----
    if ~isempty(OUTDIR)
        fC = figure('Color','w','Visible','off'); 
        plot(theta*180/pi, omega_th, 'LineWidth', 1.8);
        xlabel('\theta (deg)'); ylabel('\omega(\theta) (rad/s)');
        title('Speed vs. Crank Angle');
        xlim([0 360]); grid off; box off;
    
        % Set figure size (inches)
        w = 5.5; h = 2.99;
        set(fC,'Units','inches','Position',[1 1 w h]);
        set(fC,'PaperUnits','inches','PaperSize',[w h], ...
               'PaperPosition',[0 0 w h],'PaperPositionMode','manual');
    
        % Save high-res PNG (600 dpi)
        fn_png = fullfile(OUTDIR,'FigC_Speed.png');
        print(fC, fn_png, '-dpng','-r600');
    
        % Optional: also save vector PDF for crisp lines
        fn_pdf = fullfile(OUTDIR,'FigC_Speed.pdf');
        print(fC, fn_pdf, '-dpdf');

    end
end

%% PHASE_SWEEP  Vary displacer phase lead phi and generate:
function out = phase_sweep(phis_deg, theta, geom_base, T, gas, const, spec, OUTDIR)
%   - REQUIRED: Fig D — Work (energy) per cycle vs phase angle
%   - OPTIONAL: Jreq vs phase (design insight)
%
% Inputs:
%   phis_deg : vector of phase angles in degrees (e.g., 60:5:120)
%   theta    : angle grid [rad]
%   geom_base: baseline geometry (fields: r,L,A_p,A_d,Vh0,Vc0,Vr,phi)
%   T        : temps struct with Th, Tc, Tr
%   gas      : struct with m, R
%   const    : struct with omega
%   spec     : struct with RPM, Cs
%   OUTDIR   : output folder
%
% Outputs (struct):
%   out.Wcyc [n x 1] : work per cycle vs phase
%   out.Jreq [n x 1] : required inertia vs phase

    % Ensure column vectors / types
    phis_deg = phis_deg(:);
    nphi     = numel(phis_deg);
    theta    = ensure_col(theta); 
    
    Wcyc = zeros(nphi,1);
    Jreq = zeros(nphi,1);

    for k = 1:nphi
        % Geometry for this phase
        geom = geom_base;
        geom.phi = deg2rad(phis_deg(k));
        
        % Volumes → Pressure → Torque/Power
        [Vh, Vc, Vtot, Vp, Vs] = volumes(theta, geom);
        
        % Use effective hot/cold volumes for the pressure model (beta-type):
        Vh_eff = Vh;
        Vc_eff = Vc + Vp;             % include power-piston volume in the cold space
        
        p = pressure_schmidt(Vh_eff, Vc_eff, geom.Vr, T, gas.m, gas.R);
        % Use a throwaway OUTDIR for per-phi internals to avoid clutter
        [Trq, Wk, ~, ~, Tmean] = torque_power(theta, p, Vtot, spec.RPM, const.omega, [], false);
        
        % Store work per cycle (required for Fig D)
        Wcyc(k) = Wk;

        % Energy deviation → required inertia at this phi
        [~, ~, Jk] = energy_and_inertia(theta, Trq, Tmean, const.omega, spec.Cf, [], false);
        Jreq(k) = Jk;
    end
    
    % ---- Fig D — Work per cycle vs phase ----
    if ~isempty(OUTDIR)
        fD = figure('Color','w','Visible','off'); 
        plot(phis_deg, Wcyc, 'LineWidth', 1.8);
        xlabel('\phi (deg)'); ylabel('Work per cycle W_{cyc} (J)');
        title('Energy/Work per Cycle vs. Phase Angle');
        xlim([min(phis_deg) max(phis_deg)]); grid off; box off;
    
        % Set figure size (inches)
        w = 5.5; h = 2.99;
        set(fD,'Units','inches','Position',[1 1 w h]);
        set(fD,'PaperUnits','inches','PaperSize',[w h], ...
               'PaperPosition',[0 0 w h],'PaperPositionMode','manual');
    
        % Save high-res PNG (600 dpi)
        fn_png = fullfile(OUTDIR,'FigD_PhaseSweep.png');
        print(fD, fn_png, '-dpng','-r600');
    
        % Optional: also save vector PDF
        fn_pdf = fullfile(OUTDIR,'FigD_PhaseSweep.pdf');
        print(fD, fn_pdf, '-dpdf');
    end
    % ---- Fig E — J_{req} vs. Phase Angle ----
    % ---- !NOTE!----
    % In many textbook Stirling models, Jreq(φ) forms a U-shape: 
    %   minimum near φ ≈ 90°, maximum toward 60° and 120°. 
    % This happens when the displacer shuttle volume is large relative to 
    % the clearance volumes, so shifting φ strongly changes how much gas 
    % is hot vs cold during compression/expansion. 
    %
    % In our Appendix B geometry, the displacer swept volume is 
    % Vdisp = 3.0e-5 m^3 (peak-to-peak), while each clearance volume is 
    % Vh0 = Vc0 ≈ 2.35e-4 m^3. That means the shuttle amplitude is only 
    % ~6% of the clearance size. 
    %
    % Because the displacer is small compared to the “dead volumes,” most 
    % of the gas remains in place regardless of φ. As a result, torque ripple 
    % (ΔEpp) changes almost monotonically with φ instead of showing a strong 
    % minimum at 90°. 
    %
    % Physically: the timing (phase angle) has only a weak effect on pressure 
    % distribution when the shuttle volume is small vs clearances. Hence the 
    % Jreq vs φ curve in this model is not U-shaped but drifts steadily across 
    % the sweep range. This is a parameter-driven effect, not a coding error.
    fJ = figure('Color','w','Visible','off'); 
    plot(phis_deg, Jreq, 'LineWidth', 1.8);
    xlabel('\phi (deg)'); ylabel('Required Inertia J_{req} (kg·m^2)');
    title('Required Inertia vs. Phase Angle');
    xlim([min(phis_deg) max(phis_deg)]); 
    grid off; box off;
    % Set figure size (inches)
    w = 5.5; h = 2.99;
    set(fJ,'Units','inches','Position',[1 1 w h]);
    set(fJ,'PaperUnits','inches','PaperSize',[w h], ...
           'PaperPosition',[0 0 w h],'PaperPositionMode','manual');
    % Save high-res PNG (600 dpi)
    fn_png = fullfile(OUTDIR,'FigJ_Jreq_vs_Phase.png');
    print(fJ, fn_png, '-dpng','-r600');

    % Return data
    out.Wcyc = Wcyc;
    out.Jreq = Jreq;
end

%% SIZEFLYWHEEL_MIND  Compute *minimum outer radius* ro (and thus D)
function fly = sizeFlywheel_minD(Jreq, in, omega)
% for a uniform-section annulus with fixed rim thickness t and width w.
% The inner radius is ri = ro - t.
%
% Inputs:
%   Jreq  : required inertia [kg·m^2]
%   in    : struct with fields:
%           rho [kg/m^3], w [m], t [m],
%           yield [Pa] (optional), vtip_max [m/s] (optional)
%   omega : angular speed [rad/s]
%
% Output:
%   fly : struct (ro, ri, w, mass, J, D, v_tip, hoop_sigma, etc.)

  arguments
    Jreq (1,1) double {mustBePositive}
    in struct
    omega (1,1) double {mustBePositive}
  end

  rho = in.rho;    w = in.w;    t = in.t;
  assert(rho>0 && w>0 && t>0, 'rho,w,t must be positive.');
  yield    = getfield_with_default(in,'yield',NaN);
  vtip_max = getfield_with_default(in,'vtip_max',NaN);

  % J(ro) with fixed rim thickness t (ri = ro - t)
  J_of_ro = @(ro) 0.5*rho*pi*w*( ro.^4 - (ro - t).^4 );

  % Find smallest ro > t s.t. J(ro) >= Jreq
  % Bracket: start just above t and expand until we exceed Jreq.
  ro_lo = t * 1.001;
  ro_hi = max(t*1.2, t + 1e-3);
  while J_of_ro(ro_hi) < Jreq
      ro_hi = ro_hi * 1.25;  % expand geometrically
      if ro_hi > 10          % sanity: 20 m diameter is absurd
          error('sizeFlywheel_minD: Jreq too large for reasonable dimensions. Check upstream (Cs, omega, units, Ed).');
      end
  end

  % Root solve F(ro) = J(ro) - Jreq = 0 for ro in (ro_lo, ro_hi)
  F = @(ro) J_of_ro(ro) - Jreq;
  ro = fzero(F, [ro_lo, ro_hi]);

  % Compute remaining properties
  ri   = ro - t;
  mass = rho*pi*w*(ro^2 - ri^2);
  J    = J_of_ro(ro);
  D    = 2*ro;

  % Operating checks (simple & conservative)
  v_tip      = omega*ro;
  hoop_sigma = rho*omega^2*ro^2; % thin-ring approx; conservative at rim

  pass_v = true; if ~isnan(vtip_max), pass_v = (v_tip <= vtip_max); end
  pass_s = true; SF_yield = NaN;
  if ~isnan(yield)
    SF_yield = yield / hoop_sigma;
    pass_s   = SF_yield >= 2.0;    % target SF≥2
  end
    
    fly = struct();
    fly.rho        = rho;
    fly.w          = w;
    fly.t          = t;
    fly.ro         = ro;
    fly.ri         = ri;
    fly.D          = D;
    fly.mass       = mass;
    fly.J          = J;
    fly.Jreq_in    = Jreq;
    fly.v_tip      = v_tip;
    fly.hoop_sigma = hoop_sigma;
    fly.SF_yield   = SF_yield;
    fly.pass_tip   = pass_v;
    fly.pass_yield = pass_s;
end

% styled by chatgpt
function print_flywheel(fly) 
  fprintf('\nFlywheel design (uniform annulus):\n');

  rho   = getfield_with_default(fly,'rho',NaN);
  ro    = getfield_with_default(fly,'ro',NaN);
  ri    = getfield_with_default(fly,'ri',NaN);
  wall  = getfield_with_default(fly,'wall', ro - ri);
  b     = getfield_with_default(fly,'b', getfield_with_default(fly,'w',NaN));
  t     = getfield_with_default(fly,'t', NaN);
  SF    = getfield_with_default(fly,'SF_yield',NaN);
  
  fprintf('  Material density  ρ      = %.0f kg/m^3\n', rho);
  fprintf('  Outer radius  r_o       = %.3f m  (OD = %.1f cm)\n', ro, 2*ro*100);
  fprintf('  Inner radius  r_i       = %.3f m  (ID = %.1f cm)\n', ri, 2*ri*100);
  if ~isnan(t), fprintf('  Rim thickness  t        = %.3f m\n', t); end
  fprintf('  Wall (r_o - r_i)        = %.3f m\n', wall);
  if ~isnan(b), fprintf('  Width/thickness  b      = %.3f m\n', b); end

  fprintf('  Mass                    = %.3f kg\n', fly.mass);
  fprintf('  Inertia  J              = %.6g kg·m^2', fly.J);
  if isfield(fly,'Jreq_in'), fprintf(' (target %.6g)', fly.Jreq_in); end
  fprintf('\n');

  fprintf('  Tip speed v_tip         = %.1f m/s\n', fly.v_tip);
  fprintf('  Hoop stress (thin ring) = %.0f MPa\n', fly.hoop_sigma/1e6);
  fprintf('  Safety Factor           = %.0f\n', SF);

  if isfield(fly,'pass_tip')
    fprintf('  Tip-speed check         : %s\n', ternary(fly.pass_tip,'SAFE','> check'));
  end
  if isfield(fly,'pass_yield')
    fprintf('  Yield-stress check      : %s\n', ternary(fly.pass_yield,'SAFE','> check'));
  end
end


function val = getfield_with_default(s, field, defaultVal)
  if isfield(s, field), val = s.(field); else, val = defaultVal; end
end

function out = ternary(cond, a, b)
  if cond, out = a; else, out = b; end
end
