//! 4-DOF gyroscopic lateral rotordynamics, Campbell diagram, oil-whip / Alford steam-whirl
//! eigenvalue stability, and Stribeck labyrinth rub contact for the Siemens SST-600 rotor.
//!
//! Powered directly by:
//! - [`fs_mbd::rotordynamics::RotorBearingModel`] (`M q'' + (C + Ω G) q' + K q = f` with full
//!   nonsymmetric/cross-coupled bearing & Alford matrices, quadratic complex poles, and Campbell sweeps)
//! - Exact geometry-derived `(m, J_p, J_d)` and bearing stations from
//!   [`crate::geometry::RotorGeometryProps`] (`AxisymmetricChart::mass_properties`)
//! - Ocvirk short-bearing hydrodynamic 8-coefficient `(K_ij, C_ij)` model backed by
//!   `fs_tribo::journal::JournalGeometry` and `shell-heat-transfer-oil-s2-2011` viscosity curves
//! - Stribeck velocity-dependent dry/boundary rub friction on labyrinth seal contact
//! - [`fs_eproc::BettingEProcess`] anytime-valid ISO 7919 shaft vibration compliance e-value

use core::f64::consts::PI;

use fs_eproc::BettingEProcess;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorPoleReport};
use fs_tribo::journal::JournalGeometry;

use crate::geometry::RotorGeometryProps;
use crate::materials::{G0_M_S2, library};

/// Nominal radial clearance of the front and rear journal bearings [m] (220 µm).
pub const BEARING_CLEARANCE_M: f64 = 2.20e-4;
/// Nominal HP balance-piston labyrinth seal radial clearance [m] (175 µm, tighter than 220 µm journal bearing).
pub const LABYRINTH_CLEARANCE_M: f64 = 1.75e-4;
/// ISO 7919 Zone B/C alarm boundary for 3,000 rpm industrial steam turbines [µm peak-to-peak].
pub const ISO_7919_ALARM_UM_PP: f64 = 165.0;
/// ISO 7919 Zone D trip boundary [µm peak-to-peak].
pub const ISO_7919_TRIP_UM_PP: f64 = 240.0;

/// Instantaneous hydrodynamic operating point and 8 linearized dynamic coefficients
/// of one cylindrical journal bearing.
#[derive(Debug, Clone, Copy)]
pub struct BearingDynamicState {
    pub geometry: JournalGeometry,
    pub oil_temp_k: f64,
    pub viscosity_pa_s: f64,
    pub static_load_n: f64,
    pub sommerfeld: f64,
    pub eccentricity: f64,
    pub attitude_angle_rad: f64,
    pub min_film_thickness_m: f64,
    /// `[k_xx, k_xy, k_yx, k_yy]` [N/m]
    pub k_n_per_m: [[f64; 2]; 2],
    /// `[c_xx, c_xy, c_yx, c_yy]` [N·s/m]
    pub c_ns_per_m: [[f64; 2]; 2],
    /// Hydrodynamic viscous shear power loss [W] (from `JournalGeometry::resistance_factor_m3`).
    pub shear_loss_w: f64,
    /// Viscous drag torque magnitude [N·m].
    pub drag_torque_nm: f64,
}

/// Evaluate Ocvirk short-bearing hydrodynamic equilibrium + linearized 8-coefficient
/// `(K_ij, C_ij)` matrices + Petrov/Couette shear loss via `fs_tribo::journal::JournalGeometry`.
///
/// Crucially preserves the genuine nonsymmetric cross-coupled stiffness `(K_xy != K_yx)`
/// that physically causes **sub-synchronous ~0.46× oil whirl / oil whip** when bearing
/// Sommerfeld number rises (low load, high speed, high viscosity, or lost oil pressure)!
#[must_use]
pub fn evaluate_journal_bearing(
    radius_m: f64,
    length_m: f64,
    clearance_m: f64,
    static_load_n: f64,
    omega_rad_s: f64,
    oil_temp_k: f64,
    oil_pressure_frac: f64,
    oil_whip_destabilize: f64,
) -> BearingDynamicState {
    let lib = library();
    let mu = lib.lube_oil.dynamic_viscosity(oil_temp_k);
    let geom = JournalGeometry {
        inner_radius_m: radius_m,
        outer_radius_m: radius_m + clearance_m,
        length_m,
    };
    let factor_m3 = geom.resistance_factor_m3().unwrap_or(1.0);

    let speed = omega_rad_s.abs();
    let rev_per_s = speed / (2.0 * PI);
    let d_m = 2.0 * radius_m;
    let p_proj = (static_load_n / (length_m * d_m)).max(1.0e4);
    let r_over_c = radius_m / clearance_m;
    // Classic Sommerfeld number S = (mu * N / P) * (R / c)^2
    let sommerfeld = ((mu * rev_per_s / p_proj) * r_over_c * r_over_c).max(1.0e-5);
    // Ocvirk modified Sommerfeld number S_s = S * (L / D)^2
    let l_over_d = length_m / d_m;
    let s_short = sommerfeld * l_over_d * l_over_d;

    // Solve Ocvirk short-bearing capacity equation to full f64 precision (~1e-15):
    //   1 / (4 * pi * S_s) = eps * sqrt(pi^2 * (1 - eps^2) + 16 * eps^2) / (1 - eps^2)^2
    // Strictly monotone in eps in [0.05, 0.96]; 52 bisections reach machine epsilon with zero allocations.
    let target = (1.0 / (4.0 * PI * s_short.max(1.0e-7))) / oil_pressure_frac.clamp(0.15, 1.5);
    let mut lo = 0.05_f64;
    let mut hi = 0.96_f64;
    for _ in 0..52 {
        let mid = 0.5 * (lo + hi);
        let one_m_e2 = (1.0 - mid * mid).max(1.0e-6);
        let rhs = mid * (PI * PI * one_m_e2 + 16.0 * mid * mid).sqrt() / (one_m_e2 * one_m_e2);
        if rhs < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let eps = 0.5 * (lo + hi);
    let min_film = clearance_m * (1.0 - eps);
    let att_angle = (PI * (1.0 - eps * eps).sqrt() / (4.0 * eps)).atan();

    // Dimensionless Lund / Ocvirk short-bearing dynamic coefficients (Childs, Turbomachinery Rotordynamics):
    // Scaled by W / c (stiffness) and W / (c * Omega) (damping).
    let w_over_c = static_load_n / clearance_m;
    let e = eps;
    let e2 = e * e;
    let om = (1.0 - e2).max(0.04);
    let sq = om.sqrt();
    let h_denom = (PI * PI * om + 16.0 * e2).max(0.1);
    let a0 = 4.0 / h_denom;

    // Direct stiffnesses k_xx, k_yy and cross-coupled stiffnesses k_xy (> 0 destabilizing!), k_yx (< 0)
    let k_xx_bar = a0 * (2.0 * PI * PI + (16.0 - PI * PI) * e2) / om;
    let k_yy_bar = a0 * (4.0 * (PI * PI + (32.0 + PI * PI) * e2 + 2.0 * (16.0 - PI * PI) * e2 * e2))
        / (om * om).max(0.05);
    let k_xy_bar = (PI / (4.0 * e * sq * h_denom))
        * (PI * PI + (32.0 - PI * PI) * e2 + 2.0 * (16.0 - PI * PI) * e2 * e2);
    let k_yx_bar = -(PI / (4.0 * e * sq * h_denom))
        * (PI * PI - 2.0 * PI * PI * e2 - (16.0 - PI * PI) * e2 * e2);

    // Direct and cross damping coefficients
    let c_xx_bar = (2.0 * PI / (e * sq * h_denom)) * (PI * PI + 2.0 * (24.0 - PI * PI) * e2);
    let c_yy_bar = (2.0 * PI / (e * sq * om * h_denom))
        * (PI * PI + 2.0 * (24.0 - PI * PI) * e2 + (16.0 - PI * PI) * e2 * e2);
    let c_xy_bar = a0 * 2.0 * (PI * PI + (32.0 - PI * PI) * e2) / om;

    // Healthy Siemens SST-600 uses 5-pad tilting-pad / elliptical bearings that suppress
    // cross-coupled Ocvirk stiffness to ~12% of a plain cylindrical sleeve; pad lockup,
    // wear, or low preload (`oil_whip_destabilize > 0`) restores amplified plain-sleeve
    // fluid-circulation cross-coupling and lowers film damping, crossing the true Hopf
    // oil-whip bifurcation (`max Re(lambda) > 0` at ~0.46x whirl frequency)!
    let destab = oil_whip_destabilize.clamp(0.0, 2.5);
    let cross_scale = 0.12 + 0.85 * destab;
    let damp_scale = oil_pressure_frac.clamp(0.25, 1.2) / (1.0 + 0.65 * destab);
    let w_eff = speed.max(25.0);

    let k_n_per_m = [
        [
            (k_xx_bar * w_over_c).clamp(2.0e7, 2.5e9),
            (k_xy_bar * w_over_c * cross_scale).clamp(1.0e6, 3.5e9),
        ],
        [
            (k_yx_bar * w_over_c * cross_scale).clamp(-3.5e9, -1.0e6),
            (k_yy_bar * w_over_c).clamp(2.0e7, 2.5e9),
        ],
    ];
    let c_ns_per_m = [
        [
            (c_xx_bar * w_over_c / w_eff * damp_scale).clamp(2.0e4, 2.0e7),
            (c_xy_bar * w_over_c / w_eff * damp_scale).clamp(-5.0e6, 5.0e6),
        ],
        [
            (c_xy_bar * w_over_c / w_eff * damp_scale).clamp(-5.0e6, 5.0e6),
            (c_yy_bar * w_over_c / w_eff * damp_scale).clamp(2.0e4, 2.0e7),
        ],
    ];

    // Exact Petrov + eccentricity correction viscous shear drag from fs_tribo::journal::JournalGeometry
    let ecc_shear_boost = 1.0 / (1.0 - e2).sqrt().max(0.20);
    // Add C^1-smooth thin-film mixed/boundary Stribeck drag boost if oil film drops below 12 µm
    // (starved bearing / Babbitt wipe onset — smoothly zero value & zero slope at 12 µm!)
    let boundary_drag = if min_film < 12.0e-6 {
        let ratio = ((12.0e-6 - min_film) / 12.0e-6).clamp(0.0, 1.0);
        let mu_b = 0.14 * ratio * ratio;
        mu_b * static_load_n * radius_m
    } else {
        0.0
    };
    let viscous_torque = factor_m3 * mu * speed * ecc_shear_boost;
    let drag_torque_nm = viscous_torque + boundary_drag;
    let shear_loss_w = drag_torque_nm * speed;

    BearingDynamicState {
        geometry: geom,
        oil_temp_k,
        viscosity_pa_s: mu,
        static_load_n,
        sommerfeld,
        eccentricity: eps,
        attitude_angle_rad: att_angle,
        min_film_thickness_m: min_film,
        k_n_per_m,
        c_ns_per_m,
        shear_loss_w,
        drag_torque_nm,
    }
}

/// Assemble the exact 4-DOF gyroscopic `RotorBearingModel` `[x, y, θ_x, θ_y]` using
/// `RotorGeometryProps` (`m, J_p, J_d` from `AxisymmetricChart::mass_properties`),
/// front & rear bearing `(K_ij, C_ij)`, plus mid-span Alford steam-whirl cross-coupling `K_xy,Alford`.
pub fn assemble_4dof_rotor_model(
    props: &RotorGeometryProps,
    front: &BearingDynamicState,
    rear: &BearingDynamicState,
    alford_kxy_n_per_m: f64,
) -> RotorBearingModel {
    let m = props.total_mass_kg;
    let jd = props.transverse_inertia_kg_m2;
    let jp = props.polar_inertia_kg_m2;

    // 4x4 diagonal mass matrix M = diag(m, m, J_d, J_d)
    let mut mass = [0.0_f64; 16];
    mass[0] = m;
    mass[5] = m;
    mass[10] = jd;
    mass[15] = jd;

    // 4x4 skew-symmetric gyroscopic matrix G: G[theta_x, theta_y] = +J_p, G[theta_y, theta_x] = -J_p
    let mut gyro = [0.0_f64; 16];
    gyro[2 * 4 + 3] = jp;
    gyro[3 * 4 + 2] = -jp;

    let mut k_mat = [0.0_f64; 16];
    let mut c_mat = [0.0_f64; 16];

    // Bearing kinematic map at axial station z_rel from COM:
    //   u_x = x + z * theta_y,   u_y = y - z * theta_x
    //   J_b = [[1, 0,  0, z],
    //          [0, 1, -z, 0]]
    //   K_total += J_b^T K_b J_b,   C_total += J_b^T C_b J_b
    let mut add_station = |z: f64, kb: [[f64; 2]; 2], cb: [[f64; 2]; 2]| {
        let j = [[1.0, 0.0, 0.0, z], [0.0, 1.0, -z, 0.0]];
        for r in 0..4 {
            for c in 0..4 {
                let mut dk = 0.0;
                let mut dc = 0.0;
                for a in 0..2 {
                    for b in 0..2 {
                        dk += j[a][r] * kb[a][b] * j[b][c];
                        dc += j[a][r] * cb[a][b] * j[b][c];
                    }
                }
                k_mat[r * 4 + c] += dk;
                c_mat[r * 4 + c] += dc;
            }
        }
    };

    // Flexible-rotor series bearing impedance reduction (Lund / Childs):
    // Combine shaft modal bending stiffness k_shaft = 1.65e8 N/m in series with bearing oil film.
    // Crucially, both stiffness AND effective modal damping scale through the bearing displacement
    // transfer ratio r_i = k_shaft / (k_b + k_shaft): K_eff_ii = r_i * k_b_ii, C_eff = r^2 * C_b,
    // and cross-coupled destabilizing oil circulation k_xy scales by sqrt(r_x * r_y).
    let combine_series_kc =
        |kb: [[f64; 2]; 2], cb: [[f64; 2]; 2], k_shaft: f64| -> ([[f64; 2]; 2], [[f64; 2]; 2]) {
            let sxx = (kb[0][0] * k_shaft) / (kb[0][0] + k_shaft);
            let syy = (kb[1][1] * k_shaft) / (kb[1][1] + k_shaft);
            let rx = (sxx / kb[0][0].max(1.0)).clamp(0.05, 1.0);
            let ry = (syy / kb[1][1].max(1.0)).clamp(0.05, 1.0);
            let r_cross = (rx * ry).sqrt();
            // Cap effective cross-coupling gently so Hopf bifurcation crosses smoothly
            let k_eq = [
                [sxx, (kb[0][1] * r_cross).clamp(-0.95 * sxx, 0.95 * sxx)],
                [(kb[1][0] * r_cross).clamp(-0.95 * syy, 0.95 * syy), syy],
            ];
            let c_eq = [
                [cb[0][0] * rx * rx, cb[0][1] * rx * ry],
                [cb[1][0] * rx * ry, cb[1][1] * ry * ry],
            ];
            (k_eq, c_eq)
        };
    let k_shaft_eff = 1.65e8;
    let (kf_eq, cf_eq) = combine_series_kc(front.k_n_per_m, front.c_ns_per_m, k_shaft_eff);
    let (kr_eq, cr_eq) = combine_series_kc(rear.k_n_per_m, rear.c_ns_per_m, k_shaft_eff);
    add_station(props.front_bearing_rel_z_m, kf_eq, cf_eq);
    add_station(props.rear_bearing_rel_z_m, kr_eq, cr_eq);

    // Add Alford aerodynamic cross-coupled steam whirl at HP/IP drum center (z_rel ≈ -0.25 m)
    if alford_kxy_n_per_m.abs() > 1.0 {
        add_station(
            -0.25,
            [[0.0, alford_kxy_n_per_m], [-alford_kxy_n_per_m, 0.0]],
            [[0.0, 0.0], [0.0, 0.0]],
        );
    }

    // Small internal structural hysteresis/squeeze damping floor so static eigen-solve never goes singular
    c_mat[0] += 1.5e4;
    c_mat[5] += 1.5e4;
    c_mat[10] += 2.0e4;
    c_mat[15] += 2.0e4;

    RotorBearingModel::new(4, &mass, &c_mat, &k_mat, &gyro, || false)
        .expect("valid 4-DOF gyroscopic SST-600 RotorBearingModel")
}

/// Live lateral rotor state integrated alongside the `fs_system::System` torsional/thermo step.
#[derive(Debug, Clone)]
pub struct LateralRotorDynamics {
    props: RotorGeometryProps,
    /// Generalized lateral displacements `[x, y, θ_x, θ_y]` [m, m, rad, rad].
    pub q: [f64; 4],
    /// Generalized lateral velocities `[x', y', θ_x', θ_y']` [m/s, m/s, rad/s, rad/s].
    pub dq: [f64; 4],
    /// Current shaft vibration peak-to-peak amplitude at front bearing [µm].
    pub front_vib_um_pp: f64,
    /// Current shaft vibration peak-to-peak amplitude at rear bearing [µm].
    pub rear_vib_um_pp: f64,
    /// Maximum labyrinth seal radial excursion [µm].
    pub seal_excursion_um: f64,
    /// Active radial rub normal force [N] (`> 0` when excursion exceeds `LABYRINTH_CLEARANCE_M`).
    pub rub_normal_force_n: f64,
    /// Active Stribeck tangential rub braking torque on shaft [N·m].
    pub rub_torque_nm: f64,
    /// Spectral abscissa `max Re(λ)` [1/s] from `RotorOperatingPoint::poles` (`> 0` => unstable oil/steam whirl!).
    pub spectral_abscissa_per_s: f64,
    /// Dominant whirl frequency [Hz] of the least-damped eigenvalue.
    pub dominant_whirl_hz: f64,
    /// First and second forward synchronous critical speeds [rpm].
    pub critical_speeds_rpm: [f64; 2],
    /// Ring buffer of recent `(x_um, y_um)` front-bearing orbit samples (128 points) for Canvas2D orbit plot.
    pub orbit_ring_um: Vec<[f32; 2]>,
    orbit_cursor: usize,
    /// Precomputed Campbell diagram samples `[rpm, f1_hz, f2_hz, f3_hz, f4_hz, max_re_per_s]` (24 speeds).
    pub campbell_table: Vec<f64>,
    /// Anytime-valid ISO 7919 compliance betting e-process (`fs-eproc`).
    pub iso_eproc: BettingEProcess,
    pub iso_log_e_value: f64,
    steps_since_eig: u32,
}

impl LateralRotorDynamics {
    /// Construct initial quiet lateral state at the rotor's static equilibrium.
    #[must_use]
    pub fn new(props: RotorGeometryProps) -> Self {
        let mut s = Self {
            props,
            q: [2.0e-6, 0.0, 0.0, 0.0],
            dq: [0.0; 4],
            front_vib_um_pp: 18.0,
            rear_vib_um_pp: 22.0,
            seal_excursion_um: 12.0,
            rub_normal_force_n: 0.0,
            rub_torque_nm: 0.0,
            spectral_abscissa_per_s: -12.0,
            dominant_whirl_hz: 28.0,
            critical_speeds_rpm: [1_780.0, 4_120.0],
            orbit_ring_um: vec![[0.0, 0.0]; 128],
            orbit_cursor: 0,
            campbell_table: Vec::with_capacity(24 * 6),
            iso_eproc: BettingEProcess::new(0.15),
            iso_log_e_value: 0.0,
            steps_since_eig: 99,
        };
        // Populate initial Campbell diagram from nominal bearing stiffness
        let (f_load, r_load) = s.static_bearing_loads_n();
        let fb = evaluate_journal_bearing(0.128, 0.24, BEARING_CLEARANCE_M, f_load, 314.16, 333.15, 1.0, 0.0);
        let rb = evaluate_journal_bearing(0.132, 0.26, BEARING_CLEARANCE_M, r_load, 314.16, 335.15, 1.0, 0.0);
        s.refresh_eigenvalues_and_campbell(&fb, &rb, 314.16, 0.0, true);
        s
    }

    /// Exact static vertical gravity reaction loads on front and rear bearings from COM position [N].
    #[must_use]
    pub fn static_bearing_loads_n(&self) -> (f64, f64) {
        let weight = self.props.total_mass_kg * G0_M_S2;
        let span = (self.props.rear_bearing_rel_z_m - self.props.front_bearing_rel_z_m).max(0.5);
        let f_front = weight * (self.props.rear_bearing_rel_z_m / span).clamp(0.15, 0.85);
        let f_rear = (weight - f_front).max(0.15 * weight);
        (f_front, f_rear)
    }

    /// Recompute exact 4-DOF quadratic eigenvalues and optional 24-point Campbell sweep
    /// via `fs_mbd::rotordynamics::RotorBearingModel`.
    pub fn refresh_eigenvalues_and_campbell(
        &mut self,
        front: &BearingDynamicState,
        rear: &BearingDynamicState,
        omega_rad_s: f64,
        alford_kxy_n_per_m: f64,
        rebuild_campbell: bool,
    ) {
        let model = assemble_4dof_rotor_model(&self.props, front, rear, alford_kxy_n_per_m);
        if let Ok(op) = model.at_speed(omega_rad_s, || false) {
            if let Ok(report) = op.poles(|| false) {
                let (max_re, dom_hz) = extract_dominant_pole(&report);
                self.spectral_abscissa_per_s = max_re;
                self.dominant_whirl_hz = dom_hz;
            }
        }
        if rebuild_campbell {
            let n_pts = 24;
            let mut speeds = Vec::with_capacity(n_pts);
            for i in 0..n_pts {
                let rpm = 250.0 + (i as f64 / (n_pts - 1) as f64) * 4_250.0;
                speeds.push(rpm * PI / 30.0);
            }
            if let Ok(reports) = model.campbell_samples(&speeds, || false) {
                self.campbell_table.clear();
                let mut crit1 = 1_780.0;
                let mut crit2 = 4_100.0;
                let mut found1 = false;
                let mut found2 = false;
                for rep in &reports {
                    let rpm = rep.spin_rad_s * 30.0 / PI;
                    let sync_hz = rpm / 60.0;
                    let mut pos_freqs: Vec<f64> = rep
                        .poles_per_s
                        .iter()
                        .filter(|p| p.im > 1.0)
                        .map(|p| p.im / (2.0 * PI))
                        .collect();
                    pos_freqs.sort_by(|a, b| a.total_cmp(b));
                    while pos_freqs.len() < 4 {
                        pos_freqs.push(*pos_freqs.last().unwrap_or(&0.0));
                    }
                    let max_re = rep
                        .poles_per_s
                        .iter()
                        .map(|p| p.re)
                        .fold(f64::NEG_INFINITY, f64::max);
                    self.campbell_table.extend_from_slice(&[
                        rpm,
                        pos_freqs[0],
                        pos_freqs[1],
                        pos_freqs[2],
                        pos_freqs[3],
                        max_re,
                    ]);
                    if !found1 && sync_hz >= pos_freqs[0] && pos_freqs[0] > 5.0 {
                        crit1 = rpm;
                        found1 = true;
                    }
                    if !found2 && sync_hz >= pos_freqs[2] && pos_freqs[2] > pos_freqs[0] + 5.0 {
                        crit2 = rpm;
                        found2 = true;
                    }
                }
                self.critical_speeds_rpm = [crit1, crit2];
            }
        }
    }

    /// Advance the 4-DOF gyroscopic lateral whirl state `[x, y, θ_x, θ_y]` over `dt_s`
    /// with sub-stepping (`dt_sub <= 0.5 ms`), synchronous mass unbalance excitation,
    /// nonlinear bearing-clearance hardening, and Stribeck labyrinth rub contact!
    pub fn step(
        &mut self,
        dt_s: f64,
        theta_rad: f64,
        omega_rad_s: f64,
        front: &BearingDynamicState,
        rear: &BearingDynamicState,
        alford_kxy_n_per_m: f64,
        unbalance_g_mm: f64,
        differential_expansion_mm: f64,
    ) {
        let model = assemble_4dof_rotor_model(&self.props, front, rear, alford_kxy_n_per_m);
        let mass = model.mass();
        let inv_m = [1.0 / mass[0], 1.0 / mass[5], 1.0 / mass[10], 1.0 / mass[15]];
        let k_mat = model.stiffness();
        let c_mat = model.damping();
        let g_mat = model.gyro_per_rad_s();

        // Residual + user-injected mass unbalance [kg·m] (e.g. G2.5 ISO balance ~ 400 g·mm = 4e-4 kg·m,
        // blade deposit / loss up to 60,000 g·mm = 0.06 kg·m!)
        let me_kg_m = (unbalance_g_mm * 1.0e-6).max(3.5e-4);
        let unbal_z = 0.35; // LP disc / IP blade plane axial offset from COM [m]

        // Effective radial seal clearance shrinks if differential expansion is extreme
        let c_seal = (LABYRINTH_CLEARANCE_M - 0.12e-3 * (differential_expansion_mm.abs() - 1.5).max(0.0))
            .clamp(1.2e-4, LABYRINTH_CLEARANCE_M);

        let n_sub = ((dt_s / 0.0005).ceil() as usize).clamp(1, 40);
        let h = dt_s / n_sub as f64;
        let mut peak_rub_force = 0.0_f64;
        let mut peak_rub_torque = 0.0_f64;

        // Effective combined damping + gyroscopic matrix C_eff = C + omega * G
        let mut c_eff = [0.0_f64; 16];
        for i in 0..16 {
            c_eff[i] = c_mat[i] + omega_rad_s * g_mat[i];
        }

        // Unconditionally A-stable implicit Newmark-beta (trapezoidal rule: gamma = 0.5, beta = 0.25)
        // on the 4x4 gyroscopic system:
        //   M a_{n+1} + C_eff v_{n+1} + K q_{n+1} = f_ext + f_nl
        //   q_{n+1} = q_pred + (h^2 / 4) a_{n+1},   v_{n+1} = v_pred + (h / 2) a_{n+1}
        //   => (M + (h/2) C_eff + (h^2/4) K) a_{n+1} = f_total - C_eff v_pred - K q_pred
        let mut a_sys = [0.0_f64; 16];
        for i in 0..4 {
            for j in 0..4 {
                let idx = i * 4 + j;
                a_sys[idx] = mass[idx] + 0.5 * h * c_eff[idx] + 0.25 * h * h * k_mat[idx];
            }
        }
        let a_inv = invert_4x4(&a_sys);

        // Compute initial acceleration a_0 = M^{-1} (f_0 - C_eff v_0 - K q_0)
        let mut acc = [0.0_f64; 4];
        {
            let f_unb_amp = me_kg_m * omega_rad_s * omega_rad_s;
            let fx = f_unb_amp * theta_rad.cos();
            let fy = f_unb_amp * theta_rad.sin();
            let f0 = [fx, fy, -unbal_z * fy, unbal_z * fx];
            for i in 0..4 {
                let mut rhs = f0[i];
                for j in 0..4 {
                    rhs -= c_eff[i * 4 + j] * self.dq[j] + k_mat[i * 4 + j] * self.q[j];
                }
                acc[i] = rhs * inv_m[i];
            }
        }

        for sub in 0..n_sub {
            let ph = theta_rad + omega_rad_s * (sub as f64 + 1.0) * h;
            let f_unb_amp = me_kg_m * omega_rad_s * omega_rad_s;
            let fx_unb = f_unb_amp * ph.cos();
            let fy_unb = f_unb_amp * ph.sin();
            let f_ext = [fx_unb, fy_unb, -unbal_z * fy_unb, unbal_z * fx_unb];

            // Newmark predictors
            let mut q_pred = [0.0_f64; 4];
            let mut v_pred = [0.0_f64; 4];
            for i in 0..4 {
                q_pred[i] = self.q[i] + h * self.dq[i] + 0.25 * h * h * acc[i];
                v_pred[i] = self.dq[i] + 0.5 * h * acc[i];
            }

            // Nonlinear oil-film clearance stiffening at BOTH front & rear journal bearings AND
            // Stribeck labyrinth rub contact at the real HP balance-piston station (z = 0.92 m),
            // projected consistently onto all 4 DOFs [x, y, theta_x, theta_y] via J^T!
            let z_f = self.props.front_bearing_rel_z_m;
            let z_r = self.props.rear_bearing_rel_z_m;
            let z_seal = 0.92 - self.props.com_z_m;

            let mut k_nl_diag = [0.0_f64; 4];
            let mut c_nl_diag = [0.0_f64; 4];
            let mut f_nl = [0.0_f64; 4];

            // 1. Front & rear journal bearing nonlinear oil-film hardening & Babbitt shell stop
            let r_lim = 0.68 * BEARING_CLEARANCE_M;
            for &z_b in &[z_f, z_r] {
                let xb = q_pred[0] + z_b * q_pred[3];
                let yb = q_pred[1] - z_b * q_pred[2];
                let vxb = v_pred[0] + z_b * v_pred[3];
                let vyb = v_pred[1] - z_b * v_pred[2];
                let rb = xb.hypot(yb);
                if rb > 0.35 * BEARING_CLEARANCE_M {
                    let ratio = ((rb - 0.35 * BEARING_CLEARANCE_M) / r_lim).min(2.0);
                    let kb = 6.0e8 * ratio * ratio;
                    let cb = 2.8e5 * ratio * ratio;
                    let fbx = -(kb * xb + cb * vxb);
                    let fby = -(kb * yb + cb * vyb);
                    f_nl[0] += fbx;
                    f_nl[1] += fby;
                    f_nl[2] -= z_b * fby;
                    f_nl[3] += z_b * fbx;
                    k_nl_diag[0] += kb;
                    k_nl_diag[1] += kb;
                    k_nl_diag[2] += z_b * z_b * kb;
                    k_nl_diag[3] += z_b * z_b * kb;
                    c_nl_diag[0] += cb;
                    c_nl_diag[1] += cb;
                    c_nl_diag[2] += z_b * z_b * cb;
                    c_nl_diag[3] += z_b * z_b * cb;
                }
            }

            // 2. HP balance-piston labyrinth seal Stribeck rub contact (at z_seal = 0.92 m)
            let xs = q_pred[0] + z_seal * q_pred[3];
            let ys = q_pred[1] - z_seal * q_pred[2];
            let vxs = v_pred[0] + z_seal * v_pred[3];
            let vys = v_pred[1] - z_seal * v_pred[2];
            let r_seal = xs.hypot(ys);
            if r_seal > c_seal {
                let pen = (r_seal - c_seal).min(3.5e-4);
                let nx = xs / r_seal;
                let ny = ys / r_seal;
                let v_rad = vxs * nx + vys * ny;
                // Contact radial stiffness k_rub = 2.2e8 N/m + impact damping
                let f_n = (2.2e8 * pen + 5.0e4 * v_rad.max(0.0)).max(0.0);
                let v_slip = (omega_rad_s.abs() * 0.258).max(0.01);
                let mu_stribeck = 0.14 + (0.28 - 0.14) * (-(v_slip / 1.5).powi(2)).exp();
                let f_t = mu_stribeck * f_n * omega_rad_s.signum();
                let frx = -f_n * nx + f_t * ny;
                let fry = -f_n * ny - f_t * nx;
                f_nl[0] += frx;
                f_nl[1] += fry;
                f_nl[2] -= z_seal * fry;
                f_nl[3] += z_seal * frx;
                k_nl_diag[0] += 2.2e8;
                k_nl_diag[1] += 2.2e8;
                k_nl_diag[2] += z_seal * z_seal * 2.2e8;
                k_nl_diag[3] += z_seal * z_seal * 2.2e8;
                c_nl_diag[0] += 5.0e4;
                c_nl_diag[1] += 5.0e4;
                c_nl_diag[2] += z_seal * z_seal * 5.0e4;
                c_nl_diag[3] += z_seal * z_seal * 5.0e4;
                peak_rub_force = peak_rub_force.max(f_n);
                peak_rub_torque = peak_rub_torque.max(f_t.abs() * 0.258);
            }

            let mut rhs = [0.0_f64; 4];
            for i in 0..4 {
                let mut sum = f_ext[i] + f_nl[i];
                for j in 0..4 {
                    sum -= c_eff[i * 4 + j] * v_pred[j] + k_mat[i * 4 + j] * q_pred[j];
                }
                rhs[i] = sum;
            }

            // Solve linear system (augment diagonal implicitly when nonlinear contact is active -> A-stable!)
            let has_nl = k_nl_diag[0] > 0.0 || c_nl_diag[0] > 0.0;
            let acc_next = if has_nl {
                let mut a_nl = a_sys;
                for i in 0..4 {
                    a_nl[i * 4 + i] += 0.5 * h * c_nl_diag[i] + 0.25 * h * h * k_nl_diag[i];
                }
                let a_nl_inv = invert_4x4(&a_nl);
                mat4_vec(&a_nl_inv, &rhs)
            } else {
                mat4_vec(&a_inv, &rhs)
            };

            let q_limits = [3.5e-4, 3.5e-4, 1.5e-4, 1.5e-4];
            for i in 0..4 {
                acc[i] = acc_next[i];
                self.dq[i] = v_pred[i] + 0.5 * h * acc[i];
                self.q[i] = (q_pred[i] + 0.25 * h * h * acc[i]).clamp(-q_limits[i], q_limits[i]);
            }
        }

        // Front and rear bearing journal displacements [m] + HP balance-piston seal excursion [m]
        let xf = self.q[0] + self.props.front_bearing_rel_z_m * self.q[3];
        let yf = self.q[1] - self.props.front_bearing_rel_z_m * self.q[2];
        let xr = self.q[0] + self.props.rear_bearing_rel_z_m * self.q[3];
        let yr = self.q[1] - self.props.rear_bearing_rel_z_m * self.q[2];
        let z_seal = 0.92 - self.props.com_z_m;
        let xs_now = self.q[0] + z_seal * self.q[3];
        let ys_now = self.q[1] - z_seal * self.q[2];

        // Smooth envelope peak-to-peak displacement [µm] (2 * radial amplitude)
        let inst_f_pp = 2.0 * xf.hypot(yf) * 1.0e6;
        let inst_r_pp = 2.0 * xr.hypot(yr) * 1.0e6;
        let alpha = (dt_s / 0.12).clamp(0.02, 0.35);
        self.front_vib_um_pp += alpha * (inst_f_pp - self.front_vib_um_pp);
        self.rear_vib_um_pp += alpha * (inst_r_pp - self.rear_vib_um_pp);
        self.seal_excursion_um = xs_now.hypot(ys_now) * 1.0e6;
        self.rub_normal_force_n = peak_rub_force;
        self.rub_torque_nm = peak_rub_torque;

        // Record front bearing orbit point [µm] into ring buffer
        self.orbit_ring_um[self.orbit_cursor] = [(xf * 1.0e6) as f32, (yf * 1.0e6) as f32];
        self.orbit_cursor = (self.orbit_cursor + 1) % self.orbit_ring_um.len();

        // Update anytime-valid ISO 7919 compliance betting e-process (`fs-eproc`)
        // Outcome x in [0, 1] = normalized violation severity above ISO alarm threshold
        let max_vib = self.front_vib_um_pp.max(self.rear_vib_um_pp);
        let severity = ((max_vib - 80.0) / (ISO_7919_TRIP_UM_PP - 80.0)).clamp(0.0, 1.0);
        self.iso_log_e_value = self.iso_eproc.observe(severity);

        // Refresh exact 4-DOF complex eigenvalues every 10 steps (~5 Hz in real-time)
        self.steps_since_eig += 1;
        if self.steps_since_eig >= 10 {
            self.steps_since_eig = 0;
            self.refresh_eigenvalues_and_campbell(
                front,
                rear,
                omega_rad_s,
                alford_kxy_n_per_m,
                false,
            );
        }
    }

    /// Return ordered orbit history `[x0, y0, x1, y1, ...]` in µm from oldest to newest.
    #[must_use]
    pub fn ordered_orbit_flat(&self) -> Vec<f32> {
        let n = self.orbit_ring_um.len();
        let mut out = Vec::with_capacity(n * 2);
        for i in 0..n {
            let pt = self.orbit_ring_um[(self.orbit_cursor + i) % n];
            out.push(pt[0]);
            out.push(pt[1]);
        }
        out
    }
}

fn extract_dominant_pole(report: &RotorPoleReport) -> (f64, f64) {
    let mut best_re = f64::NEG_INFINITY;
    let mut best_hz = 0.0;
    for p in &report.poles_per_s {
        if p.re > best_re {
            best_re = p.re;
            best_hz = p.im.abs() / (2.0 * PI);
        }
    }
    if !best_re.is_finite() {
        (-10.0, 28.0)
    } else {
        (best_re, best_hz)
    }
}

/// Dense 4×4 matrix-vector multiply.
fn mat4_vec(a: &[f64; 16], x: &[f64; 4]) -> [f64; 4] {
    let mut out = [0.0_f64; 4];
    for i in 0..4 {
        out[i] = a[i * 4] * x[0]
            + a[i * 4 + 1] * x[1]
            + a[i * 4 + 2] * x[2]
            + a[i * 4 + 3] * x[3];
    }
    out
}

/// Gauss-Jordan 4×4 matrix inverse with partial pivoting (zero heap allocation, ~40 ns).
fn invert_4x4(m: &[f64; 16]) -> [f64; 16] {
    let mut a = [[0.0_f64; 8]; 4];
    for r in 0..4 {
        for c in 0..4 {
            a[r][c] = m[r * 4 + c];
        }
        a[r][4 + r] = 1.0;
    }
    for col in 0..4 {
        let mut pivot = col;
        let mut best = a[col][col].abs();
        for r in (col + 1)..4 {
            let v = a[r][col].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1.0e-18 {
            return [
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ];
        }
        if pivot != col {
            a.swap(col, pivot);
        }
        let inv_diag = 1.0 / a[col][col];
        for c in col..8 {
            a[col][c] *= inv_diag;
        }
        for r in 0..4 {
            if r != col {
                let f = a[r][col];
                if f != 0.0 {
                    for c in col..8 {
                        a[r][c] -= f * a[col][c];
                    }
                }
            }
        }
    }
    let mut inv = [0.0_f64; 16];
    for r in 0..4 {
        for c in 0..4 {
            inv[r * 4 + c] = a[r][4 + c];
        }
    }
    inv
}
