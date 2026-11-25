//! Real 2D Axisymmetric `(r, z)` Finite-Volume Thermal Conduction, Fourier Heat-Flux,
//! and Coupled Centrifugal + Pressure-Hoop + Thermal-Gradient Thermo-Elastic Stress PDE Solver
//! for the steam turbine model (`64 × 24 = 1,536` continuum cells).
//!
//! Solves the cylindrical heat equation directly on the physical cross-section:
//!
//!   ρ c_p ∂T/∂t = (1/r) ∂(r k(T) ∂T/∂r)/∂r + ∂(k(T) ∂T/∂z)/∂z + q'''(r, z)
//!
//! with real `matdb` temperature-dependent conductivity `k(T)` & heat capacity `c_p`,
//! convective boundary coupling to the live 17-stage IAPWS-IF97 steam flowpath `T_steam(z)`,
//! bearing shear heat injection at `z = 0.25, 4.62 m`, labyrinth rub heating at `z = 0.92 m`,
//! and ambient natural convection + radiation on the outer casing wall.
//!
//! Exports a flat `64 × 24 × 8` `f32` continuum tensor uploaded directly to two WebGL
//! `FloatType DataTexture`s for continuous per-pixel GPU false-color shaders, Fourier heat
//! flux & stress vector fields, and live 3D raycast surface probing.

use crate::materials::{STEFAN_BOLTZMANN, library};
use crate::stage_meanline::{StageMeanlineReport, TurbineDesignParams};

/// Axial grid cells along the turbine centerline (`z ∈ [GRID_Z_MIN_M, GRID_Z_MAX_M]`).
pub const GRID_NZ: usize = 64;
/// Radial grid cells from shaft centerline outward (`r ∈ [0.0, GRID_R_MAX_M]`).
pub const GRID_NR: usize = 24;
/// Total continuum finite-volume cells (`64 × 24 = 1,536`).
pub const GRID_CELLS: usize = GRID_NZ * GRID_NR;
/// Floats exported per cell: `[T_K, sigma_vM_MPa, P_bar, Mach, q_r_kWm2, q_z_kWm2, wetness, creep_log]`.
pub const GRID_STRIDE: usize = 8;

/// Axial domain start [m] (includes front coupling & front bearing pedestal).
pub const GRID_Z_MIN_M: f64 = -0.15;
/// Axial domain end [m] (includes full axial exhaust diffuser exit).
pub const GRID_Z_MAX_M: f64 = 5.55;
/// Radial domain outer bound [m] (covers giant 1.38 m exhaust hood outer flange).
pub const GRID_R_MAX_M: f64 = 1.52;

/// Material region tag of each `(iz, ir)` cell in the axisymmetric cross-section.
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellRegion {
    SteamOrAir = 0,
    RotorForging = 1,
    BladesAnnulus = 2,
    InnerCarrier = 3,
    OuterCasing = 4,
    ExhaustHood = 5,
    BearingPedestal = 6,
}

/// Live 2D axisymmetric `(r, z)` continuum thermal, stress, and steam-path field solver.
#[derive(Debug, Clone)]
pub struct AxisymmetricContinuum2D {
    pub region: [CellRegion; GRID_CELLS],
    pub temp_k: [f64; GRID_CELLS],
    pub sigma_vm_mpa: [f64; GRID_CELLS],
    pub p_bar: [f64; GRID_CELLS],
    pub mach: [f64; GRID_CELLS],
    pub q_r_kw_m2: [f64; GRID_CELLS],
    pub q_z_kw_m2: [f64; GRID_CELLS],
    pub wetness: [f64; GRID_CELLS],
    pub creep_log: [f64; GRID_CELLS],
    /// Peak metal von Mises stress across rotor & casing [MPa].
    pub peak_metal_stress_mpa: f64,
    /// Peak HP casing wall hoop + thermal stress [MPa].
    pub casing_hoop_stress_mpa: f64,
    /// Peak Fourier heat flux magnitude `|q|` [kW/m²].
    pub peak_heat_flux_kw_m2: f64,
    /// Scratch buffer for explicit/semi-implicit FVM sub-stepping (zero heap allocs per step).
    scratch_t: [f64; GRID_CELLS],
    /// Packed `f32` export buffer (`64 × 24 × 8`).
    export_buf: Vec<f32>,
}

impl AxisymmetricContinuum2D {
    /// Construct and classify the 64 × 24 axisymmetric cross-section of the model turbine.
    #[must_use]
    pub fn new(cold_start: bool, design: &TurbineDesignParams) -> Self {
        let mut region = [CellRegion::SteamOrAir; GRID_CELLS];
        let mut temp_k = [298.15; GRID_CELLS];
        let wall_m = design.clamped().casing_wall_mm * 1.0e-3;
        let dz = (GRID_Z_MAX_M - GRID_Z_MIN_M) / GRID_NZ as f64;
        let dr = GRID_R_MAX_M / GRID_NR as f64;

        for ir in 0..GRID_NR {
            let r = (ir as f64 + 0.5) * dr;
            for iz in 0..GRID_NZ {
                let z = GRID_Z_MIN_M + (iz as f64 + 0.5) * dz;
                let idx = ir * GRID_NZ + iz;
                let reg = classify_section_cell(z, r, wall_m);
                region[idx] = reg;
                temp_k[idx] = initial_temp_for_cell(z, r, reg, cold_start);
            }
        }

        let mut this = Self {
            region,
            temp_k,
            sigma_vm_mpa: [0.0; GRID_CELLS],
            p_bar: [1.0; GRID_CELLS],
            mach: [0.0; GRID_CELLS],
            q_r_kw_m2: [0.0; GRID_CELLS],
            q_z_kw_m2: [0.0; GRID_CELLS],
            wetness: [0.0; GRID_CELLS],
            creep_log: [-16.0; GRID_CELLS],
            peak_metal_stress_mpa: 140.0,
            casing_hoop_stress_mpa: 110.0,
            peak_heat_flux_kw_m2: 10.0,
            scratch_t: temp_k,
            export_buf: vec![0.0_f32; GRID_CELLS * GRID_STRIDE],
        };
        this.pack_export();
        this
    }

    /// Re-classify casing thickness when engineer edits `casing_wall_mm` in the Design Studio.
    pub fn update_geometry_mask(&mut self, design: &TurbineDesignParams) {
        let wall_m = design.clamped().casing_wall_mm * 1.0e-3;
        let dz = (GRID_Z_MAX_M - GRID_Z_MIN_M) / GRID_NZ as f64;
        let dr = GRID_R_MAX_M / GRID_NR as f64;
        for ir in 0..GRID_NR {
            let r = (ir as f64 + 0.5) * dr;
            for iz in 0..GRID_NZ {
                let z = GRID_Z_MIN_M + (iz as f64 + 0.5) * dz;
                self.region[ir * GRID_NZ + iz] = classify_section_cell(z, r, wall_m);
            }
        }
    }

    /// Advance the 2D cylindrical thermal conduction PDE + compute coupled thermo-elastic/centrifugal
    /// stress field, Fourier heat flux `q = -k ∇T`, pressure, Mach, wetness, and Norton creep!
    pub fn step(
        &mut self,
        dt_s: f64,
        meanline: &StageMeanlineReport,
        design: &TurbineDesignParams,
        omega_rad_s: f64,
        p_chest_bar: f64,
        t_front_brg_k: f64,
        t_rear_brg_k: f64,
        rub_power_w: f64,
        m_dot_live_kg_s: f64,
    ) {
        let lib = library();
        let dz = (GRID_Z_MAX_M - GRID_Z_MIN_M) / GRID_NZ as f64;
        let dr = GRID_R_MAX_M / GRID_NR as f64;
        let inv_dz2 = 1.0 / (dz * dz);
        let inv_dr2 = 1.0 / (dr * dr);
        let t_amb = 298.15_f64;
        let flow_h_scale = 0.20 + 0.80 * (m_dot_live_kg_s.max(0.0) / 50.0).clamp(0.0, 2.0).powf(0.8);

        // Precompute steam-path properties along all 64 axial stations once per step
        let mut z_p_bar = [0.0_f64; GRID_NZ];
        let mut z_t_k = [300.0_f64; GRID_NZ];
        let mut z_mach = [0.0_f64; GRID_NZ];
        let mut z_wet = [0.0_f64; GRID_NZ];
        for iz in 0..GRID_NZ {
            let z = GRID_Z_MIN_M + (iz as f64 + 0.5) * dz;
            let (pb, tk, ma, wet, _, _) = meanline.sample_flowpath_at_z(z);
            z_p_bar[iz] = pb;
            z_t_k[iz] = tk;
            z_mach[iz] = ma;
            z_wet[iz] = wet;
        }

        // Thermal diffusivity of CrMoV steel α = k/(ρ cp) ≈ 7.8e-6 m²/s; explicit CFL on dr = 0.063 m
        // allows dt <= ~80 s; the stiffest term (bearing-pedestal exchange, 38 kW/m³K into bronze) allows
        // ~90 s. Sub-steps of at most 0.25 s therefore keep a wide margin, and the whole `dt_s` is
        // covered (the plant calls with up to 2 s; a clamp here would drop field time).
        let dt_s = dt_s.clamp(1.0e-4, 2.0);
        let n_sub = ((dt_s / 0.25).ceil() as usize).max(2);
        let h_dt = dt_s / n_sub as f64;

        for _ in 0..n_sub {
            self.scratch_t.copy_from_slice(&self.temp_k);
            for ir in 0..GRID_NR {
                let r = (ir as f64 + 0.5) * dr;
                let r_p = (ir as f64 + 1.0) * dr;
                let r_m = (ir as f64) * dr;
                for iz in 0..GRID_NZ {
                    let idx = ir * GRID_NZ + iz;
                    let reg = self.region[idx];
                    let t_c = self.temp_k[idx];
                    let z = GRID_Z_MIN_M + (iz as f64 + 0.5) * dz;

                    if reg == CellRegion::SteamOrAir {
                        // Steam path tracks local 17-stage steam T(z); external air stays ambient
                        let tgt = if (1.05..=5.45).contains(&z) && r < steam_outer_radius(z) {
                            z_t_k[iz]
                        } else {
                            t_amb
                        };
                        self.scratch_t[idx] = t_c + (tgt - t_c) * (1.0 - (-h_dt * 12.0).exp());
                        continue;
                    }

                    // Material conductivity k(T) & volumetric heat capacity ρ c_p from matdb
                    let (k_c, rho_cp) = match reg {
                        CellRegion::RotorForging | CellRegion::BladesAnnulus => (
                            lib.rotor_steel.conductivity_at(t_c),
                            lib.rotor_steel.rho_kg_m3 * lib.rotor_steel.cp_j_kg_k,
                        ),
                        CellRegion::InnerCarrier | CellRegion::OuterCasing => (
                            lib.hp_casing_steel.conductivity_at(t_c),
                            lib.hp_casing_steel.rho_kg_m3 * lib.hp_casing_steel.cp_j_kg_k,
                        ),
                        CellRegion::ExhaustHood => (
                            lib.lp_exhaust_steel.conductivity_at(t_c),
                            lib.lp_exhaust_steel.rho_kg_m3 * lib.lp_exhaust_steel.cp_j_kg_k,
                        ),
                        CellRegion::BearingPedestal => (
                            lib.bearing_bronze.conductivity_at(t_c),
                            lib.bearing_bronze.rho_kg_m3 * lib.bearing_bronze.cp_j_kg_k,
                        ),
                        CellRegion::SteamOrAir => unreachable!(),
                    };

                    // Cylindrical finite-volume conduction fluxes across 4 faces
                    let mut lap_k = 0.0_f64;

                    // +r face
                    if ir + 1 < GRID_NR {
                        let n_idx = (ir + 1) * GRID_NZ + iz;
                        if self.region[n_idx] != CellRegion::SteamOrAir {
                            lap_k += (r_p / r) * k_c * (self.temp_k[n_idx] - t_c) * inv_dr2;
                        } else if r > 0.60 {
                            // External ambient natural convection + Stefan-Boltzmann radiation
                            let q_amb = 14.0 * (t_c - t_amb)
                                + 0.78 * STEFAN_BOLTZMANN * (t_c.powi(4) - t_amb.powi(4));
                            lap_k -= (r_p / r) * q_amb / dr;
                        } else {
                            // Internal convective steam-to-rotor boundary
                            let h_conv = 1_400.0 * flow_h_scale;
                            lap_k += (r_p / r) * h_conv * (z_t_k[iz] - t_c) / dr;
                        }
                    } else {
                        let q_amb = 14.0 * (t_c - t_amb);
                        lap_k -= q_amb / dr;
                    }

                    // -r face (r = 0 symmetry boundary has zero flux r_m = 0!)
                    if ir > 0 {
                        let n_idx = (ir - 1) * GRID_NZ + iz;
                        if self.region[n_idx] != CellRegion::SteamOrAir {
                            lap_k += (r_m / r) * k_c * (self.temp_k[n_idx] - t_c) * inv_dr2;
                        } else {
                            // Internal steam-to-casing convective heat transfer!
                            let h_conv = if z < 2.2 { 1_850.0 } else if z < 3.6 { 950.0 } else { 420.0 } * flow_h_scale;
                            lap_k += (r_m / r) * h_conv * (z_t_k[iz] - t_c) / dr;
                        }
                    }

                    // +z and -z axial conduction faces
                    if iz + 1 < GRID_NZ {
                        let n_idx = ir * GRID_NZ + (iz + 1);
                        if self.region[n_idx] != CellRegion::SteamOrAir {
                            lap_k += k_c * (self.temp_k[n_idx] - t_c) * inv_dz2;
                        }
                    }
                    if iz > 0 {
                        let n_idx = ir * GRID_NZ + (iz - 1);
                        if self.region[n_idx] != CellRegion::SteamOrAir {
                            lap_k += k_c * (self.temp_k[n_idx] - t_c) * inv_dz2;
                        }
                    }

                    // Volumetric heat sources:
                    // - Blades immersed directly in 17-stage steam flowpath
                    // - Bearing pedestals coupled to live journal bearing temperatures
                    // - Balance-piston labyrinth rub friction heat injection at z ≈ 0.92 m!
                    let mut q_vol = 0.0_f64;
                    if reg == CellRegion::BladesAnnulus {
                        let h_blade_vol = 22_000.0 * flow_h_scale;
                        q_vol += h_blade_vol * (z_t_k[iz] - t_c);
                    } else if reg == CellRegion::BearingPedestal {
                        let t_brg = if z < 2.0 { t_front_brg_k } else { t_rear_brg_k };
                        q_vol += 38_000.0 * (t_brg - t_c);
                    }
                    if rub_power_w > 10.0 && (z - 0.92).abs() < 0.22 && (0.18..=0.40).contains(&r) {
                        // Local volumetric heating at balance-piston labyrinth teeth during rub!
                        q_vol += rub_power_w / 0.045;
                    }

                    let t_next = (t_c + h_dt * (lap_k + q_vol) / rho_cp).clamp(280.0, 1_150.0);
                    self.scratch_t[idx] = t_next;
                }
            }
            self.temp_k.copy_from_slice(&self.scratch_t);
        }

        // Compute derived engineering fields: Fourier q = -k ∇T, thermo-elastic + centrifugal + hoop σ_vM,
        // steam pressure P(z), Mach(z), wetness, and Norton creep log-rate!
        let mut peak_vm = 0.0_f64;
        let mut peak_casing_hoop = 0.0_f64;
        let mut peak_q = 0.0_f64;
        let wall_m = design.clamped().casing_wall_mm * 1.0e-3;
        let omega = omega_rad_s.abs();

        for ir in 0..GRID_NR {
            let r = (ir as f64 + 0.5) * dr;
            for iz in 0..GRID_NZ {
                let idx = ir * GRID_NZ + iz;
                let reg = self.region[idx];
                let t_c = self.temp_k[idx];
                let z = GRID_Z_MIN_M + (iz as f64 + 0.5) * dz;

                // Central finite-difference temperature gradient -> exact Fourier heat flux q = -k(T) ∇T [kW/m²]
                let dt_dr = if ir == 0 {
                    (self.temp_k[(ir + 1) * GRID_NZ + iz] - t_c) / dr
                } else if ir + 1 == GRID_NR {
                    (t_c - self.temp_k[(ir - 1) * GRID_NZ + iz]) / dr
                } else {
                    (self.temp_k[(ir + 1) * GRID_NZ + iz] - self.temp_k[(ir - 1) * GRID_NZ + iz])
                        / (2.0 * dr)
                };
                let dt_dz = if iz == 0 {
                    (self.temp_k[ir * GRID_NZ + (iz + 1)] - t_c) / dz
                } else if iz + 1 == GRID_NZ {
                    (t_c - self.temp_k[ir * GRID_NZ + (iz - 1)]) / dz
                } else {
                    (self.temp_k[ir * GRID_NZ + (iz + 1)] - self.temp_k[ir * GRID_NZ + (iz - 1)])
                        / (2.0 * dz)
                };

                let k_val = if reg == CellRegion::SteamOrAir {
                    0.06
                } else {
                    lib.hp_casing_steel.conductivity_at(t_c)
                };
                let qr_kw = -k_val * dt_dr * 1.0e-3;
                let qz_kw = -k_val * dt_dz * 1.0e-3;
                self.q_r_kw_m2[idx] = qr_kw;
                self.q_z_kw_m2[idx] = qz_kw;
                if reg != CellRegion::SteamOrAir {
                    peak_q = peak_q.max(qr_kw.hypot(qz_kw));
                }

                // Flowpath pressure, Mach, and wetness
                self.p_bar[idx] = z_p_bar[iz];
                self.mach[idx] = z_mach[iz];
                self.wetness[idx] = z_wet[iz];

                // Coupled von Mises stress [MPa]:
                // 1) Rotating disc centrifugal stress σ_θ,cf = ((3+ν)/8) ρ ω² R_o² (exact solid rotating disc!)
                // 2) Blade root centrifugal + bending from 17-stage meanline
                // 3) Thick-walled Lamé pressure hoop stress σ_hoop = P_in * R_i / wall_t
                // 4) Local radial/axial thermal gradient stress σ_th = E α |∇T| L_char / (2(1-ν))
                let grad_t = dt_dr.hypot(dt_dz);
                let vm_mpa = match reg {
                    CellRegion::RotorForging => {
                        let r_outer = rotor_outer_radius(z).max(0.24);
                        let nu = lib.rotor_steel.poisson_ratio;
                        let cf_pa = ((3.0 + nu) / 8.0)
                            * lib.rotor_steel.rho_kg_m3
                            * omega
                            * omega
                            * (r_outer * r_outer - (1.0 + 3.0 * nu) / (3.0 + nu) * r * r).max(0.0);
                        let e_r = lib.rotor_steel.young_modulus_at(t_c);
                        let a_r = lib.rotor_steel.mean_alpha_at(t_c);
                        let th_pa = e_r * a_r * (grad_t * 0.14) / (2.0 * (1.0 - nu));
                        (cf_pa * 1.0e-6).hypot(th_pa * 1.0e-6)
                    }
                    CellRegion::BladesAnnulus => {
                        // Look up closest stage root/span centrifugal + bending stress!
                        let mut best_vm = 90.0_f64;
                        let mut best_dz = 99.0_f64;
                        for st in &meanline.stages {
                            let dzz = (st.z_m - z).abs();
                            if dzz < best_dz {
                                best_dz = dzz;
                                // Stress peaks at blade hub root and tapers to 0 at free tip
                                let span_f = ((r - st.r_hub_m) / (st.r_tip_m - st.r_hub_m).max(0.04))
                                    .clamp(0.0, 1.0);
                                best_vm = st.root_von_mises_mpa * (1.0 - 0.72 * span_f);
                            }
                        }
                        best_vm
                    }
                    CellRegion::InnerCarrier | CellRegion::OuterCasing => {
                        let p_local_pa = if z < 2.15 {
                            p_chest_bar * 1.0e5
                        } else {
                            z_p_bar[iz] * 1.0e5
                        };
                        let r_inner = steam_outer_radius(z);
                        // Thick-walled Lamé pressure hoop stress at inner bore
                        let hoop_pa = p_local_pa * (r_inner / wall_m.max(0.04));
                        let e_c = lib.hp_casing_steel.young_modulus_at(t_c);
                        let a_c = lib.hp_casing_steel.mean_alpha_at(t_c);
                        let nu = lib.hp_casing_steel.poisson_ratio;
                        let th_pa = e_c * a_c * (grad_t * wall_m) / (2.0 * (1.0 - nu));
                        let vm = (hoop_pa * 1.0e-6).hypot(th_pa * 1.0e-6);
                        if reg == CellRegion::OuterCasing && z < 2.2 {
                            peak_casing_hoop = peak_casing_hoop.max(vm);
                        }
                        vm
                    }
                    CellRegion::ExhaustHood => {
                        // External atmospheric vacuum crushing differential (~0.9 bar across 28 mm shell)
                        let dp_vac = (1.013e5 - z_p_bar[iz] * 1.0e5).max(0.0);
                        let hoop = dp_vac * (1.18 / 0.028) * 1.0e-6;
                        let th = 0.08 * grad_t * 0.03;
                        hoop.hypot(th)
                    }
                    CellRegion::BearingPedestal => 28.0 + 0.35 * (t_c - 320.0).abs(),
                    CellRegion::SteamOrAir => 0.0,
                };

                self.sigma_vm_mpa[idx] = vm_mpa;
                if reg != CellRegion::SteamOrAir {
                    peak_vm = peak_vm.max(vm_mpa);
                }

                let creep = if reg == CellRegion::RotorForging
                    || reg == CellRegion::BladesAnnulus
                    || reg == CellRegion::InnerCarrier
                    || reg == CellRegion::OuterCasing
                {
                    lib.rotor_creep_rate_per_s(vm_mpa * 1.0e6, t_c)
                } else {
                    0.0
                };
                self.creep_log[idx] = if creep > 1.0e-18 {
                    creep.log10().clamp(-18.0, -4.0)
                } else {
                    -18.0
                };
            }
        }

        self.peak_metal_stress_mpa = peak_vm;
        self.casing_hoop_stress_mpa = peak_casing_hoop;
        self.peak_heat_flux_kw_m2 = peak_q;
        self.pack_export();
    }

    fn pack_export(&mut self) {
        for idx in 0..GRID_CELLS {
            let b = idx * GRID_STRIDE;
            self.export_buf[b] = self.temp_k[idx] as f32;
            self.export_buf[b + 1] = self.sigma_vm_mpa[idx] as f32;
            self.export_buf[b + 2] = self.p_bar[idx] as f32;
            self.export_buf[b + 3] = self.mach[idx] as f32;
            self.export_buf[b + 4] = self.q_r_kw_m2[idx] as f32;
            self.export_buf[b + 5] = self.q_z_kw_m2[idx] as f32;
            self.export_buf[b + 6] = self.wetness[idx] as f32;
            self.export_buf[b + 7] = self.creep_log[idx] as f32;
        }
    }

    /// Zero-copy slice of the packed `64 × 24 × 8` continuum tensor for GPU DataTexture upload.
    #[must_use]
    pub fn export_flat(&self) -> &[f32] {
        &self.export_buf
    }

    /// Bilinear engineering probe at any axisymmetric `(z_m, r_m)` coordinate (for 3D hover inspection!).
    /// Returns `[temp_c, sigma_vm_mpa, p_bar, mach, q_mag_kw_m2, wetness_pct, creep_log, region_id]`.
    #[must_use]
    pub fn probe_at(&self, z_m: f64, r_m: f64) -> [f64; 8] {
        let u = ((z_m - GRID_Z_MIN_M) / (GRID_Z_MAX_M - GRID_Z_MIN_M)).clamp(0.0, 0.9999);
        let v = (r_m.abs() / GRID_R_MAX_M).clamp(0.0, 0.9999);
        let fz = u * (GRID_NZ - 1) as f64;
        let fr = v * (GRID_NR - 1) as f64;
        let iz0 = fz.floor() as usize;
        let ir0 = fr.floor() as usize;
        let iz1 = (iz0 + 1).min(GRID_NZ - 1);
        let ir1 = (ir0 + 1).min(GRID_NR - 1);
        let tz = fz - iz0 as f64;
        let tr = fr - ir0 as f64;

        let bilerp = |arr: &[f64; GRID_CELLS]| -> f64 {
            let a00 = arr[ir0 * GRID_NZ + iz0];
            let a10 = arr[ir0 * GRID_NZ + iz1];
            let a01 = arr[ir1 * GRID_NZ + iz0];
            let a11 = arr[ir1 * GRID_NZ + iz1];
            let a0 = a00 + tz * (a10 - a00);
            let a1 = a01 + tz * (a11 - a01);
            a0 + tr * (a1 - a0)
        };

        let q_r = bilerp(&self.q_r_kw_m2);
        let q_z = bilerp(&self.q_z_kw_m2);
        let reg = self.region[ir0 * GRID_NZ + iz0] as u8 as f64;
        [
            bilerp(&self.temp_k) - 273.15,
            bilerp(&self.sigma_vm_mpa),
            bilerp(&self.p_bar),
            bilerp(&self.mach),
            q_r.hypot(q_z),
            bilerp(&self.wetness) * 100.0,
            bilerp(&self.creep_log),
            reg,
        ]
    }
}

fn rotor_outer_radius(z: f64) -> f64 {
    if z < 0.255 {
        0.132
    } else if z < 0.350 {
        0.214 // Kingsbury thrust collar (geometry.rs z = 0.273..0.345)
    } else if z < 0.640 {
        0.128 // Front radial journal bearing land (FRONT_BEARING_Z_M = 0.50)
    } else if z < 0.880 {
        0.156 // Front gland seal
    } else if z < 1.185 {
        0.304 // Twin balance pistons (geometry.rs z = 0.880..1.185)
    } else if z < 1.385 {
        0.340 // Curtis impulse control wheel hub (z = 1.298)
    } else if z < 2.255 {
        0.275 // HP reaction drum
    } else if z < 3.380 {
        0.305 // IP1/IP2 reaction drums
    } else if z < 3.960 {
        0.410 // Giant LP Stage 16 integral forged disc (z = 3.775)
    } else if z < 4.210 {
        0.175 // Rear gland seal land
    } else {
        0.132 // Rear journal bearing shaft (REAR_BEARING_Z_M = 4.44, ends z = 4.78)
    }
}

fn steam_outer_radius(z: f64) -> f64 {
    if z < 2.02 {
        0.396
    } else if z < 3.35 {
        0.396 + 0.165 * ((z - 2.02) / 1.33)
    } else if z < 3.85 {
        0.561 + 0.385 * ((z - 3.35) / 0.50)
    } else {
        0.955 + 0.285 * ((z - 3.85) / 1.60).clamp(0.0, 1.0)
    }
}

fn classify_section_cell(z: f64, r: f64, wall_m: f64) -> CellRegion {
    // Bearing pedestals matched strictly to geometry.rs: front (z = 0.145..0.635), rear cone (z = 3.98..5.14)
    if (0.145..=0.635).contains(&z) && (0.135..=0.325).contains(&r) {
        return CellRegion::BearingPedestal;
    }
    if (3.980..=5.120).contains(&z) && (0.140..=0.415).contains(&r) {
        return CellRegion::BearingPedestal;
    }
    let r_rot = rotor_outer_radius(z);
    if (0.0..=4.78).contains(&z) && r <= r_rot {
        return CellRegion::RotorForging;
    }
    let r_ann_out = steam_outer_radius(z);
    // All 17 blade stages live strictly in z = 1.22..3.88 (Stage 0 z=1.298 -> Stage 16 z=3.775!)
    if (1.22..=3.88).contains(&z) && r > r_rot && r <= r_ann_out {
        return CellRegion::BladesAnnulus;
    }
    if (1.18..=3.64).contains(&z) && r > r_ann_out && r <= r_ann_out + 0.11 {
        return CellRegion::InnerCarrier;
    }
    if (0.62..=3.68).contains(&z) && r > r_ann_out + 0.11 && r <= r_ann_out + 0.11 + wall_m {
        return CellRegion::OuterCasing;
    }
    if (3.66..=5.45).contains(&z) && r > r_ann_out - 0.03 && r <= r_ann_out + 0.11 {
        return CellRegion::ExhaustHood;
    }
    CellRegion::SteamOrAir
}

fn initial_temp_for_cell(z: f64, _r: f64, reg: CellRegion, cold_start: bool) -> f64 {
    if cold_start {
        return match reg {
            CellRegion::SteamOrAir => 300.0,
            _ => 316.0,
        };
    }
    match reg {
        CellRegion::BearingPedestal => 349.0,
        CellRegion::ExhaustHood => 332.0,
        CellRegion::SteamOrAir => 300.0,
        _ => {
            // Smooth warm axial profile from ~455 °C at HP Curtis down to ~60 °C at LP exhaust
            let u = ((z - 1.10) / 2.85).clamp(0.0, 1.0);
            728.0 - 385.0 * u.powf(0.85)
        }
    }
}
