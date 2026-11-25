//! Row-by-row 17-stage IAPWS-IF97 meanline aerothermodynamics, Euler velocity triangles,
//! Ainley-Mathieson/tip-clearance/Baumann loss breakdown, and blade-root mechanical stress
//! & Norton creep integrity for the model turbine steam turbine.
//!
//! Covers every physical blade stage from the Stage 0 Curtis impulse control wheel through
//! HP (1..=5), IP1 (6..=9), IP2 (10..=13), and LP flared stages (14..=16) using exact
//! geometry radii, annulus areas, and blade counts from [`crate::geometry::BladeRowSpec`]
//! plus engineer-editable [`TurbineDesignParams`].

use core::f64::consts::PI;

use crate::geometry::BladeRowSpec;
use crate::materials::library;
use crate::steam::{eval_ph, eval_ps};

/// Number of turbine expansion stages in the model turbine flowpath (Stage 0..=16).
pub const NUM_STAGES: usize = 17;

/// Number of `f64` telemetry fields exported per stage in [`StageMeanlineReport::to_flat_vec`].
pub const STAGE_ROW_STRIDE: usize = 16;

/// Engineer-editable physical machine design parameters for the turbine model.
///
/// Mutating these parameters directly alters the 17-stage meanline aerothermodynamics,
/// tip-leakage efficiency, blade centrifugal root stress, 4-DOF rotordynamics & oil-whip
/// stability margin, casing hoop/thermal stress, AND live L2 MORPH 3D CAD blade geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TurbineDesignParams {
    /// Radial blade-tip / shroud-fin clearance [mm] (nominal 0.75 mm; range 0.35..1.80 mm).
    pub tip_clearance_mm: f64,
    /// Target HP/IP reaction stage degree of reaction `R` (nominal 0.48; range 0.15..0.65).
    pub hp_reaction: f64,
    /// LP Stage 16 last-stage free-standing blade height [m] (nominal 0.527 m; range 0.40..0.66 m).
    pub lp_last_blade_h_m: f64,
    /// Front & rear journal bearing radial clearance [µm] (nominal 220 µm; range 150..300 µm).
    pub bearing_clearance_um: f64,
    /// Tilting-pad bearing geometric preload factor `m_p` (nominal 0.32; range 0.0..0.60).
    pub bearing_preload: f64,
    /// Outer cast CrMoV HP casing wall thickness [mm] (nominal 95 mm; range 65..145 mm).
    pub casing_wall_mm: f64,
}

impl Default for TurbineDesignParams {
    fn default() -> Self {
        Self {
            tip_clearance_mm: 0.75,
            hp_reaction: 0.48,
            lp_last_blade_h_m: 0.527, // exact r_tip 0.935 - r_hub 0.408 from geometry.rs!
            bearing_clearance_um: 220.0,
            bearing_preload: 0.32,
            casing_wall_mm: 95.0,
        }
    }
}

impl TurbineDesignParams {
    /// Clamp every design variable strictly into its admitted physical engineering envelope.
    #[must_use]
    pub fn clamped(self) -> Self {
        Self {
            tip_clearance_mm: self.tip_clearance_mm.clamp(0.35, 1.80),
            hp_reaction: self.hp_reaction.clamp(0.15, 0.65),
            lp_last_blade_h_m: self.lp_last_blade_h_m.clamp(0.40, 0.66),
            bearing_clearance_um: self.bearing_clearance_um.clamp(150.0, 300.0),
            bearing_preload: self.bearing_preload.clamp(0.0, 0.60),
            casing_wall_mm: self.casing_wall_mm.clamp(65.0, 145.0),
        }
    }

    /// Fraction of the HP section's flow that bypasses the blading through the shaft-gland and
    /// diaphragm labyrinths. A fitted value. It does not depend on the blade tip clearance: tip
    /// leakage is already a loss of each row in the meanline (`tip_loss`), and counting it here
    /// as well would charge the clearance twice.
    #[must_use]
    pub fn hp_leakage_fraction(&self) -> f64 {
        0.018
    }

    /// Fraction of the IP/LP section's flow that bypasses the blading through its labyrinths
    /// (fitted; independent of the blade tip clearance, as for the HP section).
    #[must_use]
    pub fn iplp_leakage_fraction(&self) -> f64 {
        0.014
    }

    /// Thrust-balance factor: higher reaction raises rotor pressure drop and net axial thrust.
    #[must_use]
    pub fn thrust_reaction_factor(&self) -> f64 {
        (0.62 + 0.79 * self.hp_reaction).clamp(0.70, 1.25)
    }
}

/// Complete aerodynamic & mechanical operating point of one stage (`0..=16`).
#[derive(Debug, Clone, Copy, Default)]
pub struct StageState {
    pub stage_index: usize,
    /// Axial center of the rotor blade row [m].
    pub z_m: f64,
    /// Hub radius [m].
    pub r_hub_m: f64,
    /// Tip radius [m] (reflects live `lp_last_blade_h_m` on LP stages!).
    pub r_tip_m: f64,
    /// Stage exit static pressure [Pa].
    pub p_out_pa: f64,
    /// Stage exit static temperature [K].
    pub t_out_k: f64,
    /// Stage exit specific enthalpy [J/kg].
    pub h_out_j_kg: f64,
    /// Stage exit specific entropy [J/(kg·K)].
    pub s_out_j_kg_k: f64,
    /// Stage exit steam dryness quality `x ∈ [0, 1]`.
    pub quality_x: f64,
    /// Stage steam mass flow [kg/s].
    pub m_dot_kg_s: f64,
    /// Mean axial flow velocity `c_a` [m/s].
    pub c_axial_m_s: f64,
    /// Stator exit absolute tangential swirl velocity `c_{u1}` [m/s].
    pub c_swirl_in_m_s: f64,
    /// Rotor exit absolute tangential swirl velocity `c_{u2}` [m/s].
    pub c_swirl_out_m_s: f64,
    /// Peak relative/stator Mach number `Ma` in this stage.
    pub mach_max: f64,
    /// Stage degree of reaction `R ∈ [0, 1]`.
    pub reaction: f64,
    /// Stage total-to-total isentropic efficiency `η_stage`.
    pub eta_stage: f64,
    /// Stage mechanical shaft power [W].
    pub power_w: f64,
    /// Stage aerodynamic axial thrust on rotor disc [N].
    pub axial_thrust_n: f64,
    /// Blade root centrifugal stress `σ_cf` [MPa].
    pub root_centrifugal_mpa: f64,
    /// Blade root aerodynamic steam bending stress `σ_bend` [MPa].
    pub root_bending_mpa: f64,
    /// Combined blade root von Mises stress `σ_vM` [MPa].
    pub root_von_mises_mpa: f64,
    /// Stage local Norton creep strain rate [1/s].
    pub creep_rate_per_s: f64,
}

/// Complete 17-stage meanline solution across the entire flowpath.
#[derive(Debug, Clone)]
pub struct StageMeanlineReport {
    pub stages: [StageState; NUM_STAGES],
    pub total_power_w: f64,
    pub total_axial_thrust_n: f64,
    pub max_blade_stress_mpa: f64,
    pub max_mach: f64,
    pub min_quality_x: f64,
    pub lp_leaving_loss_kw: f64,
    pub min_goodman_sf: f64,
}

impl Default for StageMeanlineReport {
    fn default() -> Self {
        Self {
            stages: [StageState::default(); NUM_STAGES],
            total_power_w: 0.0,
            total_axial_thrust_n: 0.0,
            max_blade_stress_mpa: 0.0,
            max_mach: 0.0,
            min_quality_x: 1.0,
            lp_leaving_loss_kw: 0.0,
            min_goodman_sf: 2.5,
        }
    }
}

/// Isentropic efficiencies of the plant's two stage sections, read off the row-by-row
/// meanline. `None` where a section carries too little flow or expansion to define one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SectionEfficiencies {
    /// HP section (stages 0..=5), chest to extraction pressure: `(h_0 - h_5) / (h_0 - h_5s)`.
    pub hp: Option<f64>,
    /// IP/LP section (stages 6..=16), extraction to condenser pressure, total-to-static: the
    /// last stage's leaving kinetic energy is not recovered, so
    /// `(m (h_5 - h_16) - P_leave) / (m (h_5 - h_16s))`.
    pub iplp: Option<f64>,
}

impl StageMeanlineReport {
    /// Section isentropic efficiencies of this meanline solution, which started from the
    /// chest state `(p_chest_pa, h_chest_j_kg)` with mass flows `m_hp_kg_s`, `m_lp_kg_s`.
    ///
    /// A section needs at least 0.5 kg/s and 1 kJ/kg of isentropic drop for its efficiency
    /// to mean anything; below that the plant keeps the last value (its work is then ~0).
    #[must_use]
    pub fn section_efficiencies(
        &self,
        p_chest_pa: f64,
        h_chest_j_kg: f64,
        m_hp_kg_s: f64,
        m_lp_kg_s: f64,
    ) -> SectionEfficiencies {
        const MIN_FLOW_KG_S: f64 = 0.5;
        const MIN_DROP_J_KG: f64 = 1_000.0;
        let ok = |eta: f64| (eta > 0.0 && eta < 1.0).then_some(eta);

        let s0 = eval_ph(p_chest_pa, h_chest_j_kg).s_j_kg_k;
        let last_hp = &self.stages[5];
        let dh_is_hp = h_chest_j_kg - eval_ps(last_hp.p_out_pa, s0).h_j_kg;
        let hp = if m_hp_kg_s >= MIN_FLOW_KG_S && dh_is_hp >= MIN_DROP_J_KG {
            ok((h_chest_j_kg - last_hp.h_out_j_kg) / dh_is_hp)
        } else {
            None
        };

        let s5 = eval_ph(last_hp.p_out_pa, last_hp.h_out_j_kg).s_j_kg_k;
        let last = &self.stages[NUM_STAGES - 1];
        let dh_is_lp = last_hp.h_out_j_kg - eval_ps(last.p_out_pa, s5).h_j_kg;
        let iplp = if m_lp_kg_s >= MIN_FLOW_KG_S && dh_is_lp >= MIN_DROP_J_KG {
            let work = m_lp_kg_s * (last_hp.h_out_j_kg - last.h_out_j_kg)
                - self.lp_leaving_loss_kw * 1.0e3;
            ok(work / (m_lp_kg_s * dh_is_lp))
        } else {
            None
        };
        SectionEfficiencies { hp, iplp }
    }

    /// Flatten all 17 stages (`17 × 16 = 272` `f64`s) for zero-copy WASM telemetry & 3D vector fields.
    ///
    /// Per-stage layout (stride 16):
    /// `[z_m, r_hub_m, r_tip_m, p_out_bar, t_out_c, quality_x, mach_max, eta_pct,
    ///   power_mw, c_axial_m_s, c_swirl_in_m_s, c_swirl_out_m_s, reaction,
    ///   root_cf_mpa, root_vm_mpa, axial_thrust_kn]`
    #[must_use]
    pub fn to_flat_vec(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(NUM_STAGES * STAGE_ROW_STRIDE);
        for st in &self.stages {
            out.extend_from_slice(&[
                st.z_m,
                st.r_hub_m,
                st.r_tip_m,
                st.p_out_pa * 1.0e-5,
                st.t_out_k - 273.15,
                st.quality_x,
                st.mach_max,
                st.eta_stage * 100.0,
                st.power_w * 1.0e-6,
                st.c_axial_m_s,
                st.c_swirl_in_m_s,
                st.c_swirl_out_m_s,
                st.reaction,
                st.root_centrifugal_mpa,
                st.root_von_mises_mpa,
                st.axial_thrust_n * 1.0e-3,
            ]);
        }
        out
    }

    /// Smoothly interpolate axial steam flowpath `(P_bar, T_k, Ma, wetness, c_a, c_u)` at any axial `z_m`.
    #[must_use]
    pub fn sample_flowpath_at_z(&self, z_m: f64) -> (f64, f64, f64, f64, f64, f64) {
        if z_m <= self.stages[0].z_m {
            let s = &self.stages[0];
            return (
                s.p_out_pa * 1.0e-5,
                s.t_out_k,
                s.mach_max,
                1.0 - s.quality_x,
                s.c_axial_m_s,
                s.c_swirl_in_m_s,
            );
        }
        let last = &self.stages[NUM_STAGES - 1];
        if z_m >= last.z_m {
            // Diffuser decelerates axial Mach gently from Stage 16 exit into condenser hood
            let decay = ((z_m - last.z_m) / 1.5).clamp(0.0, 1.0);
            return (
                last.p_out_pa * 1.0e-5,
                last.t_out_k,
                last.mach_max * (1.0 - 0.55 * decay),
                1.0 - last.quality_x,
                last.c_axial_m_s * (1.0 - 0.45 * decay),
                last.c_swirl_out_m_s * (1.0 - 0.7 * decay),
            );
        }
        for i in 0..(NUM_STAGES - 1) {
            let a = &self.stages[i];
            let b = &self.stages[i + 1];
            if z_m >= a.z_m && z_m <= b.z_m {
                let t = ((z_m - a.z_m) / (b.z_m - a.z_m).max(1.0e-5)).clamp(0.0, 1.0);
                // Log-linear interpolation on pressure matches exponential turbine pressure drop
                let p_bar = (a.p_out_pa.max(500.0).ln() * (1.0 - t)
                    + b.p_out_pa.max(500.0).ln() * t)
                    .exp()
                    * 1.0e-5;
                let t_k = a.t_out_k + t * (b.t_out_k - a.t_out_k);
                let ma = a.mach_max + t * (b.mach_max - a.mach_max);
                let wet = (1.0 - a.quality_x) + t * ((1.0 - b.quality_x) - (1.0 - a.quality_x));
                let ca = a.c_axial_m_s + t * (b.c_axial_m_s - a.c_axial_m_s);
                let cu = a.c_swirl_in_m_s + t * (b.c_swirl_in_m_s - a.c_swirl_in_m_s);
                return (p_bar, t_k, ma, wet, ca, cu);
            }
        }
        (
            last.p_out_pa * 1.0e-5,
            last.t_out_k,
            last.mach_max,
            1.0 - last.quality_x,
            last.c_axial_m_s,
            last.c_swirl_out_m_s,
        )
    }
}

/// Solve all 17 stages row-by-row from exact IAPWS-IF97 inlet/extraction/condenser boundary states,
/// computing Euler velocity triangles, Ainley-Mathieson + tip-clearance + Baumann losses,
/// centrifugal + steam bending root stresses, and Norton creep rate.
#[must_use]
pub fn solve_17_stage_meanline(
    blade_rows: &[BladeRowSpec],
    design: &TurbineDesignParams,
    p_chest_pa: f64,
    h_chest_j_kg: f64,
    p_ext_pa: f64,
    p_cond_pa: f64,
    m_dot_hp_kg_s: f64,
    m_dot_lp_kg_s: f64,
    omega_rad_s: f64,
) -> StageMeanlineReport {
    let lib = library();
    let d = design.clamped();
    let mut report = StageMeanlineReport::default();
    report.min_quality_x = 1.0;
    report.min_goodman_sf = 10.0;

    // Extract the 17 rotor rows in axial order
    let mut rotor_rows = [BladeRowSpec {
        stage_index: 0,
        family_id: 0,
        is_rotor: true,
        z_center_m: 1.30,
        r_hub_m: 0.34,
        r_tip_m: 0.392,
        blade_count: 56,
        axial_chord_m: 0.044,
        has_shroud: true,
        r_mean_m: 0.366,
        annulus_area_m2: 0.119,
    }; NUM_STAGES];
    for row in blade_rows {
        if row.is_rotor && row.stage_index < NUM_STAGES {
            rotor_rows[row.stage_index] = *row;
        }
    }

    // Apply live LP last-stage blade height design parameter cleanly to Stage 15/16 annulus
    let dh_lp = d.lp_last_blade_h_m - 0.527;
    rotor_rows[15].r_tip_m = (rotor_rows[15].r_tip_m + 0.45 * dh_lp).clamp(0.62, 0.84);
    rotor_rows[15].r_mean_m = 0.5 * (rotor_rows[15].r_hub_m + rotor_rows[15].r_tip_m);
    rotor_rows[15].annulus_area_m2 =
        PI * (rotor_rows[15].r_tip_m.powi(2) - rotor_rows[15].r_hub_m.powi(2));
    rotor_rows[16].r_tip_m = (rotor_rows[16].r_hub_m + d.lp_last_blade_h_m).clamp(0.80, 1.08);
    rotor_rows[16].r_mean_m = 0.5 * (rotor_rows[16].r_hub_m + rotor_rows[16].r_tip_m);
    rotor_rows[16].annulus_area_m2 =
        PI * (rotor_rows[16].r_tip_m.powi(2) - rotor_rows[16].r_hub_m.powi(2));

    // Pressure-ratio distribution weights (Stage 0 Curtis carries ~2.2x enthalpy drop of a reaction stage;
    // Stages 1..=5 expand from Curtis exit to p_ext; Stages 6..=16 expand from p_ext to p_cond!)
    const HP_WEIGHTS: [f64; 6] = [1.85, 0.82, 0.84, 0.86, 0.88, 0.90];
    const IPLP_WEIGHTS: [f64; 11] = [
        0.72, 0.74, 0.76, 0.80, 0.84, 0.88, 0.94, 1.02, 1.15, 1.32, 1.65,
    ];
    let hp_w_sum: f64 = HP_WEIGHTS.iter().sum();
    let iplp_w_sum: f64 = IPLP_WEIGHTS.iter().sum();

    let p_ch = p_chest_pa.max(1_000.0);
    let p_ex = p_ext_pa.clamp(1_000.0, p_ch);
    let p_cn = p_cond_pa.clamp(800.0, p_ex);

    let ln_hp_drop = (p_ch / p_ex).max(1.0).ln();
    let ln_iplp_drop = (p_ex / p_cn).max(1.0).ln();

    let mut cur_st = eval_ph(p_ch, h_chest_j_kg);
    let mut hp_cum_w = 0.0_f64;
    let mut iplp_cum_w = 0.0_f64;
    let tip_clear_m = d.tip_clearance_mm * 1.0e-3;
    let omega = omega_rad_s.abs();

    for i in 0..NUM_STAGES {
        let rspec = rotor_rows[i];
        let is_hp = i < 6;
        let m_dot = if is_hp {
            m_dot_hp_kg_s.max(0.0)
        } else {
            m_dot_lp_kg_s.max(0.0)
        };

        let p_out = if is_hp {
            hp_cum_w += HP_WEIGHTS[i];
            (p_ch.ln() - ln_hp_drop * (hp_cum_w / hp_w_sum)).exp().max(p_ex)
        } else {
            iplp_cum_w += IPLP_WEIGHTS[i - 6];
            (p_ex.ln() - ln_iplp_drop * (iplp_cum_w / iplp_w_sum)).exp().max(p_cn)
        };

        let st_is = eval_ps(p_out, cur_st.s_j_kg_k);
        let dh_is = (cur_st.h_j_kg - st_is.h_j_kg).max(0.0);

        // Stage degree of reaction: Stage 0 Curtis impulse (~0.10), HP/IP uses design.hp_reaction,
        // LP flared stages progressively climb to ~0.58-0.65 free-vortex reaction at meanline.
        let reaction = if i == 0 {
            0.10
        } else if i < 14 {
            d.hp_reaction
        } else {
            (d.hp_reaction + 0.05 * (i - 13) as f64).clamp(0.42, 0.68)
        };

        let h_blade = (rspec.r_tip_m - rspec.r_hub_m).max(0.03);
        let u_mean = omega * rspec.r_mean_m;
        let c0 = (2.0 * dh_is).sqrt().max(5.0);
        let nu = u_mean / c0;
        // Effective blade-speed ratio matching at nominal speed (compact industrial HP drum runs
        // near nu ≈ 0.44 on Curtis and nu ≈ 0.52..0.68 across reaction/LP stages)
        let nu_opt = if i == 0 {
            0.44
        } else if i < 6 {
            0.52
        } else {
            0.62
        };

        // Ainley-Mathieson / Soderberg / reaction-profile / tip-clearance / Baumann loss model:
        // Reaction stages (i >= 1) suffer extra rotor turning/impulse profile loss far from R = 0.50.
        // The plant's two stage sections take their efficiencies from this loss model
        // (`StageMeanlineReport::section_efficiencies`).
        let base_eta = if i == 0 { 0.835 } else if i < 6 { 0.905 } else { 0.920 };
        let vel_factor = (1.0 - 0.48 * ((nu - nu_opt) / nu_opt).powi(2)).clamp(0.22, 1.0);
        let reaction_factor = if i == 0 {
            1.0
        } else {
            let r_opt = if i < 14 { 0.50 } else { 0.56 };
            (1.0 - 0.14 * (reaction - r_opt).powi(2)).clamp(0.84, 1.0)
        };
        let shroud_factor = if rspec.has_shroud { 0.48 } else { 0.95 };
        // Tip-clearance leakage penalty scales inversely with blade height h_blade
        let tip_loss = (shroud_factor * (tip_clear_m / h_blade) * 1.85).clamp(0.0, 0.18);
        // Secondary endwall loss scales with chord / height aspect ratio
        let sec_loss = (0.018 * (rspec.axial_chord_m / h_blade)).clamp(0.004, 0.05);
        // Baumann wetness penalty from real IAPWS-IF97 Region 4 moisture
        let y_mean = 0.5 * (cur_st.moisture() + st_is.moisture());
        let baumann = (1.0 - 0.85 * y_mean).clamp(0.40, 1.0);

        let eta_stage =
            (base_eta * vel_factor * reaction_factor * baumann - tip_loss - sec_loss).clamp(0.15, 0.94);
        let dh_actual = dh_is * eta_stage;
        let h_out = (cur_st.h_j_kg - dh_actual).clamp(150_000.0, 4_050_000.0);
        let st_out = eval_ph(p_out, h_out);

        // Continuity axial velocity c_a = m_dot * v / A_annulus
        let v_mid = 0.5 * (cur_st.v_m3_kg + st_out.v_m3_kg);
        let c_axial = (m_dot * v_mid / rspec.annulus_area_m2.max(0.02)).clamp(0.0, 420.0);

        // Euler turbomachinery equation: Δh_actual = u_mean * (c_u1 - c_u2)
        // With symmetric reaction R: (c_u1 + c_u2) / (2 u) = 1 - R  =>  c_u1, c_u2 uniquely determined!
        let delta_cu = if u_mean > 5.0 {
            (dh_actual / u_mean).min(950.0)
        } else {
            c0 * 0.75
        };
        let cu_mean = (1.0 - reaction) * u_mean;
        let c_u1 = cu_mean + 0.5 * delta_cu;
        let c_u2 = cu_mean - 0.5 * delta_cu;

        // Absolute stator exit velocity c_1 and relative rotor exit velocity w_2
        let c1 = c_axial.hypot(c_u1);
        let w2 = c_axial.hypot(c_u2 - u_mean);
        let a_mid = (0.5 * (cur_st.sound_speed_m_s + st_out.sound_speed_m_s)).max(150.0);
        let mach_stator = c1 / a_mid;
        let mach_rotor_rel = w2 / st_out.sound_speed_m_s.max(150.0);
        let mach_max = mach_stator.max(mach_rotor_rel);

        let power_w = m_dot * dh_actual;
        // Stage axial thrust = static pressure drop across rotor (R * ΔP) * A_annulus
        let dp_stage = (cur_st.p_pa - p_out).max(0.0);
        let axial_thrust_n = reaction * dp_stage * rspec.annulus_area_m2;

        // Mechanical blade-root stress (12% Cr stainless for HP/IP, Ti-6Al-4V grade for giant LP Stage 16):
        // Centrifugal stress at blade root: σ_cf = ρ * ω² * (r_tip² - r_hub²)/2 * taper_factor
        let (rho_blade, sy_mpa) = if i == 16 {
            (4_430.0_f64, 830.0_f64) // Ti-6Al-4V free-standing LP last-stage blade (low density, high yield!)
        } else {
            (lib.hp_blade_steel.rho_kg_m3, 680.0_f64) // 12% Cr martensitic turbine blade steel
        };
        let taper_factor = if i >= 14 { 0.54 } else { 0.84 };
        let sigma_cf_pa =
            rho_blade * omega * omega * 0.5 * (rspec.r_tip_m.powi(2) - rspec.r_hub_m.powi(2)) * taper_factor;

        // Aerodynamic cantilever bending stress at blade root from tangential & axial steam loads:
        // M_root ≈ (F_tang / N_blades) * (h_blade / 2), section modulus Z_root ≈ 0.045 * chord³
        let f_tang_blade = if u_mean > 5.0 {
            (power_w / u_mean) / f64::from(rspec.blade_count)
        } else {
            0.0
        };
        let f_ax_blade = axial_thrust_n / f64::from(rspec.blade_count);
        let f_res_blade = f_tang_blade.hypot(f_ax_blade);
        let m_root = f_res_blade * 0.52 * h_blade;
        let z_mod = (0.042 * rspec.axial_chord_m.powi(3)).max(1.0e-7);
        let sigma_bend_pa = m_root / z_mod;

        let sigma_cf_mpa = sigma_cf_pa * 1.0e-6;
        let sigma_bend_mpa = sigma_bend_pa * 1.0e-6;
        let sigma_vm_mpa = (sigma_cf_mpa + sigma_bend_mpa).hypot(0.55 * sigma_bend_mpa);
        let creep_rate = lib.rotor_creep_rate_per_s(sigma_vm_mpa * 1.0e6, st_out.t_k);

        // Temperature-derated yield strength for Goodman safety factor
        let temp_derate = (1.0 - 0.000_55 * (st_out.t_k - 300.0).max(0.0)).clamp(0.45, 1.0);
        let sf = (sy_mpa * temp_derate) / sigma_vm_mpa.max(1.0);

        report.stages[i] = StageState {
            stage_index: i,
            z_m: rspec.z_center_m,
            r_hub_m: rspec.r_hub_m,
            r_tip_m: rspec.r_tip_m,
            p_out_pa: p_out,
            t_out_k: st_out.t_k,
            h_out_j_kg: st_out.h_j_kg,
            s_out_j_kg_k: st_out.s_j_kg_k,
            quality_x: st_out.quality_x,
            m_dot_kg_s: m_dot,
            c_axial_m_s: c_axial,
            c_swirl_in_m_s: c_u1,
            c_swirl_out_m_s: c_u2,
            mach_max,
            reaction,
            eta_stage,
            power_w,
            axial_thrust_n,
            root_centrifugal_mpa: sigma_cf_mpa,
            root_bending_mpa: sigma_bend_mpa,
            root_von_mises_mpa: sigma_vm_mpa,
            creep_rate_per_s: creep_rate,
        };

        report.total_power_w += power_w;
        report.total_axial_thrust_n += axial_thrust_n;
        report.max_blade_stress_mpa = report.max_blade_stress_mpa.max(sigma_vm_mpa);
        report.max_mach = report.max_mach.max(mach_max);
        report.min_quality_x = report.min_quality_x.min(st_out.quality_x);
        report.min_goodman_sf = report.min_goodman_sf.min(sf);

        if i == NUM_STAGES - 1 {
            let c_exit = c_axial.hypot(c_u2);
            report.lp_leaving_loss_kw = 0.5 * m_dot * c_exit * c_exit * 1.0e-3;
        }

        cur_st = st_out;
    }

    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::build_turbine_geometry;
    use crate::steam::eval_pt;

    #[test]
    fn meanline_17_stages_computes_physical_expansion_mach_and_blade_stress() {
        let cat = crate::with_clock_free_cx(build_turbine_geometry);
        let live = eval_pt(11.0e6, 798.15);
        let rep = solve_17_stage_meanline(
            &cat.blade_rows,
            &TurbineDesignParams::default(),
            9.5e6,
            live.h_j_kg,
            1.28e6,
            8_500.0,
            52.0,
            42.0,
            100.0 * PI,
        );
        assert!(
            rep.total_power_w > 35.0e6 && rep.total_power_w < 68.0e6,
            "17-stage total power out of range: {} MW",
            rep.total_power_w * 1e-6
        );
        assert!(
            rep.max_mach > 0.55 && rep.max_mach < 1.8,
            "LP last-stage Mach out of range: {}",
            rep.max_mach
        );
        assert!(
            rep.max_blade_stress_mpa > 70.0 && rep.max_blade_stress_mpa < 720.0,
            "LP blade root stress out of range: {} MPa",
            rep.max_blade_stress_mpa
        );
        // Tighter tip clearance must raise HP efficiency & total power!
        let tight = TurbineDesignParams {
            tip_clearance_mm: 0.40,
            ..TurbineDesignParams::default()
        };
        let loose = TurbineDesignParams {
            tip_clearance_mm: 1.60,
            ..TurbineDesignParams::default()
        };
        let r_tight = solve_17_stage_meanline(
            &cat.blade_rows, &tight, 9.5e6, live.h_j_kg, 1.28e6, 8_500.0, 52.0, 42.0, 100.0 * PI,
        );
        let r_loose = solve_17_stage_meanline(
            &cat.blade_rows, &loose, 9.5e6, live.h_j_kg, 1.28e6, 8_500.0, 52.0, 42.0, 100.0 * PI,
        );
        assert!(
            r_tight.total_power_w > r_loose.total_power_w + 1.5e5,
            "tighter tip clearance must boost power: {} vs {}",
            r_tight.total_power_w,
            r_loose.total_power_w
        );
    }
}
