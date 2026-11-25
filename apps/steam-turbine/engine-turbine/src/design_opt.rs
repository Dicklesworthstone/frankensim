//! L4 ASCENT Multi-Variable Engineering Design Optimization for the model turbine
//! powered directly by [`fs_dfo::cma::cmaes`] (Information-Geometric Natural-Gradient CMA-ES)
//! and [`fs_dfo::nelder_mead`] local simplex polish.
//!
//! Optimizes the 6 physical machine design variables ([`TurbineDesignParams`]) directly over:
//! - Full 17-stage IAPWS-IF97 meanline expansion & Euler velocity triangles ([`solve_17_stage_meanline`])
//! - Real Ocvirk short-bearing `(K_ij, C_ij)` + 4-DOF `fs_mbd::rotordynamics::RotorBearingModel`
//!   stability eigenvalues (`max Re(λ)`)
//! - Blade centrifugal + steam bending root von Mises stress & Goodman safety factor
//! - Thick-walled Lamé pressure hoop + thermal transient stress on the cast CrMoV HP casing
//! - Net electrical Rankine efficiency, CHP utilization, and Norton creep life.

use fs_dfo::cma::{CmaParams, cmaes};
use fs_dfo::nelder_mead;

use crate::geometry::TurbineGeometryCatalog;
use crate::materials::library;
use crate::rotordynamics::{assemble_4dof_rotor_model, evaluate_journal_bearing};
use crate::stage_meanline::{TurbineDesignParams, solve_17_stage_meanline};
use crate::steam::{eval_ph, eval_pt, sat_pair_at_p};

/// Engineering optimization objective selectable in the Design Studio.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DesignObjective {
    /// Maximize net electrical efficiency subject to blade stress, dryness, & rotor stability limits.
    MaxElectricEfficiency = 0,
    /// Maximize combined CHP heat+power utilization & output power under safety constraints.
    MaxChpUtilization = 1,
    /// Maximize 4-DOF rotordynamic stability margin (`-max Re(λ)`) + creep/stress margin with high η.
    MaxStabilityAndLife = 2,
}

impl DesignObjective {
    #[must_use]
    pub fn from_u32(id: u32) -> Self {
        match id {
            1 => Self::MaxChpUtilization,
            2 => Self::MaxStabilityAndLife,
            _ => Self::MaxElectricEfficiency,
        }
    }
}

/// Complete physical evaluation of one candidate design at the current operating envelope.
#[derive(Debug, Clone, Copy)]
pub struct DesignMetrics {
    pub elec_eff_pct: f64,
    pub chp_eff_pct: f64,
    pub elec_power_mw: f64,
    pub max_blade_stress_mpa: f64,
    pub casing_stress_mpa: f64,
    pub spectral_abscissa_per_s: f64,
    pub bearing_loss_kw: f64,
    pub lp_leaving_loss_kw: f64,
    pub min_quality_x: f64,
    pub min_goodman_sf: f64,
}

/// Result of one complete `fs-dfo` (`cmaes` + `nelder_mead`) design optimization study.
#[derive(Debug, Clone)]
pub struct DesignOptimizationReport {
    pub objective: DesignObjective,
    pub baseline_params: TurbineDesignParams,
    pub optimized_params: TurbineDesignParams,
    pub baseline: DesignMetrics,
    pub optimized: DesignMetrics,
    pub total_evals: usize,
    pub cma_generations: usize,
    /// Monotone best-objective trace sampled across the CMA-ES + Nelder-Mead run (up to 32 points).
    pub history_obj: Vec<f64>,
}

/// Map normalized `[0, 1]^6` vector `u` -> physical [`TurbineDesignParams`].
fn decode_params(u: &[f64]) -> TurbineDesignParams {
    let c = |i: usize| u.get(i).copied().unwrap_or(0.5).clamp(0.0, 1.0);
    TurbineDesignParams {
        tip_clearance_mm: 0.35 + c(0) * (1.80 - 0.35),
        hp_reaction: 0.15 + c(1) * (0.65 - 0.15),
        lp_last_blade_h_m: 0.40 + c(2) * (0.66 - 0.40),
        bearing_clearance_um: 150.0 + c(3) * (300.0 - 150.0),
        bearing_preload: 0.0 + c(4) * 0.60,
        casing_wall_mm: 65.0 + c(5) * (145.0 - 65.0),
    }
    .clamped()
}

/// Map physical [`TurbineDesignParams`] -> normalized `[0, 1]^6` vector.
fn encode_params(p: &TurbineDesignParams) -> [f64; 6] {
    let p = p.clamped();
    [
        (p.tip_clearance_mm - 0.35) / (1.80 - 0.35),
        (p.hp_reaction - 0.15) / (0.65 - 0.15),
        (p.lp_last_blade_h_m - 0.40) / (0.66 - 0.40),
        (p.bearing_clearance_um - 150.0) / (300.0 - 150.0),
        p.bearing_preload / 0.60,
        (p.casing_wall_mm - 65.0) / (145.0 - 65.0),
    ]
}

/// Evaluate exact coupled 17-stage IAPWS-IF97 + 4-DOF rotordynamics + structural stress metrics
/// for a candidate [`TurbineDesignParams`] at the live operating boundary conditions.
#[must_use]
pub fn evaluate_design_candidate(
    geom: &TurbineGeometryCatalog,
    params: &TurbineDesignParams,
    p_live_pa: f64,
    t_live_k: f64,
    p_chest_pa: f64,
    p_ext_pa: f64,
    p_cond_pa: f64,
    m_dot_hp_kg_s: f64,
    m_dot_lp_kg_s: f64,
    omega_rad_s: f64,
) -> DesignMetrics {
    let d = params.clamped();
    let lib = library();
    let live = eval_pt(p_live_pa, t_live_k);
    let meanline = solve_17_stage_meanline(
        &geom.blade_rows,
        &d,
        p_chest_pa,
        live.h_j_kg,
        p_ext_pa,
        p_cond_pa,
        m_dot_hp_kg_s,
        m_dot_lp_kg_s,
        omega_rad_s,
    );

    // 4-DOF gyroscopic bearing stability & shear losses using live bearing_clearance_um & bearing_preload
    let span = (geom.rotor_props.rear_bearing_rel_z_m - geom.rotor_props.front_bearing_rel_z_m)
        .abs()
        .max(1.0);
    let w_tot = geom.rotor_props.total_mass_kg * 9.806_65;
    let f_front = w_tot * geom.rotor_props.rear_bearing_rel_z_m.abs() / span;
    let f_rear = w_tot * geom.rotor_props.front_bearing_rel_z_m.abs() / span;

    // Higher tilting-pad preload lowers effective cross-coupling (exact same law & bearing geometry
    // as `JournalBearingComponent::with_design` in components.rs!) at modest shear drag cost.
    let eff_destab = (0.32 - d.bearing_preload).max(0.0) * 1.1;
    let c_brg = d.bearing_clearance_um * 1.0e-6;
    let front = evaluate_journal_bearing(0.128, 0.24, c_brg, f_front, omega_rad_s, 333.15, 1.0, eff_destab);
    let rear = evaluate_journal_bearing(0.132, 0.26, c_brg, f_rear, omega_rad_s, 335.15, 1.0, eff_destab);

    let alford_kxy = 2.0 * (meanline.total_power_w / omega_rad_s.max(20.0)) / (2.0 * 0.42 * 0.18);
    let model = assemble_4dof_rotor_model(&geom.rotor_props, &front, &rear, alford_kxy);
    let max_re = model
        .at_speed(omega_rad_s, || false)
        .ok()
        .and_then(|op| op.poles(|| false).ok())
        .map(|rep| {
            rep.poles_per_s
                .iter()
                .map(|p| p.re)
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .unwrap_or(-2.0);

    // Preload gently raises direct damping & shear drag (unified with JournalBearingComponent!)
    let preload_shear_factor = 1.0 + 0.35 * d.bearing_preload;
    let bearing_loss_w = (front.shear_loss_w + rear.shear_loss_w) * preload_shear_factor;
    // Tip-clearance bypass steam leakage reduces net working steam flow identically to components.rs
    let hp_leak_frac = d.hp_leakage_fraction();
    let iplp_leak_frac = d.iplp_leakage_fraction();
    let hp_power_w: f64 = meanline.stages[0..6].iter().map(|s| s.power_w).sum();
    let iplp_power_w: f64 = meanline.stages[6..17].iter().map(|s| s.power_w).sum();
    let net_shaft_w = hp_power_w * (1.0 - hp_leak_frac) + iplp_power_w * (1.0 - iplp_leak_frac);
    let gen_loss_w = 0.014 * net_shaft_w + 220_000.0;
    let p_elec_w = (net_shaft_w - bearing_loss_w - gen_loss_w - meanline.lp_leaving_loss_kw * 1.0e3)
        .max(0.0);

    // Exact IAPWS-IF97 condensate hotwell enthalpy & net Rankine / CHP efficiencies
    let h_f_cond = sat_pair_at_p(p_cond_pa.clamp(1_500.0, 4.0e5))
        .liquid()
        .specific_enthalpy()
        .value();
    let q_boiler_w = (m_dot_hp_kg_s.max(0.1) * (live.h_j_kg - h_f_cond)).max(1.0);
    let m_ext = (m_dot_hp_kg_s - m_dot_lp_kg_s).max(0.0);
    let st_ext = eval_ph(p_ext_pa, meanline.stages[5].h_out_j_kg);
    let q_ext_w = (m_ext * (st_ext.h_j_kg - h_f_cond)).max(0.0);

    // Thick-walled Lamé pressure hoop + startup thermal gradient stress on outer HP casing:
    // Thin wall -> high pressure hoop stress; overly thick wall -> high thermal gradient stress!
    let wall_m = d.casing_wall_mm * 1.0e-3;
    let hoop_mpa = (p_chest_pa * 0.44 / wall_m) * 1.0e-6;
    let e_c = lib.hp_casing_steel.young_modulus_at(700.0);
    let a_c = lib.hp_casing_steel.mean_alpha_at(700.0);
    let nu = lib.hp_casing_steel.poisson_ratio;
    // Thermal ΔT across casing wall during startup scales quadratically with wall thickness (Fo = αt/L²)
    let dt_wall = 28.0 * (d.casing_wall_mm / 95.0).powf(1.65);
    let thermal_mpa = (e_c * a_c * dt_wall / (2.0 * (1.0 - nu))) * 1.0e-6;
    let casing_stress_mpa = hoop_mpa.hypot(thermal_mpa);

    DesignMetrics {
        elec_eff_pct: (p_elec_w / q_boiler_w * 100.0).clamp(0.0, 65.0),
        chp_eff_pct: ((p_elec_w + q_ext_w) / q_boiler_w * 100.0).clamp(0.0, 96.0),
        elec_power_mw: p_elec_w * 1.0e-6,
        max_blade_stress_mpa: meanline.max_blade_stress_mpa,
        casing_stress_mpa,
        spectral_abscissa_per_s: max_re,
        bearing_loss_kw: bearing_loss_w * 1.0e-3,
        lp_leaving_loss_kw: meanline.lp_leaving_loss_kw,
        min_quality_x: meanline.min_quality_x,
        min_goodman_sf: meanline.min_goodman_sf,
    }
}

/// Run FrankenSim's L4 ASCENT [`fs_dfo::cma::cmaes`] followed by [`fs_dfo::nelder_mead`] polish
/// directly in native Rust or WebAssembly to optimize the 6 design parameters!
#[must_use]
pub fn optimize_turbine_design(
    geom: &TurbineGeometryCatalog,
    baseline_params: TurbineDesignParams,
    objective: DesignObjective,
    max_evals: usize,
    p_live_pa: f64,
    t_live_k: f64,
    p_chest_pa: f64,
    p_ext_pa: f64,
    p_cond_pa: f64,
    m_dot_hp_kg_s: f64,
    m_dot_lp_kg_s: f64,
    omega_rad_s: f64,
) -> DesignOptimizationReport {
    let baseline = evaluate_design_candidate(
        geom,
        &baseline_params,
        p_live_pa,
        t_live_k,
        p_chest_pa,
        p_ext_pa,
        p_cond_pa,
        m_dot_hp_kg_s,
        m_dot_lp_kg_s,
        omega_rad_s,
    );

    let budget = max_evals.clamp(36, 360);
    let cma_budget = (budget * 3) / 4;
    let nm_budget = (budget - cma_budget).max(16);

    let mut best_seen = f64::INFINITY;
    let mut history_obj: Vec<f64> = Vec::with_capacity(32);
    let stride = (budget / 28).max(1);
    let mut eval_idx = 0usize;

    let mut score_u = |u: &[f64]| -> f64 {
        // Soft box barrier keeping CMA-ES & Nelder-Mead cleanly inside admitted design bounds
        let mut box_pen = 0.0_f64;
        for &v in u {
            if v < 0.0 {
                box_pen += 80.0 * v * v;
            } else if v > 1.0 {
                box_pen += 80.0 * (v - 1.0) * (v - 1.0);
            }
        }
        let cand = decode_params(u);
        let m = evaluate_design_candidate(
            geom,
            &cand,
            p_live_pa,
            t_live_k,
            p_chest_pa,
            p_ext_pa,
            p_cond_pa,
            m_dot_hp_kg_s,
            m_dot_lp_kg_s,
            omega_rad_s,
        );

        // Engineering hard-margin soft penalties:
        // - LP blade root von Mises stress <= 520 MPa (Ti-6Al-4V fatigue & overspeed margin!)
        // - Outer HP casing Lamé + thermal stress <= 230 MPa (CrMoV ASME Sec VIII limit!)
        // - Rotordynamic stability Re(λ) <= -5.0 /s (safe margin away from oil whip!)
        // - Tip clearance rub risk below 0.48 mm (prevents cold-start differential expansion rub!)
        let pen_blade = (m.max_blade_stress_mpa - 520.0).max(0.0).powi(2) * 0.015;
        let pen_casing = (m.casing_stress_mpa - 230.0).max(0.0).powi(2) * 0.02;
        let pen_stab = (m.spectral_abscissa_per_s + 5.0).max(0.0).powi(2) * 0.85;
        let pen_rub = (0.48 - cand.tip_clearance_mm).max(0.0).powi(2) * 45.0;
        let pen_wet = (0.88 - m.min_quality_x).max(0.0).powi(2) * 200.0;
        let penalties = box_pen + pen_blade + pen_casing + pen_stab + pen_rub + pen_wet;

        let raw_cost = match objective {
            DesignObjective::MaxElectricEfficiency => {
                // Maximize net electrical efficiency & power, minimize leaving/bearing losses
                -m.elec_eff_pct - 0.12 * m.elec_power_mw + 0.002 * m.casing_stress_mpa
            }
            DesignObjective::MaxChpUtilization => {
                // Maximize CHP utilization + electric power output
                -m.chp_eff_pct - 0.35 * m.elec_eff_pct - 0.10 * m.elec_power_mw
            }
            DesignObjective::MaxStabilityAndLife => {
                // Maximize damping margin (-Re(λ)), minimize blade & casing stress, preserve high η
                0.45 * m.spectral_abscissa_per_s
                    + 0.022 * m.max_blade_stress_mpa
                    + 0.035 * m.casing_stress_mpa
                    - 0.75 * m.elec_eff_pct
            }
        };
        let total = raw_cost + penalties;
        if total < best_seen {
            best_seen = total;
        }
        if eval_idx % stride == 0 && history_obj.len() < 32 {
            history_obj.push(best_seen);
        }
        eval_idx += 1;
        total
    };

    let x0 = encode_params(&baseline_params);
    // 1. Global natural-gradient CMA-ES exploration (`fs-dfo::cma::cmaes`)
    let cma_params = CmaParams::standard(6, 0.22, cma_budget, f64::NEG_INFINITY);
    let cma_rep = cmaes(&mut score_u, &x0, &cma_params, 0x5354_3630_3044_464F);

    // 2. Deterministic Nelder-Mead local simplex polish (`fs-dfo::nelder_mead`)
    let (nm_x, nm_f, nm_evals) = nelder_mead(&mut score_u, &cma_rep.x_best, 0.05, nm_budget, f64::NEG_INFINITY);

    let best_u = if nm_f <= cma_rep.f_best {
        nm_x
    } else {
        cma_rep.x_best
    };
    let optimized_params = decode_params(&best_u);
    let optimized = evaluate_design_candidate(
        geom,
        &optimized_params,
        p_live_pa,
        t_live_k,
        p_chest_pa,
        p_ext_pa,
        p_cond_pa,
        m_dot_hp_kg_s,
        m_dot_lp_kg_s,
        omega_rad_s,
    );
    if history_obj.last().copied().unwrap_or(f64::INFINITY) > best_seen {
        history_obj.push(best_seen);
    }

    DesignOptimizationReport {
        objective,
        baseline_params,
        optimized_params,
        baseline,
        optimized,
        total_evals: cma_rep.evals + nm_evals,
        cma_generations: cma_rep.generations,
        history_obj,
    }
}

/// Flatten [`DesignOptimizationReport`] for zero-copy WASM transfer to the Design Studio UI.
///
/// Layout (28 fixed floats followed by convergence trace):
/// `[0..6]: optimized_params (tip_mm, hp_react, lp_h_m, brg_c_um, brg_preload, wall_mm)`
/// `[6..16]: baseline metrics (elec_eff, chp_eff, mw, blade_mpa, casing_mpa, re_s, brg_kw, leave_kw, min_x, sf)`
/// `[16..26]: optimized metrics (same 10 fields)`
/// `[26]: total_evals, [27]: cma_generations, [28..]: history_obj`
#[must_use]
pub fn pack_opt_report(rep: &DesignOptimizationReport) -> Vec<f64> {
    let p = &rep.optimized_params;
    let b = &rep.baseline;
    let o = &rep.optimized;
    let mut out = Vec::with_capacity(28 + rep.history_obj.len());
    out.extend_from_slice(&[
        p.tip_clearance_mm,
        p.hp_reaction,
        p.lp_last_blade_h_m,
        p.bearing_clearance_um,
        p.bearing_preload,
        p.casing_wall_mm,
        b.elec_eff_pct,
        b.chp_eff_pct,
        b.elec_power_mw,
        b.max_blade_stress_mpa,
        b.casing_stress_mpa,
        b.spectral_abscissa_per_s,
        b.bearing_loss_kw,
        b.lp_leaving_loss_kw,
        b.min_quality_x,
        b.min_goodman_sf,
        o.elec_eff_pct,
        o.chp_eff_pct,
        o.elec_power_mw,
        o.max_blade_stress_mpa,
        o.casing_stress_mpa,
        o.spectral_abscissa_per_s,
        o.bearing_loss_kw,
        o.lp_leaving_loss_kw,
        o.min_quality_x,
        o.min_goodman_sf,
        rep.total_evals as f64,
        rep.cma_generations as f64,
    ]);
    out.extend_from_slice(&rep.history_obj);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;
    use crate::geometry::build_turbine_geometry;

    #[test]
    fn fs_dfo_cmaes_and_nelder_mead_improve_turbine_design_efficiency() {
        let cat = crate::with_clock_free_cx(build_turbine_geometry);
        // Start from a deliberately suboptimal loose-clearance / short-LP-blade baseline
        let subopt = TurbineDesignParams {
            tip_clearance_mm: 1.55,
            hp_reaction: 0.22,
            lp_last_blade_h_m: 0.43,
            bearing_clearance_um: 280.0,
            bearing_preload: 0.08,
            casing_wall_mm: 135.0,
        };
        let rep = optimize_turbine_design(
            &cat,
            subopt,
            DesignObjective::MaxElectricEfficiency,
            120,
            11.0e6,
            798.15,
            9.5e6,
            1.28e6,
            8_500.0,
            52.0,
            42.0,
            100.0 * PI,
        );
        assert!(rep.total_evals >= 60);
        assert!(
            rep.optimized.elec_eff_pct > rep.baseline.elec_eff_pct + 0.8,
            "fs-dfo must visibly improve net electric efficiency: {:.2}% -> {:.2}%",
            rep.baseline.elec_eff_pct,
            rep.optimized.elec_eff_pct
        );
        assert!(
            rep.optimized.max_blade_stress_mpa <= 540.0,
            "optimized LP blade root stress must respect structural bound: {:.1} MPa",
            rep.optimized.max_blade_stress_mpa
        );
    }
}
