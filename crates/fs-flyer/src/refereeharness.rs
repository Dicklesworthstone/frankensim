//! E10.1 referee-plane harness (bead wf-root-guzez.11.2): batch
//! RE-RUNS at pinned configs through the production solvers vs the
//! independent referees (E4.9a analytic + the E4.9b dense unsteady
//! reference), discrepancy receipts aggregated into one harness
//! receipt, and OPTIONAL correction tables under the plan §4.2
//! rules: calibration pins ONLY, holdout EVALUATED (its error is in
//! the table, never hidden), domain-bound, and bound to the
//! EvidenceRegistryId they were built against — stale or
//! out-of-domain application refuses.

use crate::aircraft::wright_rotor_v1;
use crate::referee::{
    DiscrepancyReceipt, referee_biplane_efficiency, referee_liftingline_cl, referee_momentum_thrust,
};
use crate::{Refusal, refuse};
use fs_blake3::hash_domain;
use fs_wing::nonlinear::{InfluenceOperator, StripRegime, StripSpec, solve_nonlinear};
use fs_wing::{SurfaceId, flat_surface};

/// Pinned harness air state (the Dec-17 trim class).
pub const HARNESS_RHO: f64 = 1.294;
/// Pinned speed [m/s].
pub const HARNESS_V: f64 = 13.86;
/// Pinned camber ratio.
pub const HARNESS_CAMBER: f64 = 0.05;

/// The pinned E4.9a alpha grid [rad].
pub const PINNED_ALPHAS: [f64; 3] = [0.03, 0.05, 0.07];

fn solve_mono_lift(alpha: f64) -> Result<f64, Refusal> {
    let panels = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).map_err(|e| {
        refuse(
            "harness-geometry",
            e.message.clone(),
            "the pinned monoplane",
        )
    })?;
    let strips: Vec<StripSpec> = (0..8)
        .map(|s| StripSpec {
            panel_indices: vec![s, 8 + s],
            chord_m: 1.981,
            twist_rad: 0.0,
        })
        .collect();
    let fs_v = [
        HARNESS_V * fs_math::det::cos(alpha),
        0.0,
        HARNESS_V * fs_math::det::sin(alpha),
    ];
    let closure = |_s: usize, a: f64| -> (f64, StripRegime) {
        (
            2.0 * core::f64::consts::PI * (a + 2.0 * HARNESS_CAMBER),
            StripRegime::Attached,
        )
    };
    let op = InfluenceOperator::build(&panels, fs_v, HARNESS_RHO).map_err(|e| {
        refuse(
            "harness-solve",
            e.message.clone(),
            "pinned config must solve",
        )
    })?;
    let sol = solve_nonlinear(
        &op,
        &panels,
        &strips,
        fs_v,
        HARNESS_RHO,
        &closure,
        None,
        None,
    )
    .map_err(|e| {
        refuse(
            "harness-solve",
            e.message.clone(),
            "pinned config must solve",
        )
    })?;
    Ok(sol.total_lift_n)
}

/// Batch re-run of the pinned E4.9a fixture set: one lifting-line
/// receipt per pinned alpha, plus the biplane-efficiency and
/// momentum-thrust cases. Every receipt is EXECUTED here, not read
/// from a log.
///
/// # Errors
/// Solver refusals surface as `harness-solve`; receipt validation
/// refusals pass through.
pub fn run_e49a_batch() -> Result<Vec<DiscrepancyReceipt>, Refusal> {
    let mut rows = Vec::new();
    let q = 0.5 * HARNESS_RHO * HARNESS_V * HARNESS_V;
    let s_ref = 12.29 * 1.981;
    let ar = 12.29 / 1.981;
    let ids: [&'static str; 3] = ["e101-mono-a03", "e101-mono-a05", "e101-mono-a07"];
    for (k, alpha) in PINNED_ALPHAS.iter().enumerate() {
        let production = solve_mono_lift(*alpha)?;
        let referee = q * s_ref * referee_liftingline_cl(*alpha, HARNESS_CAMBER, ar, 1.0);
        rows.push(DiscrepancyReceipt::new(
            ids[k],
            "wing_lift",
            "N",
            production,
            referee,
            "closed-form lifting line, no shared code with the panel solver",
            "formulation-band-0.15",
        )?);
    }
    // Biplane interference at the 1903 gap.
    let mono = solve_mono_lift(0.05)?;
    let mut p = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2)
        .map_err(|e| refuse("harness-geometry", e.message.clone(), "biplane"))?;
    p.extend(
        flat_surface(SurfaceId::WingUpper, 12.29, 1.981, 0.0, -1.89, 8, 2)
            .map_err(|e| refuse("harness-geometry", e.message.clone(), "biplane"))?,
    );
    let mut strips = Vec::new();
    for plane in 0..2 {
        let base = plane * 16;
        for s in 0..8 {
            strips.push(StripSpec {
                panel_indices: vec![base + s, base + 8 + s],
                chord_m: 1.981,
                twist_rad: 0.0,
            });
        }
    }
    let fs_v = [
        HARNESS_V * fs_math::det::cos(0.05),
        0.0,
        HARNESS_V * fs_math::det::sin(0.05),
    ];
    let closure = |_s: usize, a: f64| -> (f64, StripRegime) {
        (
            2.0 * core::f64::consts::PI * (a + 2.0 * HARNESS_CAMBER),
            StripRegime::Attached,
        )
    };
    let op = InfluenceOperator::build(&p, fs_v, HARNESS_RHO)
        .map_err(|e| refuse("harness-solve", e.message.clone(), "biplane"))?;
    let bi = solve_nonlinear(&op, &p, &strips, fs_v, HARNESS_RHO, &closure, None, None)
        .map_err(|e| refuse("harness-solve", e.message.clone(), "biplane"))?
        .total_lift_n;
    rows.push(DiscrepancyReceipt::new(
        "e101-biplane",
        "biplane_lift",
        "N",
        bi,
        2.0 * mono * referee_biplane_efficiency(1.89 / 12.29),
        "Prandtl/Munk biplane correction on the closed-form monoplane",
        "reported-only",
    )?);
    // BEMT vs momentum at matched power.
    let rotor = wright_rotor_v1();
    let omega = 52.0;
    let sol = fs_airscrew::bemt_solve(&rotor, HARNESS_RHO, HARNESS_V, omega)
        .map_err(|e| refuse("harness-solve", e.message.clone(), "bemt"))?;
    let area = core::f64::consts::PI * rotor.radius_m * rotor.radius_m;
    let ideal = referee_momentum_thrust(sol.torque_nm * omega, HARNESS_V, HARNESS_RHO, area)?;
    rows.push(DiscrepancyReceipt::new(
        "e101-prop-momentum",
        "thrust",
        "N",
        sol.thrust_n,
        ideal,
        "actuator-disk momentum bisection at matched power, no BEMT code",
        "below-ideal-within-0.5",
    )?);
    Ok(rows)
}

/// The aggregated harness receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct RefereeHarnessReceiptV1 {
    /// Schema id.
    pub schema: &'static str,
    /// Per-case rows.
    pub rows: Vec<DiscrepancyReceipt>,
    /// Worst |rel| across rows.
    pub worst_abs_rel: f64,
    /// E4.9b ingestion: the dense referee's Wagner-class starting
    /// ratio (canard lift at step 3 / steady) on the pinned Step case.
    pub wakeref_wagner_ratio: f64,
    /// The dense referee's series digest (its own identity).
    pub wakeref_series_digest: String,
    /// Digest over the whole receipt.
    pub receipt_digest: String,
}

/// Run the whole harness: the E4.9a batch plus the E4.9b dense
/// unsteady reference (fs-wakeref pinned Step case, free air).
///
/// # Errors
/// Batch refusals; `harness-wakeref` when the dense case refuses.
pub fn run_harness() -> Result<RefereeHarnessReceiptV1, Refusal> {
    let rows = run_e49a_batch()?;
    let case = fs_wakeref::RefereeCase {
        fixture: fs_wakeref::Fixture::Step,
        ground_z_m: None,
        v_mps: 13.0,
        alpha0_rad: 0.05,
        rho_kg_m3: HARNESS_RHO,
        convection: 1.0,
        dt_s: 1.0 / 120.0,
        n_steps: 240,
    };
    let series = fs_wakeref::run_case(&fs_wakeref::wright_geometry_v1(), &case)
        .map_err(|e| refuse("harness-wakeref", e.message.clone(), "the pinned Step case"))?;
    let steady = series.canard_lift_n[239];
    let wagner = series.canard_lift_n[3] / steady;
    let worst = rows
        .iter()
        .map(|r| r.rel_discrepancy.abs())
        .fold(0.0f64, f64::max);
    let mut b = Vec::new();
    for r in &rows {
        b.extend_from_slice(r.case_id.as_bytes());
        b.extend_from_slice(&r.rel_discrepancy.to_bits().to_le_bytes());
    }
    b.extend_from_slice(&wagner.to_bits().to_le_bytes());
    b.extend_from_slice(series.digest.as_bytes());
    let receipt_digest = hash_domain("org.frankensim.wf.referee-harness.v1", &b).to_hex();
    Ok(RefereeHarnessReceiptV1 {
        schema: "org.frankensim.wf.referee-harness.v1",
        rows,
        worst_abs_rel: worst,
        wakeref_wagner_ratio: wagner,
        wakeref_series_digest: series.digest,
        receipt_digest,
    })
}

/// Correction-knot cap.
pub const MAX_CORRECTION_KNOTS: usize = 32;

/// A §4.2 correction table: calibration pins only, holdout error
/// carried IN the table, domain-bound, registry-bound.
#[derive(Clone, Debug, PartialEq)]
pub struct CorrectionTable {
    /// Corrected quantity.
    pub quantity: &'static str,
    /// Domain lower bound (inclusive).
    pub domain_lo: f64,
    /// Domain upper bound (inclusive).
    pub domain_hi: f64,
    /// (x, referee/production) knots, ascending x.
    pub knots: Vec<(f64, f64)>,
    /// Case ids the calibration consumed (pins — never evidence).
    pub calibration_ids: Vec<&'static str>,
    /// Worst |corrected/referee − 1| on the HOLDOUT set (honest).
    pub holdout_rel_worst: f64,
    /// The registry the table was built against.
    pub built_against_registry_id: String,
}

/// Build a correction table from calibration receipts and EVALUATE
/// it on a disjoint holdout (mandatory — a table without a holdout
/// evaluation never exists).
///
/// # Errors
/// `correction-calibration-invalid` (0 knots or beyond the cap — AT
/// the cap admits; unordered x); `correction-holdout-missing`.
pub fn build_correction_table(
    quantity: &'static str,
    calibration: &[(f64, DiscrepancyReceipt)],
    holdout: &[(f64, DiscrepancyReceipt)],
    registry_id: &str,
) -> Result<CorrectionTable, Refusal> {
    if calibration.is_empty() || calibration.len() > MAX_CORRECTION_KNOTS {
        return Err(refuse(
            "correction-calibration-invalid",
            format!(
                "{} knots outside [1, {MAX_CORRECTION_KNOTS}]",
                calibration.len()
            ),
            "pin a bounded calibration set",
        ));
    }
    if calibration.windows(2).any(|w| w[1].0 <= w[0].0) {
        return Err(refuse(
            "correction-calibration-invalid",
            "calibration x not strictly ascending".into(),
            "sort and dedupe the pins",
        ));
    }
    if holdout.is_empty() {
        return Err(refuse(
            "correction-holdout-missing",
            "a correction table without a holdout evaluation".into(),
            "§4.2: holdout evaluated, always",
        ));
    }
    let knots: Vec<(f64, f64)> = calibration
        .iter()
        .map(|(x, r)| (*x, r.referee / r.production))
        .collect();
    let domain_lo = knots.first().expect("nonempty").0;
    let domain_hi = knots.last().expect("nonempty").0;
    let table = CorrectionTable {
        quantity,
        domain_lo,
        domain_hi,
        knots,
        calibration_ids: calibration.iter().map(|(_, r)| r.case_id).collect(),
        holdout_rel_worst: 0.0,
        built_against_registry_id: registry_id.to_string(),
    };
    let mut worst = 0.0f64;
    for (x, r) in holdout {
        let corrected = apply_correction(&table, *x, r.production, registry_id)?;
        worst = worst.max((corrected / r.referee - 1.0).abs());
    }
    Ok(CorrectionTable {
        holdout_rel_worst: worst,
        ..table
    })
}

/// Apply a correction. Refuses stale (wrong registry) and
/// out-of-domain applications — the DONE-WHEN hostile twins.
///
/// # Errors
/// `correction-stale`; `correction-out-of-domain` (AT the bounds
/// admits, beyond refuses).
pub fn apply_correction(
    table: &CorrectionTable,
    x: f64,
    production: f64,
    registry_id: &str,
) -> Result<f64, Refusal> {
    if registry_id != table.built_against_registry_id {
        return Err(refuse(
            "correction-stale",
            format!(
                "table built against {}, applied under {registry_id}",
                table.built_against_registry_id
            ),
            "rebuild the table under the current registry",
        ));
    }
    if !(x >= table.domain_lo && x <= table.domain_hi) {
        return Err(refuse(
            "correction-out-of-domain",
            format!("x {x} outside [{}, {}]", table.domain_lo, table.domain_hi),
            "corrections never extrapolate",
        ));
    }
    // Linear interpolation between bracketing knots.
    let mut factor = table.knots.first().expect("nonempty").1;
    for w in table.knots.windows(2) {
        let (x0, f0) = w[0];
        let (x1, f1) = w[1];
        if x >= x0 && x <= x1 {
            let t = if x1 > x0 { (x - x0) / (x1 - x0) } else { 0.0 };
            factor = f0 + t * (f1 - f0);
            break;
        }
    }
    Ok(production * factor)
}
