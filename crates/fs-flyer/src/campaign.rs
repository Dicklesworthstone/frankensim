//! V-08b1/V-08b2 campaign + A0/A1 historical-default selection (bead
//! wf-root-guzez.5.9, E4.3b3). The INDEPENDENT fs-wakeref referee
//! (V-08b1 — its own crate, its own kernel; NEVER the FOM judging
//! itself) supplies the truth series; the two candidates are:
//!
//!   A0 — quasi-steady: the instantaneous DC map of the A1 lane
//!        (zero-lag; today's fixed-control spine tier), and
//!   A1 — the shared-basis scheduled ROM (E4.3b2) marched in time.
//!
//! DECISION RULE (pre-declared here, recorded verbatim in the
//! receipt): the candidates are scored on TRANSIENT SHAPE against the
//! referee — each model's step/impulse response is normalized by its
//! OWN steady asymptote (the models legitimately differ in DC tier:
//! UVLM vs lifting-line closures), and the score is the RMS deviation
//! of that normalized shape from the referee's over the transient
//! window. Lower aggregate score wins. The E0.6 runtime-budget axis is
//! EXPLICIT NO-DATA until the browser microbench lands (recorded);
//! physics error alone decides until then. The selection gates E5.3b
//! (the historical mode) ONLY — never the fixed-control spine — and
//! the loser remains an admission-selectable identified mode.

use crate::Refusal;
use fs_blake3::hash_domain;
use fs_wakeref::{Fixture as RefFixture, RefereeCase, run_case, wright_geometry_v1};
use fs_wing::images::CertifiedGround;
use fs_wing::prescribedwake::frozen_grid_v1;
use fs_wing::rom::{A1Lti, assemble_a1_lti, wright_a1_layout_v1};
use fs_wing::romreduce::{project, reduce_shared};

/// Campaign fixtures (the overlap of the referee's registered set and
/// what both candidates express: canard step + impulse, free + ground).
pub const CAMPAIGN_FIXTURES: [(&str, bool); 4] = [
    ("step", false),
    ("step", true),
    ("impulse", false),
    ("impulse", true),
];

/// One candidate's score on one fixture.
#[derive(Clone, Debug, PartialEq)]
pub struct FixtureScore {
    /// Fixture name.
    pub fixture: &'static str,
    /// Ground case?
    pub ground: bool,
    /// Normalized-shape RMS deviation from the referee (canard-lift
    /// channel — the driven channel).
    pub shape_rms: f64,
    /// The model's own normalized shape at sample 3 (liveness: A1's
    /// lag shows here as < 1; A0 is exactly 1 — proves the comparison
    /// can discriminate regardless of who wins).
    pub early_shape: f64,
}

/// The full campaign + selection receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct SelectionReceipt {
    /// Schema id.
    pub schema: &'static str,
    /// The decision rule, verbatim.
    pub decision_rule: &'static str,
    /// A0 per-fixture scores (FULL results — the loser's too).
    pub a0: Vec<FixtureScore>,
    /// A1 per-fixture scores.
    pub a1: Vec<FixtureScore>,
    /// Aggregate scores (mean over fixtures).
    pub a0_aggregate: f64,
    /// A1 aggregate.
    pub a1_aggregate: f64,
    /// The winner ("A0" or "A1").
    pub winner: &'static str,
    /// The loser stays selectable (plan law; always true).
    pub loser_remains_selectable: bool,
    /// Runtime-budget axis status (E0.6).
    pub budget_axis: &'static str,
    /// Referee receipt digest this campaign judged against.
    pub referee_digest: String,
    /// Receipt digest.
    pub receipt_digest: String,
}

/// Receipt schema.
pub const SELECTION_SCHEMA: &str = "org.frankensim.wf.a1-selection-receipt.v1";

/// The verbatim decision rule (also stored in the receipt).
pub const DECISION_RULE: &str = "transient-shape-rms-v1: normalize each model's response by its \
     own steady asymptote; RMS deviation from the referee's normalized shape over the transient \
     window; lower mean over the fixture set wins; runtime budget axis NO-DATA until E0.6";

fn normalized_shape(series: &[f64]) -> Vec<f64> {
    // Scale: the steady asymptote for step-class responses; the PEAK
    // magnitude for impulse-class responses whose tail returns to ~0
    // (normalizing by a near-zero tail would manufacture divergence).
    let steady = *series.last().unwrap_or(&0.0);
    let peak = series.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let scale = if steady.abs() >= 0.05 * peak.max(1e-12) {
        steady
    } else if peak > 1e-12 {
        peak
    } else {
        return series.to_vec();
    };
    series.iter().map(|v| v / scale).collect()
}

fn shape_rms(model: &[f64], referee: &[f64]) -> f64 {
    let n = model.len().min(referee.len());
    let m = normalized_shape(&model[..n]);
    let r = normalized_shape(&referee[..n]);
    let sum: f64 = m.iter().zip(r.iter()).map(|(a, b)| (a - b) * (a - b)).sum();
    (sum / n as f64).sqrt()
}

/// Run the campaign and emit the selection receipt.
///
/// # Errors
/// Referee/FOM/reduction refusals pass through.
pub fn run_selection_campaign() -> Result<SelectionReceipt, Refusal> {
    let map_r = |e: fs_wakeref::Refusal| Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    };
    let map_w = |e: fs_wing::Refusal| Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    };
    let geometry = wright_geometry_v1();
    let referee_receipt = fs_wakeref::emit_v08b1_receipt(&geometry).map_err(map_r)?;
    // A1 lane: FOM at the campaign's operating points + shared reduction.
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let ground = CertifiedGround {
        z_m: 3.0,
        certificate_slope: 0.000606,
        certificate_rms_m: 0.801,
    };
    // Free-air point (h/b 10) and the ground point (h/b 0.2 — the
    // referee's flat-ground case class), pitch 0.05, neutral controls.
    let free_pt = *grid
        .points
        .iter()
        .find(|p| {
            p.h_over_b == 10.0
                && p.pitch_rad == 0.05
                && p.roll_rad == 0.0
                && p.canard_rad == 0.0
                && p.warp_rad == 0.0
                && p.convection == 1.0
        })
        .expect("registered free point");
    let gnd_pt = *grid
        .points
        .iter()
        .find(|p| {
            p.h_over_b == 0.2
                && p.pitch_rad == 0.05
                && p.roll_rad == 0.0
                && p.canard_rad == 0.0
                && p.warp_rad == 0.0
                && p.convection == 1.0
        })
        .expect("registered ground point");
    let v = 13.0;
    let rows = 120;
    let fom_free = assemble_a1_lti(&layout, &free_pt, &ground, v, rows).map_err(map_w)?;
    let fom_gnd = assemble_a1_lti(&layout, &gnd_pt, &ground, v, rows).map_err(map_w)?;
    let refs = [&fom_free, &fom_gnd];
    let red = reduce_shared(&refs).map_err(map_w)?;
    let dt = 1.0 / 120.0;
    let n_steps = 480usize;
    let mut a0_scores = Vec::new();
    let mut a1_scores = Vec::new();
    for &(fixture, ground_case) in &CAMPAIGN_FIXTURES {
        // Referee truth.
        let ref_fixture = match fixture {
            "step" => RefFixture::Step,
            _ => RefFixture::Impulse,
        };
        let ref_case = RefereeCase {
            fixture: ref_fixture,
            ground_z_m: if ground_case { Some(-2.4) } else { None },
            v_mps: v,
            alpha0_rad: 0.05,
            rho_kg_m3: 1.294,
            convection: 1.0,
            dt_s: dt,
            n_steps,
        };
        let truth = run_case(&geometry, &ref_case).map_err(map_r)?;
        // Candidate input: δ(t) matching the referee's fixture.
        let delta = |k: usize| -> f64 {
            match ref_fixture {
                RefFixture::Step => 0.1,
                _ => {
                    if k == 0 {
                        0.1
                    } else {
                        0.0
                    }
                }
            }
        };
        let fom = if ground_case { &fom_gnd } else { &fom_free };
        // A1: the reduced scheduled system at this point, marched.
        let a1_sys = project(fom, &red.basis, red.order);
        let mut x = vec![0.0f64; red.order];
        let mut a1_series = Vec::with_capacity(n_steps);
        for k in 0..n_steps {
            let u = [delta(k), 0.0];
            let mut y = 0.0; // canard-lift channel (output 1)
            y += a1_sys.d[2] * u[0] + a1_sys.d[3] * u[1];
            for j in 0..red.order {
                y += a1_sys.c[red.order + j] * x[j];
            }
            a1_series.push(y);
            let mut dx = vec![0.0f64; red.order];
            for i in 0..red.order {
                let mut s = a1_sys.b[i * 2] * u[0] + a1_sys.b[i * 2 + 1] * u[1];
                for j in 0..red.order {
                    s += a1_sys.a[i * red.order + j] * x[j];
                }
                dx[i] = s;
            }
            for i in 0..red.order {
                x[i] += dt * dx[i];
            }
        }
        // A0: the instantaneous DC map (zero lag — today's spine tier).
        let dc =
            A1Lti::dc_direct(&layout, &fom.point, &ground, v, rows, [0.1, 0.0]).map_err(map_w)?;
        let a0_series: Vec<f64> = (0..n_steps).map(|k| dc[1] * delta(k) / 0.1).collect();
        let early = |series: &[f64]| -> f64 {
            let shape = normalized_shape(series);
            shape.get(3).copied().unwrap_or(0.0)
        };
        a0_scores.push(FixtureScore {
            fixture: if fixture == "step" { "step" } else { "impulse" },
            ground: ground_case,
            shape_rms: shape_rms(&a0_series, &truth.canard_lift_n),
            early_shape: early(&a0_series),
        });
        a1_scores.push(FixtureScore {
            fixture: if fixture == "step" { "step" } else { "impulse" },
            ground: ground_case,
            shape_rms: shape_rms(&a1_series, &truth.canard_lift_n),
            early_shape: early(&a1_series),
        });
    }
    let agg = |s: &[FixtureScore]| s.iter().map(|f| f.shape_rms).sum::<f64>() / s.len() as f64;
    let a0_aggregate = agg(&a0_scores);
    let a1_aggregate = agg(&a1_scores);
    let winner = if a1_aggregate < a0_aggregate {
        "A1"
    } else {
        "A0"
    };
    let mut bytes = Vec::new();
    for s in a0_scores.iter().chain(a1_scores.iter()) {
        bytes.extend_from_slice(s.fixture.as_bytes());
        bytes.push(u8::from(s.ground));
        bytes.extend_from_slice(&s.shape_rms.to_bits().to_le_bytes());
        bytes.extend_from_slice(&s.early_shape.to_bits().to_le_bytes());
    }
    bytes.extend_from_slice(winner.as_bytes());
    bytes.extend_from_slice(referee_receipt.receipt_digest.as_bytes());
    let receipt_digest = hash_domain(SELECTION_SCHEMA, &bytes).to_hex();
    Ok(SelectionReceipt {
        schema: SELECTION_SCHEMA,
        decision_rule: DECISION_RULE,
        a0: a0_scores,
        a1: a1_scores,
        a0_aggregate,
        a1_aggregate,
        winner,
        loser_remains_selectable: true,
        budget_axis: "NO-DATA: E0.6 browser microbench not landed; physics error alone decides",
        referee_digest: referee_receipt.receipt_digest,
        receipt_digest,
    })
}
