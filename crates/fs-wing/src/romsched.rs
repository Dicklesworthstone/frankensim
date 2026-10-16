//! Projected-operator scheduling + state continuation + the V-08b2
//! receipt (bead wf-root-guzez.5.8.3, E4.3b2-iii).
//!
//! Scheduling law: the reduced matrices are interpolated IN THE SHARED
//! BASIS (the only representation in which interpolation is meaningful
//! — the 5.8.2 battery's rotation falsifier is the proof), including
//! the ground axis: images entered through the operator BEFORE the
//! projection, so NO image logic exists here (the receipt carries the
//! structural proof — this module never touches `images::`).
//!
//! State continuation: a schedule switch keeps the reduced state
//! vector — same coordinates by construction — so the output is
//! continuous across switches while the steady levels genuinely move.
//!
//! Balanced truncation and Loewner are DIAGNOSTIC-ONLY per plan; this
//! receipt records them as explicit NO-DATA (never a silent blank —
//! the NO-DATA-vs-measured-zero doctrine).

use crate::romreduce::{ReducedLti, transfer_of};
use crate::{Refusal, refuse};
use fs_blake3::hash_domain;

/// One scheduled axis sample: parameter value + the reduced system in
/// the SHARED basis at that value.
#[derive(Clone, Debug, PartialEq)]
pub struct SchedulePoint {
    /// Axis parameter (e.g. convection ratio, or h/b on the ground axis).
    pub param: f64,
    /// The reduced system at this point (shared basis).
    pub sys: ReducedLti,
}

/// A 1-D scheduled ROM along one axis.
#[derive(Clone, Debug, PartialEq)]
pub struct ScheduledRom {
    /// Samples in strictly increasing parameter order.
    pub points: Vec<SchedulePoint>,
    /// Reduced order (uniform).
    pub order: usize,
}

impl ScheduledRom {
    /// Build from samples (validated: >= 2 points, strictly increasing,
    /// uniform order).
    ///
    /// # Errors
    /// `rom-schedule-invalid`.
    pub fn new(points: Vec<SchedulePoint>) -> Result<Self, Refusal> {
        if points.len() < 2 {
            return Err(refuse(
                "rom-schedule-invalid",
                format!("{} points", points.len()),
                "a schedule needs at least two axis samples",
            ));
        }
        let order = points[0].sys.order;
        for w in points.windows(2) {
            if !(w[1].param > w[0].param) {
                return Err(refuse(
                    "rom-schedule-invalid",
                    "parameters not strictly increasing".into(),
                    "sort the axis samples",
                ));
            }
        }
        if points.iter().any(|p| p.sys.order != order) {
            return Err(refuse(
                "rom-schedule-invalid",
                "mixed reduced orders".into(),
                "one shared basis, one order",
            ));
        }
        Ok(ScheduledRom { points, order })
    }

    /// The interpolated reduced system at `param` (linear in the shared
    /// basis; clamping is REFUSED — extrapolation is not scheduling).
    ///
    /// # Errors
    /// `rom-schedule-out-of-domain` (param at the edges admits; one
    /// float past either edge refuses).
    pub fn at(&self, param: f64) -> Result<ReducedLti, Refusal> {
        let lo = self.points.first().expect("validated").param;
        let hi = self.points.last().expect("validated").param;
        if !(param.is_finite() && (lo..=hi).contains(&param)) {
            return Err(refuse(
                "rom-schedule-out-of-domain",
                format!("param {param} outside [{lo}, {hi}]"),
                "scheduling interpolates; it never extrapolates",
            ));
        }
        let idx = self
            .points
            .windows(2)
            .position(|w| param <= w[1].param)
            .expect("param in range");
        let a = &self.points[idx];
        let b = &self.points[idx + 1];
        let t = (param - a.param) / (b.param - a.param);
        let mix = |x: &[f64], y: &[f64]| -> Vec<f64> {
            x.iter()
                .zip(y.iter())
                .map(|(u, v)| (1.0 - t) * u + t * v)
                .collect()
        };
        let mut d = [0.0f64; 6];
        for (k, dk) in d.iter_mut().enumerate() {
            *dk = (1.0 - t) * a.sys.d[k] + t * b.sys.d[k];
        }
        Ok(ReducedLti {
            order: self.order,
            a: mix(&a.sys.a, &b.sys.a),
            b: mix(&a.sys.b, &b.sys.b),
            c: mix(&a.sys.c, &b.sys.c),
            d,
        })
    }
}

/// March the scheduled ROM with a per-step (input, param) path —
/// the state VECTOR CONTINUES across parameter changes (shared-basis
/// coordinates; the whole point of the 5.8.2 law).
///
/// # Errors
/// `rom-march-invalid` (dt/steps caps at cap AND cap+1); schedule
/// refusals pass through.
pub fn march_scheduled(
    rom: &ScheduledRom,
    path: &dyn Fn(usize) -> ([f64; 2], f64),
    dt_s: f64,
    n_steps: usize,
) -> Result<Vec<[f64; 3]>, Refusal> {
    if !(dt_s.is_finite() && dt_s > 0.0 && dt_s <= 0.01 && (1..=200_000).contains(&n_steps)) {
        return Err(refuse(
            "rom-march-invalid",
            format!("dt {dt_s}, steps {n_steps}"),
            "dt (0, 0.01]; steps [1, 200000]",
        ));
    }
    let r = rom.order;
    let mut x = vec![0.0f64; r];
    let mut out = Vec::with_capacity(n_steps);
    for k in 0..n_steps {
        let (u, param) = path(k);
        let sys = rom.at(param)?;
        let mut y = [0.0f64; 3];
        for (o, yo) in y.iter_mut().enumerate() {
            let mut s = sys.d[o * 2] * u[0] + sys.d[o * 2 + 1] * u[1];
            for j in 0..r {
                s += sys.c[o * r + j] * x[j];
            }
            *yo = s;
        }
        out.push(y);
        let mut dx = vec![0.0f64; r];
        for i in 0..r {
            let mut s = sys.b[i * 2] * u[0] + sys.b[i * 2 + 1] * u[1];
            for j in 0..r {
                s += sys.a[i * r + j] * x[j];
            }
            dx[i] = s;
        }
        for i in 0..r {
            x[i] += dt_s * dx[i];
        }
    }
    Ok(out)
}

/// Phase and group delay of one channel at ω (group delay via the
/// symmetric finite difference of unwrapped phase — declared).
///
/// # Errors
/// Solver refusals pass through.
pub fn phase_and_group_delay(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64; 6],
    n: usize,
    ch: usize,
    w: f64,
) -> Result<(f64, f64), Refusal> {
    let dw = 1e-4 * w.max(1e-3);
    let g0 = transfer_of(a, b, c, d, n, w - dw)?[ch];
    let g1 = transfer_of(a, b, c, d, n, w + dw)?[ch];
    let p0 = g0.1.atan2(g0.0);
    let p1 = g1.1.atan2(g1.0);
    let mut dp = p1 - p0;
    // Unwrap the single step.
    while dp > core::f64::consts::PI {
        dp -= core::f64::consts::TAU;
    }
    while dp < -core::f64::consts::PI {
        dp += core::f64::consts::TAU;
    }
    let gm = transfer_of(a, b, c, d, n, w)?[ch];
    Ok((gm.1.atan2(gm.0), -dp / (2.0 * dw)))
}

/// One receipt clause verdict.
#[derive(Clone, Debug, PartialEq)]
pub struct ClauseVerdict {
    /// Clause name (stable).
    pub clause: &'static str,
    /// Passed?
    pub passed: bool,
    /// The measured figure the verdict rests on.
    pub measure: f64,
}

/// The V-08b2 reduction receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct V08b2Receipt {
    /// Schema id.
    pub schema: &'static str,
    /// Winning reduced order.
    pub order: usize,
    /// Ladder summary (order, worst error, passed) triplets.
    pub ladder: Vec<(usize, f64, bool)>,
    /// Clause verdicts (every DONE-WHEN clause EXECUTED).
    pub clauses: Vec<ClauseVerdict>,
    /// Diagnostic-only fields (explicit NO-DATA, never silent blanks).
    pub balanced_truncation: &'static str,
    /// Loewner diagnostic status.
    pub loewner: &'static str,
    /// Receipt digest.
    pub receipt_digest: String,
}

/// Receipt schema id.
pub const V08B2_SCHEMA: &str = "org.frankensim.wf.v08b2-receipt.v1";

/// Assemble the receipt from clause verdicts + the 5.8.2 ladder.
#[must_use]
pub fn emit_v08b2_receipt(
    order: usize,
    ladder: Vec<(usize, f64, bool)>,
    clauses: Vec<ClauseVerdict>,
) -> V08b2Receipt {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(order as u64).to_le_bytes());
    for (o, e, p) in &ladder {
        bytes.extend_from_slice(&(*o as u64).to_le_bytes());
        bytes.extend_from_slice(&e.to_bits().to_le_bytes());
        bytes.push(u8::from(*p));
    }
    for c in &clauses {
        bytes.extend_from_slice(c.clause.as_bytes());
        bytes.push(u8::from(c.passed));
        bytes.extend_from_slice(&c.measure.to_bits().to_le_bytes());
    }
    let receipt_digest = hash_domain(V08B2_SCHEMA, &bytes).to_hex();
    V08b2Receipt {
        schema: V08B2_SCHEMA,
        order,
        ladder,
        clauses,
        balanced_truncation: "NO-DATA: diagnostic-only per plan; the method receipt cites rational-Krylov",
        loewner: "NO-DATA: diagnostic-only per plan; the method receipt cites rational-Krylov",
        receipt_digest,
    }
}
