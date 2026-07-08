//! Topological derivatives: the asymptotic sensitivity of compliance
//! to inserting an infinitesimal (traction-free) hole — the principled
//! nucleation mechanism. For 2D linear elasticity the classical
//! result (Garreau–Guillaume–Masmoudi / Amstutz form, plane strain):
//!
//! `DT(x) = π(λ+2μ)/(2μ(λ+μ)) · [4μ σ:ε + (λ−μ) tr(σ)·tr(ε)]`
//!
//! evaluated at the pre-hole state — positive for compliance (removing
//! material makes a loaded structure softer). CONSTANTS ARE
//! NUMERICALLY GATED: the battery punches a real hole where the
//! derivative points and checks the measured compliance change against
//! `DT·(hole area)` within a documented first-order band, so a wrong
//! sign or scale cannot ship silently.

use crate::gridsdf::GridSdf;
use std::fmt::Write as _;

/// The compliance topological derivative per unit hole AREA at a
/// point, from the local stress/strain state (Voigt: xx, yy, xy with
/// tensor shear).
#[must_use]
pub fn topological_derivative(lambda: f64, mu: f64, sigma: [f64; 3], eps: [f64; 3]) -> f64 {
    let se = sigma[0] * eps[0] + sigma[1] * eps[1] + 2.0 * sigma[2] * eps[2];
    let tr_s = sigma[0] + sigma[1];
    let tr_e = eps[0] + eps[1];
    (lambda + 2.0 * mu) / (2.0 * mu * (lambda + mu)) * (4.0 * mu * se + (lambda - mu) * tr_s * tr_e)
        / 2.0
}

/// One nucleation event (ledger row).
#[derive(Debug, Clone)]
pub struct NucleationEvent {
    /// Hole center.
    pub center: [f64; 2],
    /// Hole radius.
    pub radius: f64,
    /// The topological derivative at the center.
    pub dt_value: f64,
    /// Predicted Lagrangian gain `(ℓ − DT)·πρ²` (> 0 fired).
    pub predicted_gain: f64,
}

impl NucleationEvent {
    /// Ledger-style JSON row.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"center\":[{:.4},{:.4}],\"radius\":{:.4},\"dt\":{:.4e},\
             \"predicted_gain\":{:.4e}}}",
            self.center[0], self.center[1], self.radius, self.dt_value, self.predicted_gain
        );
        s
    }
}

/// Punch holes where the augmented Lagrangian improves: candidate
/// nodes are well inside the material (|φ| > 2ρ) AND at least `margin`
/// from every box edge (clamps and loads are not nucleation targets),
/// the criterion is
/// `gain = (ℓ − DT)·πρ² > 0`, winners are picked greedily best-first
/// with a spacing of `3ρ`, capped at `max_holes`. Each hole updates
/// φ ← max(φ, ρ − |x − c|); the caller redistances afterwards.
#[must_use]
pub fn nucleate(
    phi: &mut GridSdf,
    dt_field: &[f64],
    ell: f64,
    radius: f64,
    margin: f64,
    max_holes: usize,
) -> Vec<NucleationEvent> {
    let n = phi.n();
    let stride = n + 1;
    assert_eq!(dt_field.len(), stride * stride, "nodal DT field");
    let area = std::f64::consts::PI * radius * radius;
    let mut candidates: Vec<(usize, f64)> = (0..stride * stride)
        .filter(|&k| {
            let (i, j) = (k % stride, k / stride);
            let p = phi.pos(i, j);
            let inside_margin =
                p[0] > margin && p[0] < 1.0 - margin && p[1] > margin && p[1] < 1.0 - margin;
            inside_margin && phi.node(i, j) < -2.0 * radius && (ell - dt_field[k]) > 0.0
        })
        .map(|k| (k, (ell - dt_field[k]) * area))
        .collect();
    candidates.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .expect("finite gains")
            .then(a.0.cmp(&b.0))
    });
    let mut events: Vec<NucleationEvent> = Vec::new();
    for (k, gain) in candidates {
        if events.len() >= max_holes {
            break;
        }
        let (i, j) = (k % stride, k / stride);
        let c = phi.pos(i, j);
        if events
            .iter()
            .any(|e| (e.center[0] - c[0]).hypot(e.center[1] - c[1]) < 3.0 * radius)
        {
            continue;
        }
        for jj in 0..=n {
            for ii in 0..=n {
                let p = phi.pos(ii, jj);
                let hole = radius - (p[0] - c[0]).hypot(p[1] - c[1]);
                let v = phi.node(ii, jj);
                *phi.node_mut(ii, jj) = v.max(hole);
            }
        }
        events.push(NucleationEvent {
            center: c,
            radius,
            dt_value: dt_field[k],
            predicted_gain: gain,
        });
    }
    events
}
