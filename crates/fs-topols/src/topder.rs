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

/// Where a hole may be punched: strictly inside the design box with room for
/// the whole hole, and clear of the regions a hole must never cut (supports
/// and loaded pads). Free box edges are NOT keep-outs.
#[derive(Debug, Clone, Copy)]
pub struct NucleationRegion {
    /// Holes must keep at least this distance from EVERY box edge (the hole
    /// itself must fit in the design domain), typically `ρ + h`.
    pub containment: f64,
    /// Axis-aligned keep-out boxes `[x0, y0, x1, y1]` for hole CENTERS
    /// (e.g. a clamped strip and a load pad, each widened by `ρ + 2h`).
    pub keep_out: [[f64; 4]; 2],
}

impl NucleationRegion {
    fn admits(&self, p: [f64; 2]) -> bool {
        let c = self.containment;
        let inside = p[0] > c && p[0] < 1.0 - c && p[1] > c && p[1] < 1.0 - c;
        inside
            && self.keep_out.iter().all(|b| !(p[0] >= b[0] && p[0] <= b[2] && p[1] >= b[1] && p[1] <= b[3]))
    }
}

/// Punch holes where exchanging material for boundary material lowers the
/// objective (q61wp.16).
///
/// With the material area held by the scheduled projection, a hole of area
/// `A` at `x` is followed by restoring `A` along the interface. To first
/// order the compliance changes by `(DT(x) − Λ)·A`, where `Λ` is the
/// interface EXCHANGE RATE: the compliance cost per unit area of boundary
/// material, `σ:ε` on the interface (`2w` with `w = ½σ:ε`). A hole is
/// admitted only where `DT(x) < Λ`, ranked by `(Λ − DT)·πρ²`, greedily with
/// `3ρ` spacing, capped at `max_holes`. The candidate node must be deeper than
/// `ρ + h` in material, so the hole leaves a ligament. This replaces the
/// augmented-Lagrangian multiplier as the threshold: once the area is
/// projected, the multiplier no longer measures the boundary trade-off. The
/// old edge margin of six radii on ALL edges emptied the admissible window on
/// the canonical bracket, so holes never nucleated at all.
#[must_use]
pub fn nucleate_by_exchange(
    phi: &mut GridSdf,
    dt_field: &[f64],
    exchange_rate: f64,
    radius: f64,
    region: NucleationRegion,
    max_holes: usize,
) -> Vec<NucleationEvent> {
    let n = phi.n();
    let stride = n + 1;
    assert_eq!(dt_field.len(), stride * stride, "nodal DT field");
    if !(exchange_rate.is_finite() && exchange_rate > 0.0 && radius > 0.0) {
        return Vec::new();
    }
    let h = phi.h();
    let area = std::f64::consts::PI * radius * radius;
    let mut candidates: Vec<(usize, f64)> = (0..stride * stride)
        .filter(|&k| {
            let (i, j) = (k % stride, k / stride);
            region.admits(phi.pos(i, j))
                && phi.node(i, j) < -(radius + h)
                && dt_field[k].is_finite()
                && dt_field[k] < exchange_rate
        })
        .map(|k| (k, (exchange_rate - dt_field[k]) * area))
        .collect();
    candidates.sort_by(|a, b| b.1.partial_cmp(&a.1).expect("finite gains").then(a.0.cmp(&b.0)));
    let mut events: Vec<NucleationEvent> = Vec::new();
    for (k, gain) in candidates {
        if events.len() >= max_holes {
            break;
        }
        let (i, j) = (k % stride, k / stride);
        let c = phi.pos(i, j);
        if events.iter().any(|e| (e.center[0] - c[0]).hypot(e.center[1] - c[1]) < 3.0 * radius) {
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
        events.push(NucleationEvent { center: c, radius, dt_value: dt_field[k], predicted_gain: gain });
    }
    events
}

#[cfg(test)]
mod exchange_tests {
    use super::{NucleationRegion, nucleate_by_exchange};
    use crate::gridsdf::GridSdf;

    /// Solid square (φ = −1 everywhere) on a 32x32 lattice.
    fn solid() -> GridSdf {
        GridSdf::from_fn(32, &|_, _| -1.0)
    }

    fn region(radius: f64, h: f64) -> NucleationRegion {
        let guard = radius + 2.0 * h;
        NucleationRegion {
            containment: radius + h,
            keep_out: [[0.0, 0.0, guard, 1.0], [1.0 - 2.0 * h - guard, 0.375 - guard, 1.0, 0.625 + guard]],
        }
    }

    #[test]
    fn holes_go_only_where_dt_beats_the_exchange_rate_and_never_in_keep_outs() {
        let mut phi = solid();
        let (n, h) = (phi.n(), phi.h());
        let rho = 1.5 * h;
        let stride = n + 1;
        // DT is low (attractive) along the clamped strip, near the free top
        // edge and inside the load pad; high everywhere else.
        let mut dt = vec![10.0; stride * stride];
        for j in 0..=n {
            for i in 0..=n {
                let p = phi.pos(i, j);
                if p[0] < 0.05 || (p[1] > 0.85 && (p[0] - 0.5).abs() < 0.02) || (p[0] > 0.95 && (p[1] - 0.5).abs() < 0.05) {
                    dt[i + j * stride] = 1.0;
                }
            }
        }
        let events = nucleate_by_exchange(&mut phi, &dt, 5.0, rho, region(rho, h), 4);
        assert!(!events.is_empty(), "the free-edge band must admit a hole");
        for e in &events {
            assert!(e.center[0] > rho + 2.0 * h, "never in the clamped strip: {e:?}");
            assert!(e.center[1] > 0.85, "only the free-top candidates beat the rate: {e:?}");
            assert!(e.dt_value < 5.0 && e.predicted_gain > 0.0);
        }
        for (a, b) in events.iter().zip(events.iter().skip(1)) {
            assert!((a.center[0] - b.center[0]).hypot(a.center[1] - b.center[1]) >= 3.0 * rho);
        }
        // A rate below every DT punches nothing and leaves φ untouched.
        let mut untouched = solid();
        assert!(nucleate_by_exchange(&mut untouched, &dt, 0.5, rho, region(rho, h), 4).is_empty());
        assert!(untouched.nodes().iter().all(|v| *v == -1.0));
    }
}
