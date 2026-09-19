//! Energy-discrete integration of fs-material's EXISTING WoolFelt law.
//! This is an antiderivative/force-quadrature adapter, not another contact law.
//! Rate-independent crush/reversal memory is retained between strikes. This
//! adapter does not pretend that elastic subloops are a measured Prony fit.
use fs_material::{Uniaxial, WoolFelt};
use fs_math::det;

pub type State = <WoolFelt as Uniaxial>::State;

pub fn demonstration_law() -> Result<WoolFelt, String> {
    WoolFelt::new(4.0e5, 0.2, 2.5, 3.2, 0.25, 0.8).map_err(|e| e.to_string())
}

/// Recoverable compression energy density [J/m^3], on a committed branch.
/// The caller first updates the maximum for an accepted new loading point.
pub fn stored(law: &WoolFelt, e: f64, state: &State) -> f64 {
    let residual = law.eps_residual(state);
    if e <= residual || state.eps_max <= residual { return 0.0; }
    let span = state.eps_max - residual;
    state.sig_max * span / (law.q + 1.0)
        * det::pow((e - residual) / span, law.q + 1.0)
}

/// Integral of the trial stress, with the last COMMITTED branch held fixed.
fn primitive(law: &WoolFelt, e: f64, state: &State) -> f64 {
    if e <= state.eps_max { return stored(law, e, state); }
    let virgin = |x: f64| {
        if x <= 0.0 { 0.0 } else { law.envelope(x).0 * x / (law.p + 1.0) }
    };
    stored(law, state.eps_max, state) + virgin(e) - virgin(state.eps_max)
}

/// Held force and its derivative with respect to END overlap. Its work is
/// exactly the integral of the named trial stress over the strain increment,
/// including crossing the free gap and residual-crush point. No strain clamp.
pub fn average(law: &WoolFelt, state: &State, start_m: f64, end_m: f64,
    thickness_m: f64, area_m2: f64) -> (f64, f64) {
    let a = start_m / thickness_m;
    let b = end_m / thickness_m;
    let de = b - a;
    let (stress, slope) = if de.abs() < 1.0e-7 * (a.abs().max(b.abs()) + 1.0e-12) {
        let mid = 0.5 * (a + b);
        (law.stress(mid, state), 0.5 * law.tangent(mid, state))
    } else {
        let stress = (primitive(law, b, state) - primitive(law, a, state)) / de;
        (stress, (law.stress(b, state) - stress) / de)
    };
    (area_m2 * stress, area_m2 * slope / thickness_m)
}

/// Solve F = average(delta0, delta_free - compliance * F), retaining
/// a bracket and accepting only a force/overlap residual in SI units.
/// Returns a force; no material state is committed during trial iterations.
pub fn solve(law: &WoolFelt, state: &State, start_m: f64, free_m: f64,
    compliance_m_n: f64, thickness_m: f64, area_m2: f64) -> Result<f64, &'static str> {
    let mut hi = free_m.min(thickness_m * law.eps_densify);
    let residual = |end: f64| {
        let (f, df) = average(law, state, start_m, end, thickness_m, area_m2);
        (end + compliance_m_n * f - free_m, 1.0 + compliance_m_n * df, f)
    };
    let (rhi, _, fhi) = residual(hi);
    if rhi < -1.0e-12 { return Err("felt exceeded its admitted densification strain"); }
    if rhi.abs() <= 1.0e-13 { return Ok(fhi); }
    let mut lo = start_m.min(free_m).min(0.0) - thickness_m;
    for _ in 0..16 {
        if residual(lo).0 <= 0.0 { break; }
        lo = 2.0 * lo - thickness_m;
    }
    if residual(lo).0 > 0.0 { return Err("felt contact could not bracket its end overlap"); }
    let mut end = 0.5 * (lo + hi);
    for _ in 0..48 {
        let (r, derivative, force) = residual(end);
        if !r.is_finite() || !derivative.is_finite() || !force.is_finite() || force < 0.0 {
            return Err("nonfinite or tensile discrete felt force");
        }
        if r.abs() <= 1.0e-13 { return Ok(force); }
        if r > 0.0 { hi = end; } else { lo = end; }
        let newton = end - r / derivative;
        end = if newton > lo && newton < hi { newton } else { 0.5 * (lo + hi) };
    }
    Err("felt contact exhausted its safeguarded solve budget")
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn virgin_loop_dissipates_and_committed_crush_is_not_reset() {
        let law = demonstration_law().unwrap();
        let mut state = law.initial_state();
        let mut e0 = 0.0;
        let mut work = 0.0;
        let mut previous_loss = 0.0;
        for e in (1..=100).map(|i| i as f64 * 0.003)
            .chain((0..100).rev().map(|i| i as f64 * 0.003)) {
            let (stress, _) = average(&law, &state, e0, e, 1.0, 1.0);
            work += stress * (e - e0);
            state = law.update_state(e, &state);
            let loss = work - stored(&law, e, &state);
            assert!(loss + 1e-8 >= previous_loss);
            previous_loss = loss;
            e0 = e;
        }
        assert!(work > 0.0);
        assert!(state.eps_max > 0.29);
        assert_eq!(law.stress(law.eps_residual(&state), &state), 0.0);
    }
    #[test]
    fn implicit_contact_closes_force_residual_without_clipping() {
        let law = demonstration_law().unwrap();
        let state = law.initial_state();
        for free in [-0.001, 0.0, 0.0001, 0.002, 0.006] {
            let c = 2e-7;
            let f = solve(&law, &state, 0.0, free, c, 0.008, 1e-4).unwrap();
            let expected = average(&law, &state, 0.0, free - c * f, 0.008, 1e-4).0;
            assert!((f - expected).abs() < 1e-5);
            assert!(f >= 0.0);
        }
        assert!(solve(&law, &state, 0.0, 0.1, 1e-10, 0.008, 1e-4).is_err());
    }
}
