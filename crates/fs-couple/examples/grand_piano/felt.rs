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

/// Read-only operands of the force quadrature, for an opted-in failed-solve trace.
/// Strains are dimensionless; the two primitive values are in Pa (J/m^3).
pub(super) fn diagnostic_primitives(law: &WoolFelt, state: &State,
    start_m: f64, end_m: f64, thickness_m: f64) -> [f64; 4] {
    let a = start_m / thickness_m; let b = end_m / thickness_m;
    [a, b, primitive(law, a, state), primitive(law, b, state)]
}

fn midpoint_interval(a: f64, b: f64) -> bool {
    (b - a).abs() < 1e-7 * (a.abs().max(b.abs()) + 1e-12)
}

/// Actual quadrature and legacy threshold, evaluated only for an opted-in trace.
/// Reuse the numerical selector so the diagnostic cannot invent a stable branch.
pub(super) fn diagnostic_quadrature(law: &WoolFelt, state: &State,
    start_m: f64, end_m: f64, thickness_m: f64) -> (&'static str, bool) {
    let a = start_m / thickness_m;
    let b = end_m / thickness_m;
    let legacy = midpoint_interval(a, b);
    let mode = if same_branch_average(law, state, a, b).is_some() {
        "stable_same_branch"
    } else if legacy { "midpoint" } else { "primitive" };
    (mode, legacy)
}

/// Cancellation-free moments of one power-law branch, sigma(e) = C (e-r)^n.
/// With m=(a+b)/2-r and t=(b-a)/(2m), the mean divided by C*m^n is
/// sum_even binom(n,k)*t^k/(k+1). Its END derivative is n*C*m^(n-1)/2
/// times sum_k binom(n-1,k)*t^k/(k+1 for even k, k+2 for odd k).
/// The latter is the weighted tangent integral, not a stress subtraction.
fn same_branch_average(law: &WoolFelt, state: &State, a: f64, b: f64)
    -> Option<(f64, f64)> {
    // Retain the coincident and branch-boundary conventions of the material law.
    if a == b { return None; }
    let residual = law.eps_residual(state);
    let (origin, scale, amplitude, exponent) = if a >= state.eps_max
        && b >= state.eps_max && a > 0.0 && b > 0.0 {
        (0.0, law.eps_ref, law.f_ref, law.p)
    } else if a > residual && b > residual && a <= state.eps_max && b < state.eps_max {
        // Starting at the maximum and unloading is still one branch. Preserve
        // the material's existing END-at-maximum stress convention separately.
        (residual, state.eps_max - residual, state.sig_max, law.q)
    } else { return None; };
    let mid = 0.5 * (a - origin) + 0.5 * (b - origin);
    let t = 0.5 * (b - a) / mid;
    // Every binomial-term ratio is at most rho <= 1/8. After degree 24,
    // each unweighted tail is <= (1/8)^25/(1-1/8) < 3.1e-23.
    // Outside this domain retain the existing primitive/crossing evaluation.
    // Subnormal coordinates can lose the halved width used in t; leave those
    // on the original path rather than claiming the series-domain bound holds.
    if mid < f64::MIN_POSITIVE { return None; }
    let (mean, tangent) = power_law_moments(exponent, t)?;
    let sigma = amplitude * det::pow(mid / scale, exponent);
    // Do not derive this from sigma/mid: sigma can underflow while its
    // tangent is still representable on a very small positive interval.
    let derivative = amplitude * exponent * det::pow(mid / scale, exponent - 1.0) / scale;
    Some((sigma * mean, 0.5 * derivative * tangent))
}

/// Dimensionless mean and weighted END-tangent moments for a power law on
/// a symmetric interval. Shared with the live adapter, whose power anchors
/// use native powf instead of det::pow. None retains the distant fallback.
#[inline]
pub(crate) fn power_law_moments(exponent: f64, t: f64) -> Option<(f64, f64)> {
    if !(t.abs() * exponent.max(1.0) <= 0.125) { return None; }
    let (mut mean, mut tangent) = (1.0, 1.0);
    let (mut mean_term, mut tangent_term) = (1.0, 1.0);
    for k in 1_u32..=24 {
        let degree = f64::from(k);
        mean_term *= (exponent + 1.0 - degree) * t / degree;
        tangent_term *= (exponent - degree) * t / degree;
        if k.is_multiple_of(2) {
            mean += mean_term / (degree + 1.0);
            tangent += tangent_term / (degree + 1.0);
        } else {
            tangent += tangent_term / (degree + 2.0);
        }
        // Remaining raw terms sum to <= current_term/7, so this early stop
        // is below EPSILON/56; it is not a contact-solver tolerance change.
        if mean_term.abs().max(tangent_term.abs()) < f64::EPSILON / 8.0 { break; }
    }
    Some((mean, tangent))
}

/// Held force and its derivative with respect to END overlap. Its work is
/// exactly the integral of the named trial stress over the strain increment,
/// including crossing the free gap and residual-crush point. No strain clamp.
pub fn average(law: &WoolFelt, state: &State, start_m: f64, end_m: f64,
    thickness_m: f64, area_m2: f64) -> (f64, f64) {
    let a = start_m / thickness_m;
    let b = end_m / thickness_m;
    let de = b - a;
    let (stress, slope) = if let Some(pair) = same_branch_average(law, state, a, b) {
        pair
    } else if midpoint_interval(a, b) {
        let mid = 0.5 * (a + b);
        (law.stress(mid, state), 0.5 * law.tangent(mid, state))
    } else {
        let stress = (primitive(law, b, state) - primitive(law, a, state)) / de;
        (stress, (law.stress(b, state) - stress) / de)
    };
    (area_m2 * stress, area_m2 * slope / thickness_m)
}

/// RT-0425/M2AN Kelvin-Voigt hammer term: R_H * d(e^p)/dt. The discrete
/// difference is passive because e^p is monotone, and the conservative force
/// above remains the exact power-law energy difference. The unilateral branch
/// drops a tensile release force; its unrecovered potential is counted as loss.
pub fn average_with_rate(law: &WoolFelt, state: &State, start_m: f64, end_m: f64,
    thickness_m: f64, area_m2: f64, rate_n_s_m_p: f64, dt_s: f64) -> (f64, f64) {
    let (elastic, slope) = average(law, state, start_m, end_m, thickness_m, area_m2);
    if rate_n_s_m_p == 0.0 { return (elastic, slope); }
    let power = |e: f64| if e <= 0.0 { 0.0 } else { det::pow(e, law.p) };
    let viscous = rate_n_s_m_p * (power(end_m) - power(start_m)) / dt_s;
    let total = elastic + viscous;
    if total <= 0.0 { return (0.0, 0.0); }
    let tangent = if end_m > 0.0 {
        rate_n_s_m_p * law.p * det::pow(end_m, law.p - 1.0) / dt_s
    } else { 0.0 };
    (total, slope + tangent)
}

/// Solve F = average(delta0, delta_free - compliance * F), retaining
/// a bracket and accepting only a force/overlap residual in SI units.
/// Returns a force; no material state is committed during trial iterations.
pub fn solve(law: &WoolFelt, state: &State, start_m: f64, free_m: f64,
    compliance_m_n: f64, thickness_m: f64, area_m2: f64) -> Result<f64, &'static str> {
    solve_with_rate(law, state, start_m, free_m, compliance_m_n,
        thickness_m, area_m2, 0.0, 1.0)
}

pub fn solve_with_rate(law: &WoolFelt, state: &State, start_m: f64, free_m: f64,
    compliance_m_n: f64, thickness_m: f64, area_m2: f64,
    rate_n_s_m_p: f64, dt_s: f64) -> Result<f64, &'static str> {
    if !rate_n_s_m_p.is_finite() || rate_n_s_m_p < 0.0 || !dt_s.is_finite() || dt_s <= 0.0 {
        return Err("hammer rate coefficient and time step must be finite positive");
    }
    let mut hi = free_m.min(thickness_m * law.eps_densify);
    let residual = |end: f64| {
        let (f, df) = average_with_rate(law, state, start_m, end,
            thickness_m, area_m2, rate_n_s_m_p, dt_s);
        (end + compliance_m_n * f - free_m, 1.0 + compliance_m_n * df, f)
    };
    // Stop at the residual's own rounding floor, not at a fixed overlap: with the
    // R_H term the tangent is R_H p e^(p-1)/dt, so 1e-13 m of overlap could still
    // be 1e-4 N of force, more than the contact and stretching iterations allow.
    let floor = |end: f64, force: f64|
        16.0 * f64::EPSILON * (end.abs() + (compliance_m_n * force).abs() + free_m.abs());
    let (rhi, _, fhi) = residual(hi);
    if rhi < -1.0e-12 { return Err("felt exceeded its admitted densification strain"); }
    if rhi.abs() <= floor(hi, fhi) { return Ok(fhi); }
    let mut lo = start_m.min(free_m).min(0.0) - thickness_m;
    for _ in 0..16 {
        if residual(lo).0 <= 0.0 { break; }
        lo = 2.0 * lo - thickness_m;
    }
    if residual(lo).0 > 0.0 { return Err("felt contact could not bracket its end overlap"); }
    let mut end = 0.5 * (lo + hi);
    for _ in 0..64 {
        let (r, derivative, force) = residual(end);
        if !r.is_finite() || !derivative.is_finite() || !force.is_finite() || force < 0.0 {
            return Err("nonfinite or tensile discrete felt force");
        }
        // A bracket collapsed to adjacent doubles is the root at this precision.
        if r.abs() <= floor(end, force) || hi - lo <= 4.0 * f64::EPSILON * lo.abs().max(hi.abs()) {
            return Ok(force);
        }
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
    fn felt_quadrature_tiny_positive_width_retains_tangent() {
        let law = WoolFelt::new(4.0, 1.0, 2.0, 2.0, 0.0, 2.0).unwrap();
        let a = 1e-200_f64;
        let b = a.next_up();
        let (_, slope) = average(&law, &law.initial_state(), a, b, 1.0, 1.0);
        // The force underflows but this quadratic's exact END derivative does not.
        let expected = 4.0 * (a + 2.0*b) / 3.0;
        assert!(slope.is_finite() && slope > 0.0, "lost positive tangent: {slope:e}");
        assert!((slope / expected - 1.0).abs() < 2e-14,
            "tiny-interval slope {slope:e}, analytic reference {expected:e}");
    }

    // G1: Night Orbit sample 114109, contact 108 (MIDI 67), sweeps 2 and 32.
    // Independent 100-digit Decimal primitive/secant and endpoint derivative;
    // inputs recovered as binary64, including committed sig_max and rounded strains.
    // These are arithmetic references, not a recording or a material-fit claim.
    fn captured_quadrature() -> (WoolFelt, State, [(f64, f64, f64, f64); 5]) {
        let law = WoolFelt::new(10784971.283384772, 0.2, 2.5116355500000003,
            3.2116355500000004, 0.25, 0.8).unwrap();
        let state = State { eps_max: 0.05313385293183455, sig_max: 386343.9999062737 };
        let a = 0.0003418749043622657;
        let b2 = 0.00034187507285155427;
        let b32 = 0.00034187507285085426;
        (law, state, [
            (a, b2, 6.9111815901134025, 45597.24015008384),
            (b2, a, 6.9111815901134025, 45597.21688026879),
            (a, b32, 6.911181590081485, 45597.24014989049),
            (b32, a, 6.911181590081485, 45597.21688017211),
            (a, a, 6.91117390747077, 45597.19361045862),
        ])
    }

    #[test]
    fn felt_quadrature_captured_force_matches_decimal() {
        let (law, state, cases) = captured_quadrature();
        for (start, end, expected, _) in cases {
            let (force, _) = average(&law, &state, start, end,
                0.007413793103448275, 3.3333333333333335e-5);
            assert!((force - expected).abs() < 1e-12,
                "start={start:.17e} end={end:.17e}: force={force:.17e}, reference={expected:.17e}");
        }
    }

    #[test]
    fn felt_quadrature_captured_derivative_matches_decimal() {
        let (law, state, cases) = captured_quadrature();
        for (start, end, _, expected) in cases {
            let (_, slope) = average(&law, &state, start, end,
                0.007413793103448275, 3.3333333333333335e-5);
            assert!((slope - expected).abs() < 1e-8,
                "start={start:.17e} end={end:.17e}: slope={slope:.17e}, reference={expected:.17e}");
        }
    }

    #[test]
    fn felt_quadrature_cubic_interval_boundaries() {
        // Closed-form cubic integral and endpoint derivative: no primitive subtraction.
        // Check virgin and crushed branches, both directions, equal/adjacent values,
        // both sides of the old midpoint switch and of the new series domain.
        for crushed in [false, true] {
            let law = WoolFelt::new(4.0, 1.0, if crushed { 2.0 } else { 3.0 },
                3.0, if crushed { 0.25 } else { 0.0 }, 2.0).unwrap();
            let state = if crushed { State { eps_max: 1.0, sig_max: 4.0 } }
                else { law.initial_state() };
            let origin = if crushed { 0.25 } else { 0.0 };
            let coefficient = if crushed { 4.0 / (0.75 * 0.75 * 0.75) } else { 4.0 };
            let a = 0.5_f64;
            let series_edge = |t: f64| origin + (a - origin) * (1.0 + t) / (1.0 - t);
            for b in [a, a.next_up(), a * (1.0 + 0.999e-7), a * (1.0 + 1.001e-7),
                series_edge((1.0 - 1e-8) / 24.0), series_edge((1.0 + 1e-8) / 24.0)] {
                for (start, end) in [(a, b), (b, a)] {
                    let (x, y) = (start - origin, end - origin);
                    let force = coefficient * (x*x*x + x*x*y + x*y*y + y*y*y) / 4.0;
                    let slope = coefficient * (x*x + 2.0*x*y + 3.0*y*y) / 4.0;
                    let (actual_force, actual_slope) = average(&law, &state, start, end, 1.0, 1.0);
                    assert!((actual_force - force).abs() < 1e-12,
                        "crushed={crushed} {start:.17e}->{end:.17e}: force {actual_force} vs {force}");
                    assert!((actual_slope - slope).abs() < 1e-11,
                        "crushed={crushed} {start:.17e}->{end:.17e}: slope {actual_slope} vs {slope}");
                }
            }
        }
    }

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
    #[test]
    fn published_rate_term_is_passive_on_loading_and_unloading() {
        let law = WoolFelt::new(2.0e6, 0.2, 2.5, 2.5, 0.0, 0.8).unwrap();
        let state = law.update_state(0.3, &law.initial_state());
        let (area, thickness, rate, dt) = (1e-4, 0.008, 6.47e5, 1.0/192_000.0);
        for (start, end) in [(0.0,0.001), (0.001,0.0011), (0.0011,0.0009), (0.0009,0.0)] {
            let (force, tangent) = average_with_rate(&law,&state,start,end,
                thickness,area,rate,dt);
            let elastic = average(&law,&state,start,end,thickness,area).0;
            let viscous = rate*(end.max(0.0_f64).powf(law.p)-start.max(0.0_f64).powf(law.p))/dt;
            assert!((force-(elastic+viscous).max(0.0)).abs()<1e-8);
            assert!(force>=0.0 && tangent>=0.0);
            let delta_u = area*thickness*(stored(&law,end/thickness,&state)
                -stored(&law,start/thickness,&state));
            assert!(force*(end-start)-delta_u >= -1e-10);
        }
        let force=solve_with_rate(&law,&state,0.0,0.001,2e-7,
            thickness,area,rate,dt).unwrap();
        let end=0.001-2e-7*force;
        assert!((force-average_with_rate(&law,&state,0.0,end,thickness,area,rate,dt).0).abs()<1e-5);
    }
    #[test]
    fn stiff_rate_root_resolves_force_to_the_rounding_floor() {
        // With R_H the tangent reaches R_H p e^(p-1)/dt; a fixed 1e-13 m overlap stop left
        // force errors that stalled the string-stretching iteration (1e-10 relative) above it.
        let law = WoolFelt::new(2.0e6, 0.2, 2.5, 2.5, 0.0, 0.8).unwrap();
        let state = law.initial_state();
        let (area, thickness, rate, dt) = (1e-4, 0.008, 6.47e5, 1.0/192_000.0);
        let mut worst = 0.0_f64;
        for k in 1..=400 {
            let free = 5e-6 * f64::from(k);
            for (start, c) in [(0.0, 2e-7), (0.5 * free, 2e-7), (0.0, 3e-6), (0.9 * free, 5e-8)] {
                let force = solve_with_rate(&law, &state, start, free, c, thickness, area, rate, dt).unwrap();
                let end = free - c * force;
                let (expected, tangent) = average_with_rate(&law, &state, start, end, thickness, area, rate, dt);
                let floor = 16.0 * f64::EPSILON * (end.abs() + (c * force).abs() + free.abs());
                let defect = (force - expected).abs() - (1e-9 + 1e-11 * expected.abs() + floor * tangent);
                worst = worst.max(defect);
            }
        }
        assert!(worst <= 0.0, "force exceeds the resolvable defect by {worst:e} N");
    }
}

/// Simultaneous sites on one physical hammer; uses the same discrete felt law.
#[path = "hammer_contact.rs"]
pub mod block;
