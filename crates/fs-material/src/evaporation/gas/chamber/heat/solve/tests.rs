//! The joint liquid heat search at an enthalpy datum where one represented
//! enthalpy moves the residual by more than the search's own gate.
use super::*;

const MASS: f64 = 8e-6;
const K: f64 = 1.478e-5;
const DATUM: f64 = 100_449_999.0;

fn temperature(h: f64) -> Result<f64, EvaporationError> {
    Ok(300.0 + (h - DATUM) / 4000.0)
}

fn at(h: f64) -> LiquidSample {
    LiquidSample {
        h,
        t: 300.0 + (h - DATUM) / 4000.0,
        q: MASS * (h - DATUM),
    }
}

fn response(tolerance: f64, span: f64) -> Response {
    Response {
        index: 0,
        h0: DATUM,
        t0: 300.0,
        mass: MASS,
        k: K,
        lo: at(DATUM),
        hi: at(DATUM + span),
        temperature_tolerance: tolerance,
        maximum_iterations: 128,
    }
}

// A gas temperature whose root of the residual lies exactly between two
// adjacent represented enthalpies near 1.0045e8 J/kg, one ulp (1.49e-8 J/kg)
// apart: each of them is then 0.5*ulp*(1/cp + m/k) = 4.03e-9 K from the root.
fn gas_between(h: f64, next: f64) -> f64 {
    let level = |h: f64| (300.0 + (h - DATUM) / 4000.0) + MASS * (h - DATUM) / K;
    0.5 * (level(h) + level(next))
}

fn search(
    tolerance: f64,
    span: f64,
    gas_t: f64,
) -> Result<(LiquidSample, f64), GasEvaporationError> {
    let mut work = 0usize;
    liquid_response(
        &response(tolerance, span),
        gas_t,
        &temperature,
        1e-9,
        None,
        &mut work,
        &mut || false,
    )
}

// The better endpoint of the adjacent pair (h, next), and its residual.
fn better_endpoint(tolerance: f64, gas_t: f64, h: f64, next: f64) -> f64 {
    let probe = response(tolerance, 6.0);
    let residual = |h: f64| probe.residual(at(h), gas_t).unwrap();
    assert!(residual(h) < 0.0 && residual(next) > 0.0);
    residual(h).abs().min(residual(next).abs())
}

#[test]
fn a_one_ulp_bracket_takes_its_better_endpoint_within_the_declared_tolerance_only() {
    let h = DATUM + 2.1;
    let next = h.next_up();
    let gas_t = gas_between(h, next);
    // Both endpoints sit outside 0.25 of the 1e-8 K tolerance and inside the
    // tolerance, so only the one-ulp collapse rule can accept either.
    let better = better_endpoint(1e-8, gas_t, h, next);
    assert!(better > 0.25e-8 && better < 1e-8, "{better:e}");

    let (accepted, error) = search(1e-8, 6.0, gas_t).unwrap();
    assert!(accepted.h == h || accepted.h == next);
    let probe = response(1e-8, 6.0);
    assert_eq!(probe.residual(accepted, gas_t).unwrap().abs(), better);
    assert!(error <= 1e-9);

    // The planted negative: the same one-ulp bracket against a declared
    // tolerance below its better endpoint still refuses.
    assert!(better > 2e-9);
    assert!(matches!(
        search(2e-9, 6.0, gas_t),
        Err(GasEvaporationError::HeatNotConverged { .. })
    ));
}

#[test]
fn a_secant_that_rounds_onto_a_wide_bracket_end_bisects_on_to_adjacent_enthalpies() {
    // Starting from datum+4 the safeguarded secant rounds onto an end of the
    // bracket while it is still three ulps wide (residuals -4.0e-9 and
    // +2.0e-8 K at that point). That bracket is not a collapse: the search
    // bisects on, closes on the adjacent pair around the root, and takes the
    // better endpoint within the declared tolerance exactly as above.
    let h = (DATUM + 1.0).next_up();
    let next = h.next_up();
    let gas_t = gas_between(h, next);
    let better = better_endpoint(1e-8, gas_t, h, next);
    let (accepted, _) = search(1e-8, 4.0, gas_t).unwrap();
    assert!(accepted.h == h || accepted.h == next);
    let probe = response(1e-8, 4.0);
    assert_eq!(probe.residual(accepted, gas_t).unwrap().abs(), better);
    assert!(matches!(
        search(2e-9, 4.0, gas_t),
        Err(GasEvaporationError::HeatNotConverged { .. })
    ));
}
