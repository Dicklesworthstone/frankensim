//! Bounded inventory algebra for the original well-mixed membrane stage.
//!
//! This module supplies no saturation or membrane constitutive law. Its inputs
//! are capacities from the shared material phase owner and the native donor
//! coefficients at one candidate solvent transfer.

/// Solve d = a min(n_left-d, capacity_left)
///         - b min(n_right+d, capacity_right).
///
/// A missing capacity denotes an ordinary aqueous cell with no solid phase.
/// The residual has slope 1 + a I_left_unsaturated + b I_right_unsaturated:
/// strictly positive, with at most two phase boundaries. Enumerating its four
/// possible active sets is bounded work; at most three have nonempty intervals.
/// Every endpoint uses one shared transfer, with positive implicit recovery of
/// a nearly depleted donor. No crystal crosses the membrane, no amount is
/// clipped, and neither an iterative inner solve nor heap storage is needed.
pub(super) fn finite_stock_exchange(
    initial: [f64; 2],
    capacities: [Option<f64>; 2],
    coefficients: [f64; 2],
) -> Result<[f64; 2], &'static str> {
    if initial
        .iter()
        .chain(&coefficients)
        .any(|x| !x.is_finite() || *x < 0.0)
        || capacities
            .iter()
            .flatten()
            .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err("invalid finite-stock membrane inventory or coefficient");
    }
    let [nl, nr] = initial;
    let [a, b] = coefficients;
    let scale = 1.0_f64.max(a).max(b);
    let a = a / scale;
    let b = b / scale;
    let z = 1.0 / scale;
    for saturated in [[false, false], [true, false], [false, true], [true, true]] {
        if (saturated[0] && capacities[0].is_none()) || (saturated[1] && capacities[1].is_none()) {
            continue;
        }
        let dl = if saturated[0] {
            capacities[0].unwrap()
        } else {
            nl
        };
        let dr = if saturated[1] {
            capacities[1].unwrap()
        } else {
            nr
        };
        let denominator =
            z + if saturated[0] { 0.0 } else { a } + if saturated[1] { 0.0 } else { b };
        let transfer = (a * dl - b * dr) / denominator;
        if !transfer.is_finite() {
            continue;
        }
        let endpoint = if transfer > 0.5 * nl && !saturated[0] {
            // Rearrange the left inventory equation into positive terms.
            let from_left = z + if saturated[1] { 0.0 } else { b };
            let left = (from_left / denominator) * nl + (b / denominator) * dr;
            [left, nr + (nl - left)]
        } else if -transfer > 0.5 * nr && !saturated[1] {
            let from_right = z + if saturated[0] { 0.0 } else { a };
            let right = (from_right / denominator) * nr + (a / denominator) * dl;
            [nl + (nr - right), right]
        } else {
            [nl - transfer, nr + transfer]
        };
        if endpoint.iter().any(|x| !x.is_finite() || *x < 0.0) {
            continue;
        }
        if (0..2).all(|side| matches_phase(endpoint[side], capacities[side], saturated[side])) {
            return Ok(endpoint);
        }
    }
    Err("finite-stock implicit membrane endpoint is not representable")
}

fn matches_phase(total: f64, capacity: Option<f64>, saturated: bool) -> bool {
    match capacity {
        None => !saturated,
        Some(capacity) => {
            // Adjacent affine pieces share the same physical endpoint. Admit
            // roundoff at their common boundary without clipping either stock.
            // The tolerance is local to this phase capacity, never a remote
            // reservoir or a combined total that could erase a dissolved trace.
            let roundoff = 8.0 * f64::EPSILON * total.max(capacity);
            if saturated {
                total >= capacity || capacity - total <= roundoff
            } else {
                total <= capacity || total - capacity <= roundoff
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3e-14 * expected.abs().max(1.0),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn finite_reservoir_exhausts_and_receiver_can_precipitate() {
        // Saturated left, dilute right: d = 2*3 - (1+d).
        let first = finite_stock_exchange([10.0, 1.0], [Some(3.0), Some(5.0)], [2.0, 1.0]).unwrap();
        near(first[0], 7.5);
        near(first[1], 3.5);

        // Donor falls below saturation within this very step:
        // d = 4*(4-d) - 0.5*(1+d), giving d=31/11.
        let exhausted =
            finite_stock_exchange([4.0, 1.0], [Some(3.0), Some(5.0)], [4.0, 0.5]).unwrap();
        near(exhausted[0], 13.0 / 11.0);
        near(exhausted[1], 42.0 / 11.0);

        // Both endpoints saturated: d=4*3 - 0.5*2 = 11.
        let crystallizing =
            finite_stock_exchange([20.0, 0.0], [Some(3.0), Some(2.0)], [4.0, 0.5]).unwrap();
        near(crystallizing[0], 9.0);
        near(crystallizing[1], 11.0);
    }

    #[test]
    fn orientation_reversal_and_aqueous_receiver_preserve_stock() {
        let forward = finite_stock_exchange([8.0, 1.0], [Some(2.0), None], [3.0, 0.25]).unwrap();
        let reverse = finite_stock_exchange([1.0, 8.0], [None, Some(2.0)], [0.25, 3.0]).unwrap();
        assert_eq!(forward, [reverse[1], reverse[0]]);
        near(forward[0] + forward[1], 9.0);
        near(forward[0], 3.4);
        near(forward[1], 5.6);
    }

    #[test]
    fn strong_exchange_retains_positive_trace_and_tiny_exchange_is_resolved() {
        let trace = finite_stock_exchange([1.0, 0.0], [Some(0.1), None], [1e250, 0.0]).unwrap();
        assert!(trace[0] > 0.0);
        assert!((trace[0] / 1e-250 - 1.0).abs() < 1e-14);
        assert_eq!(trace[1], 1.0);

        let tiny = finite_stock_exchange([1e20, 0.0], [Some(1.0), None], [1e-12, 0.0]).unwrap();
        assert_eq!(tiny[0], 1e20);
        assert_eq!(tiny[1], 1e-12);
    }

    #[test]
    fn exact_phase_boundaries_and_closed_membrane_are_admitted() {
        for coefficients in [[0.0, 0.0], [1.0, 1.0], [1e200, 1e200]] {
            let endpoint =
                finite_stock_exchange([3.0, 3.0], [Some(3.0), Some(3.0)], coefficients).unwrap();
            assert_eq!(endpoint, [3.0, 3.0]);
        }
        let endpoint =
            finite_stock_exchange([5.0, 0.0], [Some(3.0), Some(2.0)], [1.0, 0.5]).unwrap();
        assert_eq!(endpoint, [3.0, 2.0]);
    }

    #[test]
    fn original_inventory_equations_hold_across_all_phase_combinations() {
        // Independent residual/balance oracle across both directions, phase
        // transitions, asymmetric coefficients, and absent solid phases.
        for nl in [0.0, 0.5, 3.0, 20.0] {
            for nr in [0.0, 0.75, 2.0, 30.0] {
                for left in [None, Some(0.25), Some(3.0)] {
                    for right in [None, Some(0.5), Some(2.0)] {
                        for a in [0.0, 0.1, 1.0, 100.0] {
                            for b in [0.0, 0.2, 1.0, 80.0] {
                                let end =
                                    finite_stock_exchange([nl, nr], [left, right], [a, b]).unwrap();
                                let dl = left.map_or(end[0], |s: f64| end[0].min(s));
                                let dr = right.map_or(end[1], |s: f64| end[1].min(s));
                                let transfer = a * dl - b * dr;
                                let tolerance = 2e-12 * (nl + nr + transfer.abs()).max(1.0);
                                assert!((end[0] + end[1] - nl - nr).abs() <= tolerance);
                                assert!((end[0] - nl + transfer).abs() <= tolerance);
                                assert!((end[1] - nr - transfer).abs() <= tolerance);
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn invalid_inputs_refuse_without_inventing_a_state() {
        assert!(finite_stock_exchange([1.0, -1.0], [None, None], [1.0, 1.0]).is_err());
        assert!(finite_stock_exchange([1.0, 0.0], [Some(0.0), None], [1.0, 1.0]).is_err());
        assert!(finite_stock_exchange([1.0, 0.0], [None, None], [f64::INFINITY, 1.0]).is_err());
    }
}
