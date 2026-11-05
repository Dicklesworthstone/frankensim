//! Positive completion of the unresolved fast diffusion spectrum.
//!
//! For a positive spectral measure dR(lambda), its force impedance is
//! integral s/(s+lambda) dR. The first two inverse-rate moments are the
//! added inertia M = integral dR/lambda and the leading dissipative
//! coefficient N = integral dR/lambda^2. One physical mass/dashpot with
//! lambda_eff=M/N and R_eff=M^2/N retains BOTH moments and is passive.
//! We apply this to the final log cell and its entire infinite fast tail,
//! so completing the continuum adds no coordinates to the coupled solver.

use super::*;
use fs_aeroac::bessel::{hankel0_abs_squared, hankel1_abs_squared};

// Eight-point Gauss-Legendre integration, with positive weights. The mapped
// integrands are smooth at u=0: u/B_n(x/u) -> pi*x/2 and u^3/B_n(x/u) -> 0.
// Nodes never reach the endpoint, so no infinite Bessel argument is evaluated.
const NODES: [f64; 4] = [
    0.183_434_642_495_649_8, 0.525_532_409_916_329,
    0.796_666_477_413_626_7, 0.960_289_856_497_536_3,
];
const WEIGHTS: [f64; 4] = [
    0.362_683_783_378_362, 0.313_706_645_877_887_3,
    0.222_381_034_453_374_5, 0.101_228_536_290_376_3,
];

fn integrate(lo: f64, hi: f64, kernel: &impl Fn(f64) -> [f64; 2]) -> [f64; 2] {
    let mid = 0.5 * (lo + hi);
    let half = 0.5 * (hi - lo);
    let mut sum = [0.0; 2];
    for (node, weight) in NODES.iter().zip(WEIGHTS) {
        for sign in [-1.0, 1.0] {
            let point = kernel(mid + sign * half * node);
            for i in 0..2 { sum[i] += half * weight * point[i]; }
        }
    }
    sum
}

fn refine(lo: f64, hi: f64, coarse: [f64; 2], depth: usize,
    kernel: &impl Fn(f64) -> [f64; 2], gate: &CancelGate)
    -> Result<[f64; 2], AcousticRealizeError>
{
    if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
    let mid = 0.5 * (lo + hi);
    let left = integrate(lo, mid, kernel);
    let right = integrate(mid, hi, kernel);
    let fine = [left[0] + right[0], left[1] + right[1]];
    if fine.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(invalid("cylinder fast-tail moments are not representable"));
    }
    // Resolve the cold quadrature far more tightly than the user-controlled
    // log grid. This is an integration tolerance, not a constitutive constant
    // or a certificate for the resulting finite-frequency approximation.
    if (0..2).all(|i| (fine[i] - coarse[i]).abs() <= 1e-10 * fine[i]) {
        return Ok(fine);
    }
    if depth == 20 {
        return Err(invalid("cylinder fast-tail moment quadrature did not converge"));
    }
    let left = refine(lo, mid, left, depth + 1, kernel, gate)?;
    let right = refine(mid, hi, right, depth + 1, kernel, gate)?;
    Ok([left[0] + right[0], left[1] + right[1]])
}

/// Dimensionless inertia and resistance weights for [x, infinity).
/// Keep the dimensional inertia calculation independent of fluid viscosity.
pub(super) fn weights(x: f64, rotary: bool, gate: &CancelGate)
    -> Result<(f64, f64), AcousticRealizeError>
{
    if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
    let norm = if rotary { hankel1_abs_squared } else { hankel0_abs_squared };
    // Above 18, inverse wavenumber makes both integrands regular at infinity.
    // Below 18, log-wavenumber panels resolve arbitrarily many slow decades;
    // plain inverse-wavenumber subdivision would need depth proportional to
    // -log(x), and could reject an otherwise valid one-interval realization.
    let split = x.max(18.0);
    let reciprocal = |u: f64| {
        let first = u / norm(split / u);
        [first, first * u * u]
    };
    let mut moments = refine(0.0, 1.0, integrate(0.0, 1.0, &reciprocal),
        0, &reciprocal, gate)?;
    if x < split {
        let ratio_squared = (x / split).powi(2);
        moments[0] *= ratio_squared;
        moments[1] *= ratio_squared * ratio_squared;
        let log_min = det::ln(x);
        let span = det::ln(split) - log_min;
        let panels = span.ceil() as usize;
        let logarithmic = |log_k: f64| {
            let relative_squared = det::exp(2.0 * (log_min - log_k));
            let first = relative_squared / norm(det::exp(log_k));
            [first, first * relative_squared]
        };
        for panel in 0..panels {
            let lo = log_min + span * panel as f64 / panels as f64;
            let hi = log_min + span * (panel + 1) as f64 / panels as f64;
            let part = refine(lo, hi, integrate(lo, hi, &logarithmic),
                0, &logarithmic, gate)?;
            for i in 0..2 { moments[i] += part[i]; }
        }
    }
    let [first, second] = moments;
    // lambda = (mu/rho/a^2) k^2; u=x/k. The normalized positive
    // moments avoid integrating on an infinite interval or subtracting tails.
    let mass = first / x / x;
    let resistance = first * (first / second);
    if !mass.is_finite() || mass <= 0.0 || !resistance.is_finite() || resistance <= 0.0 {
        return Err(invalid("cylinder fast-tail mass/dashpot is not representable"));
    }
    Ok((mass, resistance))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn low_final_cutoffs_and_one_cell_spectra_do_not_exhaust_tail_subdivision() {
        let gate = CancelGate::new_clock_free();
        for x in [1e-8, 1e-5, 0.001, 0.1, 18.0, 3000.0] {
            for rotary in [false, true] {
                let (mass, resistance) = weights(x, rotary, &gate).unwrap();
                assert!(mass.is_finite() && mass > 0.0);
                assert!(resistance.is_finite() && resistance > 0.0);
                // Every represented rate was at least x^2. Matching inverse
                // moments must leave the effective rate within that support.
                assert!(resistance / mass >= x*x);
            }
        }
    }

    #[test]
    fn fast_tail_preparation_observes_the_original_cancellation_gate() {
        let gate = CancelGate::new_clock_free();
        gate.request();
        for rotary in [false, true] {
            assert!(matches!(weights(0.001, rotary, &gate), Err(AcousticRealizeError::Cancelled)));
        }
    }
}
