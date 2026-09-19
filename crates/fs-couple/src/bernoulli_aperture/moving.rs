//! Continuous fixed-opening evaluation for the dynamic midpoint island.
//!
//! The legacy fs-phs helper suppresses |dp| < 1e-12 Pa. That discontinuity
//! leaves some small moving-face flows without a junction root. Use the same
//! owner's Bernoulli law at unit pressure and its sqrt-pressure homogeneity;
//! this removes the numerical cutoff here without changing massless clients.

use crate::acoustic_realize::AcousticRealizeError;
use fs_math::det;

/// Homogeneous evaluation of the owner's continuous two-sided Bernoulli law.
/// Inputs are admitted by the dynamic-junction callers.
pub(crate) fn volume_flow(width: f64, opening: f64, dp: f64, rho: f64) -> f64 {
    let reference = fs_phs::bernoulli_volume_flow(width, opening, 1.0, rho);
    let flow = reference * det::sqrt(dp.abs());
    if dp < 0.0 { -flow } else { flow }
}

/// Solve the fixed-opening characteristic junction algebraically.
/// With d = Pmouth - 2 Pminus - Z Umoving and a = Z Ujet(1 Pa),
/// sqrt(|dp|) solves s^2 + a s = |d|. Rationalization avoids cancellation
/// at rest; a scaled hypotenuse avoids squaring a large coefficient.
#[allow(clippy::too_many_arguments)]
pub(crate) fn characteristic_pressure(
    width: f64, opening: f64, rho: f64, zc: f64,
    p_minus: f64, p_m: f64, moving_flow: f64,
) -> Result<f64, AcousticRealizeError> {
    let bad = || AcousticRealizeError::Reed {
        what: "continuous moving-aperture junction requires finite physical inputs and results",
    };
    if ![width, opening, rho, zc, p_minus, p_m, moving_flow]
        .iter().all(|v| v.is_finite())
        || width <= 0.0 || opening < 0.0 || rho <= 0.0 || zc <= 0.0
    {
        return Err(bad());
    }
    let no_jet = p_minus + zc * moving_flow;
    let no_drop = p_m - p_minus;
    let drive = no_drop - no_jet;
    let a = zc * fs_phs::bernoulli_volume_flow(width, opening, 1.0, rho);
    if ![no_jet, no_drop, drive, a].iter().all(|v| v.is_finite()) {
        return Err(bad());
    }
    if a == 0.0 || drive == 0.0 {
        return Ok(no_jet);
    }
    let b = 2.0 * det::sqrt(drive.abs());
    let scale = a.max(b);
    let an = a / scale;
    let bn = b / scale;
    let s = (2.0 * (drive.abs() / scale)) / (an + det::sqrt(an * an + bn * bn));
    let dp = drive.signum() * s * s;
    let jet_pressure = drive.signum() * a * s;
    // These are algebraically identical. Choose the subtraction/addition with
    // the smaller operand scale to avoid destroying a small outgoing pressure.
    let outgoing = if no_drop.abs().max(dp.abs()) <= no_jet.abs().max(jet_pressure.abs()) {
        no_drop - dp
    } else {
        no_jet + jet_pressure
    };
    if ![dp, jet_pressure, outgoing].iter().all(|v| v.is_finite()) {
        return Err(bad());
    }
    Ok(outgoing)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_signed_flows_have_roots_below_the_legacy_dead_zone() {
        let (width, opening, rho, zc) = (0.013, 4e-4, 1.2, 1e6);
        for magnitude in [1e-6_f64, 1e-12, 1e-20, 1e-40, 1e-200] {
            for sign in [-1.0, 1.0] {
                let dp = sign * magnitude;
                let jet = width * opening * sign * (2.0 * magnitude / rho).sqrt();
                let moving = -dp / zc - jet;
                let outgoing = characteristic_pressure(width, opening, rho, zc, 0.0, 0.0, moving).unwrap();
                assert!((outgoing + dp).abs() <= 2e-12 * magnitude);
                let actual_jet = volume_flow(width, opening, -outgoing, rho);
                let wave = outgoing / zc;
                let scale = actual_jet.abs() + moving.abs() + wave.abs();
                assert!((actual_jet + moving - wave).abs() <= 2e-12 * scale);
                if magnitude < 1e-12 {
                    assert_eq!(fs_phs::bernoulli_volume_flow(width, opening, dp, rho), 0.0);
                    assert!(actual_jet.abs() > 0.0, "continuous response must not acquire a dead zone");
                }
            }
        }
    }

    #[test]
    fn algebraic_junction_matches_manufactured_forward_and_reverse_flows() {
        for opening in [0.0, 1e-5, 4e-4] {
            for dp in [-3000.0_f64, -250.0, 250.0, 3000.0] {
                for incoming in [-200.0, 0.0, 200.0] {
                    for moving in [-1e-6, 0.0, 1e-6] {
                        let jet = 0.013 * opening * dp.signum() * (2.0 * dp.abs() / 1.2).sqrt();
                        let expected = incoming + 1e6 * (jet + moving);
                        let mouth = dp + expected + incoming;
                        let actual = characteristic_pressure(0.013, opening, 1.2, 1e6, incoming, mouth, moving).unwrap();
                        assert!((actual - expected).abs() <= 2e-12 * (1.0 + expected.abs()));
                    }
                }
            }
        }
    }

    #[test]
    fn equal_pressure_and_closed_opening_do_not_invent_motion_or_flow() {
        assert_eq!(characteristic_pressure(0.013, 4e-4, 1.2, 1e6, 0.0, 0.0, 0.0).unwrap(), 0.0);
        assert_eq!(characteristic_pressure(0.013, 0.0, 1.2, 1e6, 75.0, 1000.0, 2e-7).unwrap(), 75.2);
        for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(characteristic_pressure(0.013, 4e-4, 1.2, 1e6, 0.0, 0.0, bad).is_err());
        }
    }
}
