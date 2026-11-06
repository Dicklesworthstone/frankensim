//! Continuous fixed-opening evaluation for the dynamic midpoint island.
//!
//! The legacy fs-phs helper suppresses |dp| < 1e-12 Pa. That discontinuity
//! leaves some small moving-face flows without a junction root. Use the same
//! owner's Bernoulli law at unit pressure and its sqrt-pressure homogeneity;
//! this removes the numerical cutoff here without changing massless clients.

use super::viscous::{RecoveredJunction, StapleRecovery};
use crate::acoustic_realize::AcousticRealizeError;
use fs_math::det;

/// Homogeneous evaluation of the owner's continuous two-sided Bernoulli law.
/// Inputs are admitted by the dynamic-junction callers.
pub(crate) fn volume_flow(width: f64, opening: f64, dp: f64, rho: f64) -> f64 {
    let reference = fs_phs::bernoulli_volume_flow(width, opening, 1.0, rho);
    let flow = reference * det::sqrt(dp.abs());
    if dp < 0.0 { -flow } else { flow }
}

fn invalid() -> AcousticRealizeError {
    AcousticRealizeError::Reed {
        what: "continuous moving-aperture junction requires finite physical inputs and results",
    }
}

// Return both terms of dp + a*sign(dp)*sqrt(abs(dp)) = drive.
// Shared by a single moving face and condensed independent pressure faces.
pub(super) fn pressure_partition(drive: f64, a: f64) -> Result<(f64, f64), AcousticRealizeError> {
    if !drive.is_finite() || !a.is_finite() || a < 0.0 { return Err(invalid()); }
    if a == 0.0 || drive == 0.0 { return Ok((drive, 0.0)); }
    let b = 2.0 * det::sqrt(drive.abs());
    let scale = a.max(b);
    let an = a / scale;
    let bn = b / scale;
    let s = (2.0 * (drive.abs() / scale)) / (an + det::sqrt(an * an + bn * bn));
    let dp = drive.signum() * s * s;
    let jet_pressure = drive.signum() * a * s;
    if ![dp, jet_pressure].iter().all(|v| v.is_finite()) { return Err(invalid()); }
    Ok((dp, jet_pressure))
}

/// Same Bernoulli pressure-drop solution after a caller has condensed an
/// affine mechanical swept-flow compliance into the positive impedance.
/// Returns dp directly so a small drop is not reconstructed by subtracting
/// two nearly equal upstream/downstream pressures.
pub(crate) fn fixed_opening_pressure_drop(
    width: f64, opening: f64, rho: f64, impedance: f64, drive: f64,
) -> Result<f64, AcousticRealizeError> {
    if ![width, opening, rho, impedance].iter().all(|x| x.is_finite())
        || width <= 0.0 || opening < 0.0 || rho <= 0.0 || impedance <= 0.0 {
        return Err(invalid());
    }
    let a = impedance * fs_phs::bernoulli_volume_flow(width, opening, 1.0, rho);
    Ok(pressure_partition(drive, a)?.0)
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
    let bad = invalid;
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
    let (dp, jet_pressure) = pressure_partition(drive, a)?;
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

/// [`characteristic_pressure`] with the staple's pressure recovery ([`StapleRecovery`], Almeida, Vergez and Caussé 2007
/// eqs. 14-17) on the inviscid Bernoulli law: the junction the oboe cards use, which have no `slit_length` and where
/// Almeida's own model neglects the viscous loss along the reed (p. 9). With `s = sqrt((dp)_c)` the constriction drop,
/// `d = Pmouth - 2 Pminus - Z Umoving` and `a = Z Ujet(1 Pa)`: forward flow solves `f s^2 + a s = d` with
/// `f = 1 - C_P (S/S_c)^2`, and reversed flow (no recovery) is exactly [`characteristic_pressure`]'s `s^2 + a s = |d|`.
/// Returns the outgoing wave and the pieces the blade and the energy account need, like
/// `ViscousSlit::characteristic_junction_with_recovery`, to which it reduces when that slit's viscous term is negligible.
///
/// # Errors
/// As [`characteristic_pressure`], and [`STAPLE_RECOVERY_DOMAIN`](super::viscous::STAPLE_RECOVERY_DOMAIN) where forward
/// flow meets `C_P (S/S_c)^2 >= 1`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn characteristic_junction_with_recovery(
    recovery: StapleRecovery, width: f64, opening: f64, rho: f64, zc: f64,
    p_minus: f64, p_m: f64, moving_flow: f64,
) -> Result<RecoveredJunction, AcousticRealizeError> {
    let bad = invalid;
    if ![width, opening, rho, zc, p_minus, p_m, moving_flow]
        .iter().all(|v| v.is_finite())
        || width <= 0.0 || opening < 0.0 || rho <= 0.0 || zc <= 0.0
    {
        return Err(bad());
    }
    let no_jet = p_minus + zc * moving_flow;
    let no_drop = p_m - p_minus;
    let drive = no_drop - no_jet;
    let unit_flow = fs_phs::bernoulli_volume_flow(width, opening, 1.0, rho);
    let a = zc * unit_flow;
    if ![no_jet, no_drop, drive, a].iter().all(|v| v.is_finite()) {
        return Err(bad());
    }
    if a == 0.0 || drive == 0.0 {
        // No flow, so no recovery: the whole mouth-to-bore difference `drive` falls across the shut reed, and it is that
        // difference the blade feels.
        return Ok(RecoveredJunction {
            outgoing_pa: no_jet, flow_m3_s: 0.0, constriction_drop_pa: drive, recovered_pa: 0.0,
        });
    }
    let factor = if drive > 0.0 { recovery.forward_factor(width * opening)? } else { 1.0 };
    // sqrt(|(dp)_c|) = s solves f s^2 + a s = |d|: s = 2|d| / (a + sqrt(a^2 + 4 f |d|)), scaled as `pressure_partition` does.
    let b = 2.0 * det::sqrt(factor * drive.abs());
    let scale = a.max(b);
    let (an, bn) = (a / scale, b / scale);
    let s = (2.0 * (drive.abs() / scale)) / (an + det::sqrt(an * an + bn * bn));
    let sign = drive.signum();
    let constriction_drop = sign * s * s;
    let flow = sign * unit_flow * s;
    let jet_pressure = sign * a * s;
    let recovered = recovery.recovered_pressure_pa(rho, flow);
    let dp = constriction_drop - recovered;
    let outgoing = if no_drop.abs().max(dp.abs()) <= no_jet.abs().max(jet_pressure.abs()) {
        no_drop - dp
    } else {
        no_jet + jet_pressure
    };
    if ![flow, constriction_drop, recovered, dp, jet_pressure, outgoing].iter().all(|v| v.is_finite()) {
        return Err(bad());
    }
    Ok(RecoveredJunction {
        outgoing_pa: outgoing, flow_m3_s: flow, constriction_drop_pa: constriction_drop, recovered_pa: recovered,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::viscous::{STAPLE_RECOVERY_DOMAIN, ViscousSlit};

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

    // ---- staple pressure recovery on the inviscid law (Almeida, Vergez and Caussé 2007, eqs. 14-17) -------------------
    // Fig. 14 and Fig. 15: S_c = 4e-6 m^2, C_P = 0.8.
    const S_C: f64 = 4e-6;

    #[test]
    fn staple_junction_on_the_inviscid_law_is_the_landed_law_when_the_slit_viscosity_is_negligible() {
        let recovery = StapleRecovery::new(0.8, S_C).unwrap();
        // 12 mu L / h^2 is about 1e-12 Pa s/m at the smallest opening below: 1e-14 of the Bernoulli head.
        let slit = ViscousSlit::new(1e-12, 1e-9).unwrap();
        let (width, rho) = (0.0056, 1.2);
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()));
        let mut forward = 0;
        for opening in [0.0, 5e-5, 1.5e-4, 3.57e-4, 6e-4] {
            for zc in [3e5, 6e7, 5e9] {
                for incoming in [-300.0, 0.0, 1000.0] {
                    for mouth in [-4000.0, -50.0, 50.0, 800.0, 7400.0, 30000.0] {
                        for moving in [-8e-5, 0.0, 6e-5] {
                            let a = characteristic_junction_with_recovery(
                                recovery, width, opening, rho, zc, incoming, mouth, moving).unwrap();
                            let b = slit.characteristic_junction_with_recovery(
                                recovery, width, opening, rho, zc, incoming, mouth, moving).unwrap();
                            assert!(close(a.outgoing_pa, b.outgoing_pa)
                                && close(a.constriction_drop_pa, b.constriction_drop_pa)
                                && close(a.recovered_pa, b.recovered_pa)
                                && (a.flow_m3_s - b.flow_m3_s).abs()
                                    <= 1e-9 * (1e-9 + a.flow_m3_s.abs().max(b.flow_m3_s.abs())),
                                "opening={opening} zc={zc} in={incoming} m={mouth} mv={moving}: {a:?} vs {b:?}");
                            if a.flow_m3_s > 0.0 { forward += 1; }
                        }
                    }
                }
            }
        }
        assert!(forward > 100, "{forward}");
    }

    #[test]
    fn staple_junction_with_c_p_zero_or_reversed_flow_is_characteristic_pressure_bit_for_bit() {
        let none = StapleRecovery::new(0.0, S_C).unwrap();
        let some = StapleRecovery::new(0.8, S_C).unwrap();
        for opening in [0.0, 1e-4, 3.57e-4, 8e-4] {
            for mouth in [-3000.0, -50.0, 800.0, 7400.0] {
                for (incoming, zc, moving) in [(75.0, 5.9e7, 1e-7), (-30.0, 3e5, -2e-8), (0.0, 6e7, 0.0)] {
                    let plain = characteristic_pressure(0.0056, opening, 1.2, zc, incoming, mouth, moving).unwrap();
                    let zero = characteristic_junction_with_recovery(
                        none, 0.0056, opening, 1.2, zc, incoming, mouth, moving).unwrap();
                    assert_eq!(zero.outgoing_pa.to_bits(), plain.to_bits());
                    // Reversed flow (a negative drive) has no recovery whatever C_P, and no area limit either: the
                    // opening 8e-4 is above S_c / sqrt(C_P) = 4.47 mm^2 and is only admitted there.
                    if mouth - incoming - (incoming + zc * moving) < 0.0 {
                        let reversed = characteristic_junction_with_recovery(
                            some, 0.0056, opening, 1.2, zc, incoming, mouth, moving).unwrap();
                        assert_eq!(reversed.outgoing_pa.to_bits(), plain.to_bits());
                        assert_eq!(reversed.recovered_pa, 0.0);
                    }
                }
            }
        }
    }

    #[test]
    fn staple_junction_returns_the_manufactured_forward_and_reversed_flows() {
        let recovery = StapleRecovery::new(0.8, S_C).unwrap();
        let (width, rho, zc, incoming, moving): (f64, f64, f64, f64, f64) = (0.0056, 1.2, 5.9e7, 75.0, 1e-7);
        for opening in [5e-5, 1.5e-4, 3.57e-4] {
            let area = width * opening;
            // s is the signed root of the constriction drop: (dp)_c = sign * s^2, q = sign * A sqrt(2/rho) s.
            for s in [-60.0_f64, -3.0, 3.0, 30.0, 70.0] {
                let q = area * (2.0 / rho).sqrt() * s;
                let constriction = s.signum() * s * s;
                let recovered = if q > 0.0 { 0.8 * 0.5 * rho * (q / S_C) * (q / S_C) } else { 0.0 };
                let expected = incoming + zc * (q + moving);
                let mouth = (constriction - recovered) + expected + incoming;
                let j = characteristic_junction_with_recovery(
                    recovery, width, opening, rho, zc, incoming, mouth, moving).unwrap();
                assert!((j.outgoing_pa - expected).abs() <= 1e-10 * (1.0 + expected.abs()), "{j:?} vs {expected}");
                assert!((j.flow_m3_s - q).abs() <= 1e-9 * q.abs(), "{j:?} vs {q}");
                assert!((j.constriction_drop_pa - constriction).abs() <= 1e-9 * constriction.abs(), "{j:?}");
                assert!((j.recovered_pa - recovered).abs() <= 1e-9 * recovered.abs().max(1e-12), "{j:?}");
            }
        }
    }

    #[test]
    fn staple_junction_refuses_by_name_and_a_shut_channel_hands_the_blade_the_whole_difference() {
        let recovery = StapleRecovery::new(0.8, S_C).unwrap();
        let width = 0.0056;
        let limit_opening = S_C / 0.8_f64.sqrt() / width; // C_P (S/S_c)^2 = 1
        let junction = |opening: f64, mouth: f64| characteristic_junction_with_recovery(
            recovery, width, opening, 1.2, 1e7, 0.0, mouth, 0.0);
        assert!(junction(limit_opening * (1.0 - 1e-9), 500.0).is_ok());
        for opening in [limit_opening * (1.0 + 1e-9), 1.5 * limit_opening] {
            match junction(opening, 500.0) {
                Err(AcousticRealizeError::Reed { what }) => assert_eq!(what, STAPLE_RECOVERY_DOMAIN),
                other => panic!("expected the staple domain refusal, got {other:?}"),
            }
        }
        // Reversed flow: admitted at the same opening, and exactly the original junction.
        let reversed = junction(1.5 * limit_opening, -500.0).unwrap();
        let original = characteristic_pressure(width, 1.5 * limit_opening, 1.2, 1e7, 0.0, -500.0, 0.0).unwrap();
        assert_eq!(reversed.outgoing_pa.to_bits(), original.to_bits());
        // A shut channel: no flow, no recovery, and the blade feels the whole mouth-to-bore difference.
        let shut = characteristic_junction_with_recovery(recovery, width, 0.0, 1.2, 6e7, 75.0, 7400.0, 1e-7).unwrap();
        assert_eq!(shut.outgoing_pa, 75.0 + 6e7 * 1e-7);
        assert_eq!((shut.flow_m3_s, shut.recovered_pa), (0.0, 0.0));
        assert_eq!(shut.constriction_drop_pa, (7400.0 - 75.0) - (75.0 + 6e7 * 1e-7));
        for bad in [f64::NAN, f64::INFINITY] {
            assert!(junction(3e-4, bad).is_err());
        }
    }
}
