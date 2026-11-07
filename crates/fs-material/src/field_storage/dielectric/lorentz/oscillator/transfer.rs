//! Original Lorentz reservoir channels for coherent, settled harmonic fields.

use super::{LorentzPole, Result, finite};
use fs_math::{c64::C64, det, qd::Qd};

/// Transfer from electric field to the original normalized polarization
/// reservoirs, retaining phase under `exp(+i omega t)`.
///
/// For a coherent real field, apply each channel to every electric-field
/// phasor before summing. The instantaneous bound-charge energy density is
/// `(restoring(t)^2 + kinetic(t)^2)/2`, and the frictional heating power
/// density is `heat(t)^2`, summed over vector components. Products must retain
/// the cross-frequency terms; summing individual tone powers loses beats.
/// These are field transfers, not squared magnitudes or cycle averages.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LorentzHarmonicTransfer {
    /// Transfer to `P/sqrt(delta_epsilon)`, the restoring-energy coordinate.
    pub restoring: C64,
    /// Transfer to `Pdot/(omega0*sqrt(delta_epsilon))`, the kinetic coordinate.
    pub kinetic: C64,
    /// Transfer to `sqrt(gamma)*kinetic`, whose real square is heat rate.
    pub heat: C64,
}

impl LorentzPole {
    /// Prepare physical storage and dissipation channels once per frequency.
    ///
    /// This is the settled spectral response of the original oscillator. It
    /// retains current inertia and phase, including above resonance where
    /// susceptibility's real part can be negative but storage stays positive.
    /// No switch-on transient, boundary loss, or independent-tone absorption
    /// approximation is inferred. Use [`Self::prepare_harmonic`] when an
    /// initial polarization/current state and its transient are required.
    ///
    /// Frequency is finite and nonnegative [rad/s]. The shared susceptibility
    /// owner supplies the physical response and its range/pole refusals. An
    /// exactly resonant lossless pole has no finite settled response. Zero
    /// damping gives exactly zero heat, and zero frequency gives static
    /// restoring storage with no kinetic storage or heat. No allocation or
    /// quadrature is needed.
    pub fn harmonic_transfer(
        self,
        angular_frequency_rad_s: f64,
    ) -> Result<LorentzHarmonicTransfer> {
        let omega = angular_frequency_rad_s;
        let susceptibility = self.susceptibility_f_m(omega)?;
        if omega == 0.0 {
            return Ok(LorentzHarmonicTransfer {
                restoring: C64::from_re(self.root_increment),
                kinetic: C64::ZERO,
                heat: C64::ZERO,
            });
        }
        let normalized = |response: C64, numerator: &[f64], denominator: &[f64]| -> Result<C64> {
            Ok(C64::new(
                scaled_component(response.re, numerator, denominator)?,
                scaled_component(response.im, numerator, denominator)?,
            ))
        };
        let restoring = normalized(susceptibility, &[], &[self.root_increment])?;
        let derivative = C64::new(-susceptibility.im, susceptibility.re);
        let normalization = [self.omega, self.root_increment];
        let kinetic = normalized(derivative, &[omega], &normalization)?;
        let heat = if self.damping == 0.0 {
            C64::ZERO
        } else {
            // Form heat from the original response, not a rounded kinetic
            // component that could underflow before sqrt(gamma) restores it.
            normalized(
                derivative,
                &[omega, det::sqrt(self.damping)],
                &normalization,
            )?
        };
        Ok(LorentzHarmonicTransfer {
            restoring,
            kinetic,
            heat,
        })
    }
}

// Separate binary exponents before multiplying dimensional coefficients.
// All source factors are positive and finite. Mantissas stay near unity;
// the shared math owner applies the final power of two. Only the final
// unrepresentable component can underflow; intermediate factors cannot erase
// a representable kinetic/heat component or spuriously overflow it.
fn scaled_component(value: f64, numerator: &[f64], denominator: &[f64]) -> Result<f64> {
    if value == 0.0 {
        return Ok(value);
    }
    let (mut mantissa, mut exponent) = binary_parts(value);
    for &factor in numerator {
        let (m, e) = binary_parts(factor);
        mantissa *= m;
        exponent += e;
    }
    for &factor in denominator {
        let (m, e) = binary_parts(factor);
        mantissa /= m;
        exponent -= e;
    }
    finite(Qd::from_f64(mantissa).scale_power_of_two(exponent).to_f64())
}

// Exact nonzero finite x = mantissa * 2^exponent, with |mantissa| in [1,2).
fn binary_parts(x: f64) -> (f64, i32) {
    let (value, correction) = if x.abs() < f64::MIN_POSITIVE {
        (x * 18_014_398_509_481_984.0, -54)
    } else {
        (x, 0)
    };
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + correction;
    let mantissa = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (1023_u64 << 52));
    (mantissa, exponent)
}

#[cfg(test)]
mod tests {
    use super::super::{LorentzError, LorentzState};
    use super::*;

    fn near(actual: f64, expected: f64, scale: f64) {
        assert!(
            (actual - expected).abs() <= 2e-12 * scale.max(1e-300),
            "actual {actual:.17e}, expected {expected:.17e}, scale {scale:.17e}"
        );
    }

    fn phase(omega: f64, time: f64) -> C64 {
        let (sine, cosine) = (omega * time).sin_cos();
        C64::new(cosine, sine)
    }

    fn real_at(tones: &[(f64, C64)], time: f64) -> f64 {
        tones
            .iter()
            .map(|&(omega, amplitude)| (amplitude * phase(omega, time)).re)
            .sum()
    }

    // Independent analytic product-to-sum integral. Both sum and difference
    // frequencies remain, including the DC limit and finite beat windows.
    fn window(omega: f64, start: f64, duration: f64) -> C64 {
        let half_angle = omega * duration * 0.5;
        let sinc = if half_angle == 0.0 {
            1.0
        } else {
            half_angle.sin() / half_angle
        };
        phase(omega, start + 0.5 * duration).scale(duration * sinc)
    }

    fn product_integral(a: &[(f64, C64)], b: &[(f64, C64)], start: f64, duration: f64) -> f64 {
        let mut value = 0.0;
        for &(wa, ca) in a {
            for &(wb, cb) in b {
                value += 0.5
                    * (ca * cb.conj() * window(wa - wb, start, duration)
                        + ca * cb * window(wa + wb, start, duration))
                    .re;
            }
        }
        value
    }

    #[test]
    fn g1_channels_retain_original_ode_phase_and_passive_cycle_heat() {
        let (delta, natural): (f64, f64) = (0.7, 7.0);
        for damping in [0.0, 0.03, 14.0, 40.0] {
            let pole = LorentzPole::try_new(delta, natural, damping).unwrap();
            for omega in [0.0, 0.03, 1.0, 4.0, 7.0, 16.0, 100.0] {
                if damping == 0.0 && omega == natural {
                    continue;
                }
                let transfer = pole.harmonic_transfer(omega).unwrap();
                let p = transfer.restoring.scale(delta.sqrt());
                let j = transfer.kinetic.scale(natural * delta.sqrt());
                let expected_j = C64::new(-p.im, p.re).scale(omega);
                near(j.re, expected_j.re, j.abs());
                near(j.im, expected_j.im, j.abs());
                let ode = p.scale(natural * natural - omega * omega) + j.scale(damping);
                near(ode.re, delta * natural * natural, delta * natural * natural);
                near(ode.im, 0.0, delta * natural * natural);
                let dissipated = 0.5 * transfer.heat.norm_sq();
                near(dissipated, 0.5 * omega * -p.im, dissipated);
                let stored = 0.25 * (transfer.restoring.norm_sq() + transfer.kinetic.norm_sq());
                assert!(stored.is_finite() && stored > 0.0);
                if omega > natural {
                    assert!(p.re < 0.0);
                }
                if damping == 0.0 {
                    assert_eq!(transfer.heat, C64::ZERO);
                }
            }
        }
    }

    #[test]
    fn g1_finite_window_matches_actual_harmonic_propagator() {
        let pole = LorentzPole::try_new(0.7, 7.0, 0.4).unwrap();
        let field = C64::new(1.2, -0.9);
        for omega in [2.0, 7.0, 19.0] {
            let transfer = pole.harmonic_transfer(omega).unwrap();
            let p = (transfer.restoring * field).scale(0.7_f64.sqrt());
            let j = (transfer.kinetic * field).scale(7.0 * 0.7_f64.sqrt());
            let h = transfer.heat * field;
            let duration = 0.83;
            let frame = pole
                .prepare_harmonic(omega, duration)
                .unwrap()
                .advance(
                    LorentzState {
                        polarization_c_m2: [p.re],
                        current_a_m2: [j.re],
                    },
                    [field.re],
                    [-field.im],
                )
                .unwrap();
            near(
                frame.state.polarization_c_m2[0],
                (p * phase(omega, duration)).re,
                p.abs(),
            );
            near(
                frame.state.current_a_m2[0],
                (j * phase(omega, duration)).re,
                j.abs(),
            );
            near(
                frame.heat_j_m3,
                product_integral(&[(omega, h)], &[(omega, h)], 0.0, duration),
                frame.heat_j_m3,
            );
            near(
                frame.supplied_work_j_m3,
                product_integral(&[(omega, field)], &[(omega, j)], 0.0, duration),
                frame.supplied_work_j_m3.abs().max(1.0),
            );
        }
    }

    #[test]
    fn g1_coherent_two_tone_window_retains_beat_storage_and_heat() {
        let pole = LorentzPole::try_new(0.7, 7.0, 0.4).unwrap();
        let fields = [(6.8, C64::new(1.2, -0.9)), (7.15, C64::new(-0.3, 0.8))];
        let mut restoring = [(0.0, C64::ZERO); 2];
        let mut kinetic = restoring;
        let mut current = restoring;
        let mut heat = restoring;
        for (i, &(omega, field)) in fields.iter().enumerate() {
            let transfer = pole.harmonic_transfer(omega).unwrap();
            restoring[i] = (omega, transfer.restoring * field);
            kinetic[i] = (omega, transfer.kinetic * field);
            current[i] = (omega, kinetic[i].1.scale(7.0 * 0.7_f64.sqrt()));
            heat[i] = (omega, transfer.heat * field);
        }
        let energy =
            |time| 0.5 * (real_at(&restoring, time).powi(2) + real_at(&kinetic, time).powi(2));
        let (start, duration) = (0.17, 3.21);
        let work = product_integral(&fields, &current, start, duration);
        let dissipated = product_integral(&heat, &heat, start, duration);
        assert!(dissipated > 0.0);
        near(
            energy(start + duration) - energy(start) + dissipated,
            work,
            work.abs(),
        );
        let incoherent_heat: f64 = heat
            .iter()
            .map(|&tone| product_integral(&[tone], &[tone], start, duration))
            .sum();
        assert!((dissipated - incoherent_heat).abs() > 0.1 * dissipated);
        let first = 0.71;
        near(
            dissipated,
            product_integral(&heat, &heat, start, first)
                + product_integral(&heat, &heat, start + first, duration - first),
            dissipated,
        );
        near(
            work,
            product_integral(&fields, &current, start, first)
                + product_integral(&fields, &current, start + first, duration - first),
            work.abs(),
        );
    }

    #[test]
    fn g1_extreme_normalization_preserves_representable_channels() {
        let transfer = LorentzPole::try_new(1e300, 1e-50, 1e-100)
            .unwrap()
            .harmonic_transfer(1e100)
            .unwrap();
        near(transfer.restoring.re, -1e-150, 1e-150);
        assert_eq!(transfer.restoring.im, 0.0);
        near(transfer.kinetic.re, 1e-200, 1e-200);
        near(transfer.kinetic.im, -1.0, 1.0);
        near(transfer.heat.re, 1e-250, 1e-250);
        near(transfer.heat.im, -1e-50, 1e-50);

        let transfer = LorentzPole::try_new(1e300, 1e100, 1e100)
            .unwrap()
            .harmonic_transfer(1e100)
            .unwrap();
        near(transfer.restoring.im, -1e150, 1e150);
        near(transfer.kinetic.re, 1e150, 1e150);
        near(transfer.heat.re, 1e200, 1e200);

        let transfer = LorentzPole::try_new(1e-100, 1e100, 1e200)
            .unwrap()
            .harmonic_transfer(1e-100)
            .unwrap();
        assert_eq!(transfer.kinetic.re, 0.0);
        near(transfer.kinetic.im, 1e-250, 1e-250);
        near(transfer.heat.re, 1e-250, 1e-250);
        near(transfer.heat.im, 1e-150, 1e-150);

        let omega = f64::from_bits(1);
        let transfer = LorentzPole::try_new(1e100, 1.0, 0.0)
            .unwrap()
            .harmonic_transfer(omega)
            .unwrap();
        near(transfer.kinetic.im, omega * 1e50, omega * 1e50);
    }

    #[test]
    fn g1_static_lossless_and_invalid_frequency_contracts() {
        let pole = LorentzPole::try_new(0.7, 7.0, 0.0).unwrap();
        let transfer = pole.harmonic_transfer(0.0).unwrap();
        assert_eq!(transfer.restoring, C64::from_re(det::sqrt(0.7)));
        assert_eq!(transfer.kinetic, C64::ZERO);
        assert_eq!(transfer.heat, C64::ZERO);
        assert_eq!(
            pole.harmonic_transfer(7.0),
            Err(LorentzError::SingularResonance)
        );
        for omega in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                pole.harmonic_transfer(omega),
                Err(LorentzError::Invalid(_))
            ));
        }
        assert!(
            LorentzPole::try_new(1e300, 1.0, 1e-300)
                .unwrap()
                .harmonic_transfer(1.0)
                .is_err()
        );
        assert_eq!(pole.harmonic_transfer(9.0), pole.harmonic_transfer(9.0));
    }
}
