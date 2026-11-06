//! Free-jet screening composed from the original cylindrical capillary owner.
//!
//! The caller supplies a real cylindrical radius, volume flow, flight length,
//! physical wavenumber band and common initial eigenmode amplitude. A planar
//! pour's connected-component count supplies none of these quantities.
//! Continuity gives U=Q/(pi R^2); the native temporal growth is observed over
//! the material transit time L/U. The report keeps logarithmic amplitudes so
//! strongly amplified candidates remain comparable without exponent overflow.
//!
//! This is a uniform, un-stretched cylinder with a prescribed
//! infinitesimal normal mode, not an inlet-forced spatial stability problem,
//! a nozzle boundary-value solution or nonlinear pinch-off/drop-size dynamics.
//! Liquid/exterior assumptions and the dispersion remain in fs-lbm::jet.
//! The selected owner admits either the full inviscid cylindrical relation or
//! an explicitly bounded viscous slender-jet relation. Source-bound owners
//! retain the original material alongside the study. No runtime dispatch or
//! duplicate free-flight calculation is introduced.

use fs_lbm::jet::{CapillaryJet, JetBandMaximum, JetDispersion, JetError};

/// An authored SI free-flight study; no liquid or disturbance defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeJetStudy<J = CapillaryJet> {
    /// Original cylindrical physics owner and its explicit approximation.
    pub jet: J,
    /// Positive liquid volume flow [m^3/s], through the declared cylinder.
    pub volume_flow_m3_s: f64,
    /// Nonnegative distance travelled by the material [m].
    pub flight_length_m: f64,
    /// Closed physical wavenumber band [rad/m], within the model's domain.
    pub wavenumbers_rad_m: [f64; 2],
    /// Common initial amplitude/radius of each normal mode, strictly in (0,1).
    /// This is not an arbitrary displacement with unspecified initial velocity.
    pub initial_relative_amplitude: f64,
    /// Authored amplitude/radius threshold, between the initial value and 1.
    /// Even below 1, linear validity requires small disturbances. Crossing this
    /// threshold reports a linear extrapolation, never a resolved breakup.
    pub threshold_relative_amplitude: f64,
}

/// Worst linear amplification across the entire requested wavenumber band.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeJetReport {
    /// Native continuous-band maximizer; no sampled wavelengths are skipped.
    pub fastest: JetBandMaximum,
    /// Velocity fixed by volume continuity [m/s].
    pub mean_velocity_m_s: f64,
    /// Authored distance divided by the same velocity [s].
    pub residence_time_s: f64,
    /// Natural log of the modal amplification factor: growth times transit time.
    pub log_amplification: f64,
    /// Natural log of the final amplitude/radius, without overflowing exp.
    pub log_relative_amplitude: f64,
    /// Linear time to the authored threshold [s]; None for a neutral/oscillatory
    /// band below threshold. An initially reached threshold returns Some(0).
    pub linear_threshold_time_s: Option<f64>,
    /// Same linear threshold expressed as advected distance [m].
    pub linear_threshold_distance_m: Option<f64>,
    /// Whether the authored transit reaches the linear threshold.
    pub threshold_reached: bool,
}

impl<J: JetDispersion> FreeJetStudy<J> {
    /// Evaluate one uniform cylindrical segment without allocating or sampling.
    ///
    /// No implied material calibration, gravity stretching, surrounding-gas
    /// inertia, arbitrary disturbance evolution or nonlinear breakup is added.
    /// Invalid inputs and unrepresentable derived quantities refuse the study.
    pub fn evaluate(self) -> Result<FreeJetReport, JetError> {
        if !self.volume_flow_m3_s.is_finite() || self.volume_flow_m3_s <= 0.0 {
            return Err(JetError::InvalidInput("jet volume flow"));
        }
        if !self.flight_length_m.is_finite() || self.flight_length_m < 0.0 {
            return Err(JetError::InvalidInput("jet flight length"));
        }
        let initial = self.initial_relative_amplitude;
        let threshold = self.threshold_relative_amplitude;
        if !initial.is_finite() || initial <= 0.0 || initial >= 1.0 {
            return Err(JetError::InvalidInput("initial relative jet amplitude"));
        }
        if !threshold.is_finite() || threshold < initial || threshold >= 1.0 {
            return Err(JetError::InvalidInput("relative jet amplitude threshold"));
        }
        let fastest = self.jet.fastest_mode(self.wavenumbers_rad_m)?;
        let mean_velocity_m_s = fs_math::det::exp(
            fs_math::det::ln(self.volume_flow_m3_s)
                - fs_math::det::ln(core::f64::consts::PI)
                - 2.0 * fs_math::det::ln(self.jet.radius_m()),
        );
        represented(mean_velocity_m_s, true)?;
        let residence_time_s = self.flight_length_m / mean_velocity_m_s;
        represented(residence_time_s, self.flight_length_m > 0.0)?;
        let growth = fastest.mode.growth_per_s;
        let log_amplification = growth * residence_time_s;
        represented(log_amplification, growth > 0.0 && residence_time_s > 0.0)?;
        let log_initial = fs_math::det::ln(initial);
        let log_relative_amplitude = log_initial + log_amplification;
        represented(log_relative_amplitude, false)?;
        let needed = fs_math::det::ln(threshold) - log_initial;
        represented(needed, threshold > initial)?;
        let linear_threshold_time_s = if threshold == initial {
            Some(0.0)
        } else if growth > 0.0 {
            let time = needed / growth;
            represented(time, true)?;
            Some(time)
        } else {
            None
        };
        let linear_threshold_distance_m = match linear_threshold_time_s {
            Some(time) => {
                let distance = time * mean_velocity_m_s;
                represented(distance, time > 0.0)?;
                Some(distance)
            }
            None => None,
        };
        Ok(FreeJetReport {
            fastest,
            mean_velocity_m_s,
            residence_time_s,
            log_amplification,
            log_relative_amplitude,
            linear_threshold_time_s,
            linear_threshold_distance_m,
            threshold_reached: log_amplification >= needed,
        })
    }
}

fn represented(value: f64, nonzero: bool) -> Result<(), JetError> {
    if !value.is_finite() || (nonzero && value == 0.0) {
        Err(JetError::NonRepresentable)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn study() -> FreeJetStudy {
        FreeJetStudy {
            jet: CapillaryJet::inviscid(0.001, 1000.0, 0.072).unwrap(),
            volume_flow_m3_s: core::f64::consts::PI * 1e-6 * 0.5,
            flight_length_m: 0.05,
            wavenumbers_rad_m: [0.0, 2000.0],
            initial_relative_amplitude: 1e-4,
            threshold_relative_amplitude: 0.01,
        }
    }

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 3e-12 * expected.abs().max(f64::MIN_POSITIVE),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g1_cylindrical_growth_reaches_the_actual_flight_score() {
        let report = study().evaluate().unwrap();
        near(report.mean_velocity_m_s, 0.5);
        near(report.residence_time_s, 0.1);
        let expected_growth = 0.343_338_701_861_299_3 * 72_000.0_f64.sqrt();
        near(report.fastest.mode.growth_per_s, expected_growth);
        near(report.log_amplification, 0.1 * expected_growth);
        near(
            report.linear_threshold_time_s.unwrap(),
            100.0_f64.ln() / expected_growth,
        );
        near(
            report.linear_threshold_distance_m.unwrap(),
            0.5 * report.linear_threshold_time_s.unwrap(),
        );
        assert!(report.threshold_reached);
    }

    #[test]
    fn g3_flow_and_length_change_transit_without_refitting_the_mode() {
        let base = study();
        let report = base.evaluate().unwrap();
        let faster = FreeJetStudy {
            volume_flow_m3_s: 4.0 * base.volume_flow_m3_s,
            ..base
        }
        .evaluate()
        .unwrap();
        assert_eq!(faster.fastest, report.fastest);
        near(faster.log_amplification, report.log_amplification / 4.0);
        near(
            faster.linear_threshold_distance_m.unwrap(),
            4.0 * report.linear_threshold_distance_m.unwrap(),
        );
        assert!(!faster.threshold_reached);
        let longer = FreeJetStudy {
            flight_length_m: 2.0 * base.flight_length_m,
            ..base
        }
        .evaluate()
        .unwrap();
        near(longer.log_amplification, 2.0 * report.log_amplification);
    }

    #[test]
    fn g0_neutral_stable_and_zero_flight_have_no_invented_breakup() {
        for candidate in [
            FreeJetStudy {
                wavenumbers_rad_m: [1200.0, 2000.0],
                ..study()
            },
            FreeJetStudy {
                jet: CapillaryJet::inviscid(0.001, 1000.0, 0.0).unwrap(),
                ..study()
            },
        ] {
            let report = candidate.evaluate().unwrap();
            assert_eq!(report.log_amplification, 0.0);
            assert_eq!(report.linear_threshold_time_s, None);
            assert_eq!(report.linear_threshold_distance_m, None);
            assert!(!report.threshold_reached);
        }
        let no_flight = FreeJetStudy {
            flight_length_m: 0.0,
            ..study()
        }
        .evaluate()
        .unwrap();
        assert_eq!(no_flight.log_amplification, 0.0);
        assert!(!no_flight.threshold_reached);
        let already_at_threshold = FreeJetStudy {
            threshold_relative_amplitude: study().initial_relative_amplitude,
            ..study()
        }
        .evaluate()
        .unwrap();
        assert_eq!(already_at_threshold.linear_threshold_time_s, Some(0.0));
        assert!(already_at_threshold.threshold_reached);
    }

    #[test]
    fn g0_large_scores_remain_logarithmic_and_bad_studies_refuse() {
        let large = FreeJetStudy {
            flight_length_m: 50.0,
            ..study()
        }
        .evaluate()
        .unwrap();
        assert!(large.log_relative_amplitude > 1000.0);
        assert!(large.log_relative_amplitude.is_finite());
        for candidate in [
            FreeJetStudy { volume_flow_m3_s: 0.0, ..study() },
            FreeJetStudy { flight_length_m: -1.0, ..study() },
            FreeJetStudy { initial_relative_amplitude: 0.0, ..study() },
            FreeJetStudy { threshold_relative_amplitude: 1.0, ..study() },
            FreeJetStudy { threshold_relative_amplitude: 1e-6, ..study() },
            FreeJetStudy { wavenumbers_rad_m: [2000.0, 100.0], ..study() },
        ] {
            assert!(candidate.evaluate().is_err());
        }
    }

    #[test]
    fn g1_viscous_model_reduces_the_actual_flight_amplification() {
        use fs_lbm::jet::viscous::ViscousSlenderJet;

        let base = study();
        let inviscid = FreeJetStudy {
            wavenumbers_rad_m: [0.0, 300.0],
            ..base
        }.evaluate().unwrap();
        let viscous_study = FreeJetStudy {
            jet: ViscousSlenderJet::new(0.001, 1000.0, 0.072, 0.001, 0.3).unwrap(),
            volume_flow_m3_s: base.volume_flow_m3_s,
            flight_length_m: base.flight_length_m,
            wavenumbers_rad_m: [0.0, 300.0],
            initial_relative_amplitude: base.initial_relative_amplitude,
            threshold_relative_amplitude: base.threshold_relative_amplitude,
        };
        let viscous = viscous_study.evaluate().unwrap();
        near(viscous.mean_velocity_m_s, inviscid.mean_velocity_m_s);
        near(viscous.residence_time_s, inviscid.residence_time_s);
        let omega = (0.072_f64 / (1000.0 * 0.001_f64.powi(3))).sqrt();
        let oh = 0.001 / (0.072_f64 * 0.001 / 1000.0).sqrt();
        let expected_growth = omega / (8.0_f64.sqrt() + 6.0 * oh);
        near(viscous.fastest.mode.growth_per_s, expected_growth);
        near(viscous.log_amplification, expected_growth * 0.1);
        assert!(viscous.log_amplification < inviscid.log_amplification);
        assert!(
            viscous.linear_threshold_distance_m.unwrap()
                > inviscid.linear_threshold_distance_m.unwrap()
        );
        assert!(!viscous.threshold_reached);
        assert!(FreeJetStudy {
            wavenumbers_rad_m: [0.0, 400.0],
            ..viscous_study
        }.evaluate().is_err());
    }
}
