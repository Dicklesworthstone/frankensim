//! Implicit characteristic termination using the original digital impedance.
//!
//! For pressure waves a (incident) and b (reflected), p=a+b and U=(a-b)/Zc.
//! The accepted impedance state gives p=d*U+h at the CURRENT sample, hence
//! U=(2*a-h)/(Zc+d). This closes the direct feedthrough algebraically instead
//! of introducing a sample of delay or fitting a second reflection operator.

use super::{DigitalFilter, DigitalFilterState, DiscretizeError, check_digital_passivity};
use crate::passivity::{GlobalPassivityBudget, GlobalPassivityReport, GlobalPassivityVerdict};
use crate::scattering::ScatteringResponse;

/// Refusal of impedance admission or an atomic boundary sample.
#[derive(Debug, Clone, PartialEq)]
pub enum ImpedanceScatteringError {
    /// Invalid characteristic impedance or unrepresentable boundary arithmetic.
    Invalid(&'static str),
    /// The original digital impedance refused coefficients or a sample.
    Discrete(DiscretizeError),
    /// The actual digital coefficients were not globally certified passive.
    Uncertified(GlobalPassivityReport),
}

impl core::fmt::Display for ImpedanceScatteringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "impedance scattering port: {reason}"),
            Self::Discrete(error) => write!(f, "impedance scattering port: {error}"),
            Self::Uncertified(report) => write!(
                f,
                "impedance scattering port: digital passivity {:?} ({:?})",
                report.verdict, report.unresolved
            ),
        }
    }
}
impl std::error::Error for ImpedanceScatteringError {}
impl From<DiscretizeError> for ImpedanceScatteringError {
    fn from(error: DiscretizeError) -> Self {
        Self::Discrete(error)
    }
}

/// One accepted pressure/volume-flow characteristic boundary sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImpedanceScatteringFrame {
    /// Total pressure at the load [Pa].
    pub pressure_pa: f64,
    /// Volume flow INTO the load [m³/s].
    pub volume_flow_m3_s: f64,
    /// Reflected pressure wave leaving the load [Pa].
    pub reflected_pressure_pa: f64,
    /// Signed supplied work p*U*T [J]. Negative means energy returned by the
    /// load. This is neither resistive loss nor a stored-energy estimate.
    pub supplied_work_j: f64,
}

/// Same-sample reflection from a globally passive digital impedance.
///
/// The supplied [`DigitalFilter`] is the pressure/volume-flow impedance in
/// Pa·s/m³, NOT a reflection filter. Its coefficients and DF-II memories remain
/// the only constitutive operator. No reflection refit, gain clamp, propagation
/// delay or extra time integrator is introduced. Preparation certifies the
/// actual stored coefficients; stepping and [`Self::response`] allocate nothing.
///
/// The coefficient certificate is not a rounded runtime-work proof or evidence
/// of the physical validity of fitted data. The original lossless Tustin
/// differentiator is supported; a DC integrator not admitted by the digital
/// checker is refused rather than approximated with leakage.
#[derive(Debug, Clone)]
pub struct ImpedanceScatteringPort {
    filter: DigitalFilter,
    state: DigitalFilterState,
    certificate: GlobalPassivityReport,
    characteristic_impedance: f64,
    scale: f64,
    normalized_denominator: f64,
    reflection_direct: f64,
    history_weight: f64,
}

impl ImpedanceScatteringPort {
    /// Admit a pressure/volume-flow impedance and start with zero memory.
    ///
    /// The reference impedance must be positive and finite. An active or
    /// unresolved full-circle certificate refuses construction. Admission
    /// rechecks the supplied filter, including a filter from a public mutable
    /// realization report; it never trusts a stale certificate.
    pub fn new(
        filter: DigitalFilter,
        characteristic_impedance: f64,
        budget: GlobalPassivityBudget,
    ) -> Result<Self, ImpedanceScatteringError> {
        if !(characteristic_impedance > 0.0 && characteristic_impedance.is_finite()) {
            return Err(ImpedanceScatteringError::Invalid(
                "characteristic impedance must be positive and finite",
            ));
        }
        let certificate = check_digital_passivity(&filter, budget)?;
        if certificate.verdict != GlobalPassivityVerdict::Certified {
            return Err(ImpedanceScatteringError::Uncertified(certificate));
        }
        let direct = filter
            .sections
            .iter()
            .fold(filter.direct, |sum, s| sum + s.b[0]);
        if !(direct >= 0.0 && direct.is_finite()) {
            return Err(ImpedanceScatteringError::Invalid(
                "instantaneous impedance is not representable and nonnegative",
            ));
        }
        // Scaling avoids overflowing d+Zc even for a valid large impedance.
        let scale = direct.max(characteristic_impedance);
        let normalized_load = direct / scale;
        let normalized_reference = characteristic_impedance / scale;
        let normalized_denominator = normalized_load + normalized_reference;
        let state = filter.zero_state();
        Ok(Self {
            filter,
            state,
            certificate,
            characteristic_impedance,
            scale,
            normalized_denominator,
            reflection_direct: (normalized_load - normalized_reference) / normalized_denominator,
            history_weight: normalized_reference / normalized_denominator,
        })
    }

    /// Immutable coefficients actually stepped by the port.
    #[must_use]
    pub fn filter(&self) -> &DigitalFilter {
        &self.filter
    }

    /// Admission certificate bound to the immutable owned coefficients.
    #[must_use]
    pub fn certificate(&self) -> &GlobalPassivityReport {
        &self.certificate
    }

    /// Characteristic pressure/volume-flow impedance [Pa·s/m³].
    #[must_use]
    pub fn characteristic_impedance(&self) -> f64 {
        self.characteristic_impedance
    }

    /// Affine reflected pressure for the next accepted incident sample.
    /// A physical source/junction can include this relation in its implicit
    /// solve, then call [`Self::push`] once with the solved incident pressure.
    pub fn response(&self) -> Result<ScatteringResponse, ImpedanceScatteringError> {
        let history = self.filter.validate_step(&self.state, 0.0)? * self.history_weight;
        if !history.is_finite() {
            return Err(ImpedanceScatteringError::Invalid(
                "reflection history is not finite",
            ));
        }
        Ok(ScatteringResponse {
            direct: self.reflection_direct,
            history,
        })
    }

    /// Solve the current boundary and commit only after ALL outputs and the
    /// original filter's candidate memories have been checked. Refusal leaves
    /// accepted history unchanged. No heap allocation occurs in this method.
    pub fn push(
        &mut self,
        incident_pressure_pa: f64,
    ) -> Result<ImpedanceScatteringFrame, ImpedanceScatteringError> {
        if !incident_pressure_pa.is_finite() {
            return Err(ImpedanceScatteringError::Invalid(
                "incident pressure is not finite",
            ));
        }
        let history = self.filter.validate_step(&self.state, 0.0)?;
        let volume_flow_m3_s = (2.0 * (incident_pressure_pa / self.scale) - history / self.scale)
            / self.normalized_denominator;
        let pressure_pa = self.filter.validate_step(&self.state, volume_flow_m3_s)?;
        let reflected_pressure_pa = pressure_pa - incident_pressure_pa;
        let supplied_work_j = pressure_pa * (volume_flow_m3_s * self.filter.t_s);
        if !reflected_pressure_pa.is_finite() || !supplied_work_j.is_finite() {
            return Err(ImpedanceScatteringError::Invalid(
                "boundary reflection/work is not finite",
            ));
        }
        self.filter.commit_step(&mut self.state, volume_flow_m3_s);
        Ok(ImpedanceScatteringFrame {
            pressure_pa,
            volume_flow_m3_s,
            reflected_pressure_pa,
            supplied_work_j,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::discretize::{bilinear, reflectance};
    use crate::model::{PoleTerm, RationalModel};
    use fs_math::{c64::C64, det};

    fn port(model: &RationalModel, reference: f64) -> ImpedanceScatteringPort {
        ImpedanceScatteringPort::new(
            bilinear(model, 1.0 / 48_000.0, 0.0).unwrap(),
            reference,
            GlobalPassivityBudget::default(),
        )
        .unwrap()
    }

    #[test]
    fn resistance_limits_have_current_sample_reflection_and_wave_work() {
        for resistance in [0.0, 0.25, 1.0, 4.0, 1000.0] {
            let model = RationalModel {
                terms: vec![],
                d: resistance,
                e: 0.0,
            };
            let mut load = port(&model, 1.0);
            let response = load.response().unwrap();
            let frame = load.push(0.75).unwrap();
            let expected = 0.75 * (resistance - 1.0) / (resistance + 1.0);
            assert!((frame.reflected_pressure_pa - expected).abs() < 1e-14);
            assert!((response.incoming(0.75).unwrap() - expected).abs() < 1e-14);
            assert!((frame.pressure_pa - (0.75 + frame.reflected_pressure_pa)).abs() < 1e-14);
            let wave_work = (0.75 * 0.75 - frame.reflected_pressure_pa.powi(2)) / 48_000.0;
            assert!((frame.supplied_work_j - wave_work).abs() < 1e-18);
            assert!(frame.supplied_work_j >= 0.0);
        }
    }

    #[test]
    fn reactive_reflection_preserves_original_filter_and_returns_stored_energy() {
        let model = RationalModel {
            terms: vec![],
            d: 0.0,
            e: 0.01,
        };
        let mut load = port(&model, 1.0);
        let filter = load.filter().clone();
        let mut original = filter.zero_state();
        let mut returned = false;
        let mut work = 0.0;
        for n in 0..2000 {
            let incident = if n < 32 { 1.0 } else { 0.0 };
            let response = load.response().unwrap();
            let frame = load.push(incident).unwrap();
            let original_pressure = filter.step(&mut original, frame.volume_flow_m3_s).unwrap();
            assert_eq!(frame.pressure_pa, original_pressure);
            assert!(
                (response.incoming(incident).unwrap() - frame.reflected_pressure_pa).abs() < 1e-12
            );
            let characteristic_flow =
                (incident - frame.reflected_pressure_pa) / load.characteristic_impedance();
            assert!((characteristic_flow - frame.volume_flow_m3_s).abs() < 1e-12);
            returned |= frame.supplied_work_j < -1e-10;
            work += frame.supplied_work_j;
            assert!(
                work >= -1e-14,
                "zero-state lossless load created net energy"
            );
        }
        assert!(
            returned,
            "reactive history must return energy after excitation ends"
        );
    }

    #[test]
    fn harmonic_reflection_matches_complex_impedance_without_sample_lag() {
        let model = RationalModel {
            terms: vec![PoleTerm::Real {
                pole: -600.0,
                residue: 900.0,
            }],
            d: 0.25,
            e: 0.0001,
        };
        for cycles in [17, 103, 997, 3021] {
            let mut load = port(&model, 2.0);
            let period = 8192;
            let theta = core::f64::consts::TAU * f64::from(cycles) / f64::from(period);
            let expected = reflectance(load.filter().eval(theta * 48_000.0), 2.0);
            let mut response = C64::ZERO;
            for n in 0..(period * 4) {
                let phase = theta * f64::from(n);
                let frame = load.push(det::cos(phase)).unwrap();
                if n >= period * 3 {
                    response = response
                        + C64::new(det::cos(phase), -det::sin(phase))
                            .scale(frame.reflected_pressure_pa * 2.0 / f64::from(period));
                }
            }
            assert!(
                (response - expected).abs() < 2e-10,
                "{cycles}: {response:?} != {expected:?}"
            );
            assert!(response.abs() <= 1.0 + 1e-12);
        }
    }

    #[test]
    fn refused_samples_leave_original_history_unchanged() {
        let model = RationalModel {
            terms: vec![PoleTerm::Real {
                pole: -40.0,
                residue: 20.0,
            }],
            d: 1.0,
            e: 0.0,
        };
        let mut load = port(&model, 1.0);
        load.push(1.0).unwrap();
        let mut untouched = load.clone();
        for invalid in [f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(load.push(invalid).is_err());
            assert_eq!(load.state, untouched.state);
            assert_eq!(load.response().unwrap(), untouched.response().unwrap());
        }
        for incident in [0.0, -0.3, 0.8, 0.0] {
            assert_eq!(
                load.push(incident).unwrap(),
                untouched.push(incident).unwrap()
            );
        }
    }

    #[test]
    fn actual_active_or_unresolved_coefficients_are_refused() {
        let filter = DigitalFilter {
            sections: vec![],
            direct: -1.0,
            t_s: 0.001,
            prewarp: 0.0,
        };
        assert!(matches!(
            ImpedanceScatteringPort::new(filter, 1.0, GlobalPassivityBudget::default()),
            Err(ImpedanceScatteringError::Uncertified(_))
        ));
        let model = RationalModel {
            terms: vec![PoleTerm::Real {
                pole: -40.0,
                residue: 20.0,
            }],
            d: 1.0,
            e: 0.0,
        };
        let filter = bilinear(&model, 0.001, 0.0).unwrap();
        assert!(matches!(
            ImpedanceScatteringPort::new(
                filter.clone(),
                1.0,
                GlobalPassivityBudget {
                    max_order: 0,
                    ..GlobalPassivityBudget::default()
                }
            ),
            Err(ImpedanceScatteringError::Uncertified(_))
        ));
        for reference in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                ImpedanceScatteringPort::new(
                    filter.clone(),
                    reference,
                    GlobalPassivityBudget::default()
                )
                .is_err()
            );
        }
    }

    #[test]
    fn scaled_characteristic_solve_does_not_overflow_large_matching_impedances() {
        let filter = DigitalFilter {
            sections: vec![],
            direct: f64::MAX,
            t_s: 0.001,
            prewarp: 0.0,
        };
        let mut load =
            ImpedanceScatteringPort::new(filter, f64::MAX, GlobalPassivityBudget::default())
                .unwrap();
        assert_eq!(load.response().unwrap().direct, 0.0);
        let frame = load.push(1.0).unwrap();
        assert!(frame.reflected_pressure_pa.abs() < 1e-14);
        assert!(frame.volume_flow_m3_s > 0.0);
    }
}
