//! Static adequacy of the exact reduced image consumed by a one-mode valve.
//!
//! The pressure loading vector `l` is force per pascal. The observation `g`
//! is mean face opening per displacement. With unit-opening modes psi_i,
//! `A_i = psi_i^T l`, pressure volume compliance is sum A_i^2/k_i,
//! and opening per pascal is sum A_i/k_i. The exact comparators are
//! `l^T K^-1 l`, `g^T K^-1 l`, and `g^T K^-1 g` from the SAME solid pencil.
//! This admits static linear truncation only, not transient, contact, high-
//! frequency, pressure-field, damping-spectrum or full-instrument accuracy.

use fs_blake3::DomainHasher;
use super::{ReduceError, ReducedMode};

/// Explicit acceptance criteria for the single mode that will actually play.
/// No default promotes an arbitrary reduction tolerance to physical accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SingleModeValveAdmission {
    /// Required `(1/k_0) / (g^T K^-1 g)`, in `(0, 1]`.
    pub minimum_opening_compliance_fraction: f64,
    /// Required `(A_0^2/k_0) / (l^T K^-1 l)`, in `(0, 1]`.
    pub minimum_pressure_compliance_fraction: f64,
    /// Required relative agreement of `A_0/k_0` with `g^T K^-1 l`.
    /// Must be finite in `[0, 1)`. Cancellation or a non-closing exact static
    /// response cannot be hidden by taking an absolute value of that response.
    pub maximum_static_opening_relative_error: f64,
}

/// Immutable static responses, preserving signed pressure/opening coupling.
/// A pressure load orthogonal to opening may still store elastic energy;
/// the two response measures are deliberately not conflated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PressureCompliance {
    exact_mean_force_m_per_n: f64,
    exact_volume_m3_per_pa: f64,
    exact_opening_m_per_pa: f64,
    retained_volume_m3_per_pa: f64,
    retained_opening_m_per_pa: f64,
    first_mean_force_m_per_n: f64,
    first_volume_m3_per_pa: f64,
    first_opening_m_per_pa: f64,
}

fn invalid(what: &'static str) -> ReduceError { ReduceError::Invalid { what } }

impl PressureCompliance {
    pub(super) fn from_modes(
        modes: &[ReducedMode],
        areas_m2: &[f64],
        exact_mean_force_m_per_n: f64,
        exact_volume_m3_per_pa: f64,
        exact_opening_m_per_pa: f64,
    ) -> Result<Self, ReduceError> {
        if modes.is_empty() || modes.len() != areas_m2.len()
            || !exact_mean_force_m_per_n.is_finite() || exact_mean_force_m_per_n <= 0.0
            || !exact_volume_m3_per_pa.is_finite() || exact_volume_m3_per_pa < 0.0
            || !exact_opening_m_per_pa.is_finite()
        {
            return Err(invalid("pressure compliance requires aligned modes and finite physical static responses"));
        }
        let mut result = Self {
            exact_mean_force_m_per_n,
            exact_volume_m3_per_pa,
            exact_opening_m_per_pa,
            retained_volume_m3_per_pa: 0.0,
            retained_opening_m_per_pa: 0.0,
            first_mean_force_m_per_n: 0.0,
            first_volume_m3_per_pa: 0.0,
            first_opening_m_per_pa: 0.0,
        };
        for (index, (mode, &area)) in modes.iter().zip(areas_m2).enumerate() {
            if !mode.stiffness_n_m.is_finite() || mode.stiffness_n_m <= 0.0 || !area.is_finite() {
                return Err(invalid("pressure compliance requires positive finite stiffness and finite area"));
            }
            let opening = area / mode.stiffness_n_m;
            let volume = area * opening;
            result.retained_opening_m_per_pa += opening;
            result.retained_volume_m3_per_pa += volume;
            if index == 0 {
                result.first_mean_force_m_per_n = 1.0 / mode.stiffness_n_m;
                result.first_opening_m_per_pa = opening;
                result.first_volume_m3_per_pa = volume;
            }
        }
        if !result.values().iter().all(|v| v.is_finite())
            || result.first_mean_force_m_per_n <= 0.0
        {
            return Err(invalid("pressure compliance reduction left the finite representable domain"));
        }
        Ok(result)
    }

    /// Exact mean opening per unit force conjugate to mean opening [m/N].
    #[must_use]
    pub const fn exact_mean_force_m_per_n(&self) -> f64 { self.exact_mean_force_m_per_n }
    /// Exact pressure-work volume compliance [m^3/Pa]. Zero denotes no free load.
    #[must_use]
    pub const fn exact_volume_m3_per_pa(&self) -> f64 { self.exact_volume_m3_per_pa }
    /// Exact signed mean opening per pascal toward opening [m/Pa].
    #[must_use]
    pub const fn exact_opening_m_per_pa(&self) -> f64 { self.exact_opening_m_per_pa }
    /// Pressure-work compliance represented by ALL exported modes [m^3/Pa].
    #[must_use]
    pub const fn retained_volume_m3_per_pa(&self) -> f64 { self.retained_volume_m3_per_pa }
    /// Signed mean opening represented by ALL exported modes [m/Pa].
    #[must_use]
    pub const fn retained_opening_m_per_pa(&self) -> f64 { self.retained_opening_m_per_pa }
    /// Pressure-work compliance of the mode-zero playing image [m^3/Pa].
    #[must_use]
    pub const fn first_volume_m3_per_pa(&self) -> f64 { self.first_volume_m3_per_pa }
    /// Signed mean opening of the mode-zero playing image [m/Pa].
    #[must_use]
    pub const fn first_opening_m_per_pa(&self) -> f64 { self.first_opening_m_per_pa }
    /// Fraction of exact opening-force compliance carried by the playing mode.
    /// A numerical value is not an interval-certified bound; it is not clamped.
    #[must_use]
    pub fn first_opening_compliance_fraction(&self) -> f64 {
        self.first_mean_force_m_per_n / self.exact_mean_force_m_per_n
    }
    /// Fraction of exact pressure compliance; undefined for a zero free load.
    #[must_use]
    pub fn first_pressure_compliance_fraction(&self) -> Option<f64> {
        (self.exact_volume_m3_per_pa > 0.0)
            .then(|| self.first_volume_m3_per_pa / self.exact_volume_m3_per_pa)
    }
    /// Relative error against the exact signed pressure/opening response.
    /// Undefined when that response vanishes, rather than divided by an epsilon.
    #[must_use]
    pub fn first_static_opening_relative_error(&self) -> Option<f64> {
        (self.exact_opening_m_per_pa != 0.0).then(|| {
            ((self.first_opening_m_per_pa - self.exact_opening_m_per_pa)
                / self.exact_opening_m_per_pa).abs()
        })
    }

    /// Check the frozen, linear one-mode approximation without changing it.
    ///
    /// # Errors
    /// Invalid policies, non-closing pressure orientation, negligible coupling
    /// or inadequate compliance refuse; no area/stiffness rescaling repairs it.
    pub fn admit_single_mode(&self, admission: SingleModeValveAdmission) -> Result<(), ReduceError> {
        for floor in [admission.minimum_opening_compliance_fraction,
            admission.minimum_pressure_compliance_fraction]
        {
            if !floor.is_finite() || floor <= 0.0 || floor > 1.0 {
                return Err(invalid("single-mode compliance floors must lie in (0, 1]"));
            }
        }
        let tolerance = admission.maximum_static_opening_relative_error;
        if !tolerance.is_finite() || !(0.0..1.0).contains(&tolerance) {
            return Err(invalid("single-mode static opening tolerance must lie in [0, 1)"));
        }
        if self.exact_opening_m_per_pa <= 0.0 || self.first_opening_m_per_pa <= 0.0 {
            return Err(invalid("single-mode beating reed requires a positive exact and reduced closing compliance"));
        }
        let opening_fraction = self.first_opening_compliance_fraction();
        if !opening_fraction.is_finite() || opening_fraction < admission.minimum_opening_compliance_fraction {
            return Err(invalid("playing mode carries insufficient opening-force compliance"));
        }
        let pressure_fraction = self.first_pressure_compliance_fraction()
            .ok_or_else(|| invalid("single-mode pressure port has zero static compliance"))?;
        if !pressure_fraction.is_finite() || pressure_fraction < admission.minimum_pressure_compliance_fraction {
            return Err(invalid("playing mode carries insufficient pressure-work compliance"));
        }
        let error = self.first_static_opening_relative_error()
            .ok_or_else(|| invalid("single-mode exact opening response vanishes"))?;
        if !error.is_finite() || error > tolerance {
            return Err(invalid("playing mode misrepresents the exact pressure-driven static opening"));
        }
        Ok(())
    }

    fn values(&self) -> [f64; 8] {
        [self.exact_mean_force_m_per_n, self.exact_volume_m3_per_pa,
         self.exact_opening_m_per_pa, self.retained_volume_m3_per_pa,
         self.retained_opening_m_per_pa, self.first_mean_force_m_per_n,
         self.first_volume_m3_per_pa, self.first_opening_m_per_pa]
    }
    pub(super) fn append_identity(&self, hash: &mut DomainHasher) {
        for value in self.values() { hash.update(&value.to_bits().to_le_bytes()); }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn modes() -> Vec<ReducedMode> {
        [1.0, 100.0].map(|k: f64| ReducedMode {
            frequency_hz: k.sqrt() / core::f64::consts::TAU,
            mass_kg: 1.0, stiffness_n_m: k, damping_n_s_m: 0.0,
            participation_sqrt_kg: 1.0,
        }).to_vec()
    }
    fn policy() -> SingleModeValveAdmission {
        SingleModeValveAdmission {
            minimum_opening_compliance_fraction: 0.98,
            minimum_pressure_compliance_fraction: 0.98,
            maximum_static_opening_relative_error: 0.02,
        }
    }
    #[test]
    fn complete_bank_does_not_certify_first_mode_for_a_different_load() {
        // K=diag(1,100), M=I, g=(1,1), l=(0.01,10).
        // Mode zero supplies 99% of opening-force compliance, but only
        // 0.01% of the pressure-work compliance. All modes retain 100%.
        let c = PressureCompliance::from_modes(&modes(), &[0.01, 10.0], 1.01, 1.0001, 0.11).unwrap();
        assert!(c.first_opening_compliance_fraction() > 0.99);
        assert!(c.first_pressure_compliance_fraction().unwrap() < 0.0001);
        assert!((c.retained_volume_m3_per_pa() - 1.0001).abs() < 1e-14);
        assert!(c.admit_single_mode(policy()).is_err());
    }
    #[test]
    fn parallel_load_passes_and_sign_reversal_does_not_become_a_closing_reed() {
        let c = PressureCompliance::from_modes(&modes(), &[1.0, 1.0], 1.01, 1.01, 1.01).unwrap();
        c.admit_single_mode(policy()).unwrap();
        let reversed = PressureCompliance::from_modes(&modes(), &[-1.0, -1.0], 1.01, 1.01, -1.01).unwrap();
        assert_eq!(c.retained_volume_m3_per_pa(), reversed.retained_volume_m3_per_pa());
        assert_eq!(c.retained_opening_m_per_pa(), -reversed.retained_opening_m_per_pa());
        assert!(reversed.admit_single_mode(policy()).is_err());
    }
    #[test]
    fn signed_modal_cancellation_requires_the_static_opening_gate() {
        let c = PressureCompliance::from_modes(&modes(), &[1.0, -10.0], 1.01, 2.0, 0.9).unwrap();
        let mut p = policy();
        p.minimum_pressure_compliance_fraction = 0.4;
        assert!(c.admit_single_mode(p).is_err()); // 11.1% opening error
        p.maximum_static_opening_relative_error = 0.12;
        c.admit_single_mode(p).unwrap();
    }
    #[test]
    fn zero_ports_nonfinite_data_and_invalid_policies_refuse() {
        let zero = PressureCompliance::from_modes(&modes(), &[0.0, 0.0], 1.01, 0.0, 0.0).unwrap();
        assert_eq!(zero.first_pressure_compliance_fraction(), None);
        assert_eq!(zero.first_static_opening_relative_error(), None);
        assert!(zero.admit_single_mode(policy()).is_err());
        let c = PressureCompliance::from_modes(&modes(), &[1.0, 1.0], 1.01, 1.01, 1.01).unwrap();
        for bad in [0.0, -1.0, 1.01, f64::NAN, f64::INFINITY] {
            let mut p = policy(); p.minimum_pressure_compliance_fraction = bad;
            assert!(c.admit_single_mode(p).is_err());
            p = policy(); p.minimum_opening_compliance_fraction = bad;
            assert!(c.admit_single_mode(p).is_err());
        }
        for bad in [-1.0, 1.0, f64::NAN, f64::INFINITY] {
            let mut p = policy(); p.maximum_static_opening_relative_error = bad;
            assert!(c.admit_single_mode(p).is_err());
        }
        assert!(PressureCompliance::from_modes(&modes(), &[1.0], 1.01, 1.0, 1.0).is_err());
        assert!(PressureCompliance::from_modes(&modes(), &[f64::NAN, 1.0], 1.01, 1.0, 1.0).is_err());
        assert!(PressureCompliance::from_modes(&modes(), &[1.0, 1.0], 1.01, -1.0, 1.0).is_err());
    }
}
