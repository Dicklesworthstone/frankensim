//! Conservative, finite-lifetime relaxation of an authored excitation store.
//!
//! A fixed lifetime and radiative branching yield define `dU/dt = P - U/tau`.
//! The constant-pump solution and both integrated decay channels are evaluated
//! analytically, so a short optical lifetime does not impose thermal substeps.
//! All joules refer to the caller's emitting-state energy convention. Radiative
//! decay releases that energy BEFORE any photon-spectrum or Stokes partition;
//! it is not automatically escaped optical energy. Pump vibrational relaxation,
//! photon energies, reabsorption and heat transport belong to the caller.
//!
//! Quantum yield counts photons, not joules (IUPAC Gold Book,
//! <https://doi.org/10.1351/goldbook.Q04991>). For one fixed-energy emitting
//! state the radiative branching probability is `kr/(kr+knr)` and the lifetime
//! is `1/(kr+knr)`. Adding a nonradiative channel changes BOTH consistently.
//! No material constants, saturation or stimulated emission are inferred.
//! [`transfer`] adds explicitly authored singlet/triplet crossing and reverse
//! crossing with separate thermal exchanges; [`recycling`] integrates repeated
//! equal-energy photon absorption between stores. Coefficients must remain fixed
//! over each step; a thermal owner must resolve changing-temperature rates.
//!
//! Scalar work is bounded and allocation-free. Elementary functions use the
//! existing deterministic math owner; a refusal returns no candidate state.

use core::fmt;
use fs_math::{decay, eft::two_sum};

/// Positive singlet/triplet kinetics, delayed emission and transfer heat.
#[path = "luminescence_transfer.rs"]
pub mod transfer;

/// Transient equal-energy photon recycling and explicit escape/heat accounts.
#[path = "luminescence_recycling.rs"]
pub mod recycling;

/// Invalid physical inputs or an unrepresentable candidate leave state untouched.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LuminescenceError {
    /// The caller supplied a nonfinite or out-of-domain physical quantity.
    Invalid(&'static str),
    /// A required finite physical result cannot be represented by binary64.
    Unrepresentable(&'static str),
}

impl fmt::Display for LuminescenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "luminescence: {self:?}")
    }
}
impl std::error::Error for LuminescenceError {}

/// One linear emitting-state reservoir with constant rates over a step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExcitationDecay {
    lifetime_s: f64,
    quantum_yield: f64,
}

/// Integrated channels on one accepted interval, all in joules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExcitationStep {
    /// Excitation energy remaining at the end of the interval.
    pub stored_energy_j: f64,
    /// Energy supplied to this store during the interval.
    pub pumped_energy_j: f64,
    /// Emitting-state energy released through radiative transitions, before
    /// photon-spectrum conversion; this is not necessarily escaped radiation.
    pub radiative_decay_energy_j: f64,
    /// Emitting-state energy released as local nonradiative heat.
    pub nonradiative_heat_j: f64,
    /// Represented `old_store + pumped - new_store - radiative - heat`.
    /// This signed arithmetic residual is reported, never hidden in a channel.
    pub residual_j: f64,
}

impl ExcitationDecay {
    /// Admit a positive finite lifetime and a branching yield in `[0, 1]`.
    pub fn try_new(lifetime_s: f64, quantum_yield: f64) -> Result<Self, LuminescenceError> {
        if !lifetime_s.is_finite()
            || lifetime_s <= 0.0
            || !quantum_yield.is_finite()
            || !(0.0..=1.0).contains(&quantum_yield)
        {
            return Err(LuminescenceError::Invalid("lifetime or quantum yield"));
        }
        Ok(Self {
            lifetime_s,
            quantum_yield,
        })
    }

    /// Total excited-state lifetime, seconds.
    pub const fn lifetime_s(&self) -> f64 {
        self.lifetime_s
    }
    /// Radiative branching probability per relaxed excitation.
    pub const fn quantum_yield(&self) -> f64 {
        self.quantum_yield
    }

    /// Add a supplied nonradiative decay rate in inverse seconds. The original
    /// radiative rate is preserved; total lifetime and yield both decrease.
    pub fn with_nonradiative_rate(self, rate_s: f64) -> Result<Self, LuminescenceError> {
        if !nonnegative(rate_s) {
            return Err(LuminescenceError::Invalid("nonradiative rate"));
        }
        if rate_s == 0.0 {
            return Ok(self);
        }
        let inverse_rate = 1.0 / rate_s;
        let (lifetime, fraction) = if self.lifetime_s <= inverse_rate {
            let fraction = 1.0 / (1.0 + rate_s * self.lifetime_s);
            (self.lifetime_s * fraction, fraction)
        } else {
            let ratio = inverse_rate / self.lifetime_s;
            (inverse_rate / (1.0 + ratio), ratio / (1.0 + ratio))
        };
        let quantum_yield = self.quantum_yield * fraction;
        if lifetime == 0.0 || (self.quantum_yield > 0.0 && quantum_yield == 0.0) {
            return Err(LuminescenceError::Unrepresentable(
                "quenched lifetime or yield",
            ));
        }
        Self::try_new(lifetime, quantum_yield)
    }

    /// Advance constant excitation power analytically, including zero-duration
    /// and pump-off intervals. Inputs are finite nonnegative joules, watts and
    /// seconds. The returned channels integrate the WHOLE interval, not just
    /// the endpoint emission rate. The caller owns acceptance and cancellation.
    pub fn advance(
        &self,
        stored_energy_j: f64,
        pump_power_w: f64,
        duration_s: f64,
    ) -> Result<ExcitationStep, LuminescenceError> {
        if ![stored_energy_j, pump_power_w, duration_s]
            .into_iter()
            .all(nonnegative)
        {
            return Err(LuminescenceError::Invalid("store, pump or duration"));
        }
        let pumped = finite(pump_power_w * duration_s)?;
        let exposure = duration_s / self.lifetime_s;
        let lost = decay::decayed_fraction(exposure);
        // An underflowing h/tau can still release representable energy from a
        // large store. Rescale before division in this rare linear-limit case.
        let old_release = if exposure < f64::MIN_POSITIVE {
            small_release(stored_energy_j, duration_s, self.lifetime_s)
        } else {
            stored_energy_j * lost
        };
        let (pump_store, pump_release) = if exposure >= 1.0 {
            // P*tau is bounded by P*h here, even when h/tau overflows.
            let remaining = (pump_power_w * self.lifetime_s) * lost;
            (remaining, pumped - remaining)
        } else {
            let released = if exposure < f64::MIN_POSITIVE {
                0.5 * small_release(pumped, duration_s, self.lifetime_s)
            } else {
                pumped * pump_release_fraction(exposure, lost)
            };
            (pumped - released, released)
        };
        let stored = finite(decay::attenuate(stored_energy_j, exposure) + pump_store)?;
        // Split before summing, avoiding an unnecessary overflow of the total
        // released energy when its two individually reportable channels fit.
        let radiative =
            finite(old_release * self.quantum_yield + pump_release * self.quantum_yield)?;
        let dark = 1.0 - self.quantum_yield;
        let heat = finite(old_release * dark + pump_release * dark)?;
        let (storage_change, e1) = two_sum(stored_energy_j, -stored);
        let (net_pump, e2) = two_sum(pumped, -radiative);
        let (balance, e3) = two_sum(storage_change, net_pump);
        let (balance, e4) = two_sum(balance, -heat);
        let residual = balance + ((e1 + e2) + (e3 + e4));
        Ok(ExcitationStep {
            stored_energy_j: stored,
            pumped_energy_j: pumped,
            radiative_decay_energy_j: radiative,
            nonradiative_heat_j: heat,
            residual_j: finite(residual)?,
        })
    }
}

/// An authored thermally activated NONRADIATIVE channel `A exp(-Ea/(k_B T))`.
/// The physical source owner supplies and enforces its valid temperature domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArrheniusQuenching {
    prefactor_rate_s: f64,
    activation_energy_j: f64,
}

impl ArrheniusQuenching {
    /// Admit finite nonnegative prefactor [1/s] and activation energy [J].
    pub fn try_new(
        prefactor_rate_s: f64,
        activation_energy_j: f64,
    ) -> Result<Self, LuminescenceError> {
        if !nonnegative(prefactor_rate_s) || !nonnegative(activation_energy_j) {
            return Err(LuminescenceError::Invalid("Arrhenius prefactor or barrier"));
        }
        Ok(Self {
            prefactor_rate_s,
            activation_energy_j,
        })
    }
    /// Supplied attempt-rate prefactor, inverse seconds.
    pub const fn prefactor_rate_s(&self) -> f64 {
        self.prefactor_rate_s
    }
    /// Supplied activation energy per excitation, joules.
    pub const fn activation_energy_j(&self) -> f64 {
        self.activation_energy_j
    }

    /// Evaluate the supplied law at finite positive absolute temperature [K].
    pub fn rate_at_temperature(&self, temperature_k: f64) -> Result<f64, LuminescenceError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(LuminescenceError::Invalid("absolute temperature"));
        }
        const BOLTZMANN_J_K: f64 = 1.380_649e-23;
        // Divide first below 1 K to avoid underflowing k_B*T prematurely.
        let exposure = if temperature_k < 1.0 {
            (self.activation_energy_j / temperature_k) / BOLTZMANN_J_K
        } else {
            self.activation_energy_j / (BOLTZMANN_J_K * temperature_k)
        };
        finite(decay::attenuate(self.prefactor_rate_s, exposure))
    }
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64) -> Result<f64, LuminescenceError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(LuminescenceError::Unrepresentable("energy or rate"))
    }
}

// 1 - (1-exp(-z))/z. The ninth-order Taylor polynomial avoids cancellation
// of the O(z) result; at z<0.1 its omitted term is <2.6e-18 absolute.
fn pump_release_fraction(z: f64, lost: f64) -> f64 {
    if z < 0.1 {
        z * (0.5
            + z * (-1.0 / 6.0
                + z * (1.0 / 24.0
                    + z * (-1.0 / 120.0
                        + z * (1.0 / 720.0
                            + z * (-1.0 / 5040.0
                                + z * (1.0 / 40320.0 + z * (-1.0 / 362880.0 + z / 3628800.0))))))))
    } else {
        1.0 - lost / z
    }
}

// Used only when h/tau is subnormal: its nonlinear correction cannot affect
// a finite binary64 result. Select an ordering that retains a normal scale.
fn small_release(energy: f64, h: f64, tau: f64) -> f64 {
    if energy == 0.0 || h == 0.0 {
        return 0.0;
    }
    let rate = energy / tau;
    if rate >= f64::MIN_POSITIVE {
        rate * h
    } else {
        (energy * h) / tau
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
            "actual={actual:e}, expected={expected:e}"
        );
    }

    #[test]
    fn g1_constant_pump_matches_independent_closed_form_and_balance() {
        let model = ExcitationDecay::try_new(2.0, 0.75).unwrap();
        let step = model.advance(3.0, 4.0, 1.0).unwrap();
        let store = 8.0 - 5.0 * (-0.5_f64).exp();
        let released = 7.0 - store;
        close(step.stored_energy_j, store, 2e-15);
        close(step.pumped_energy_j, 4.0, 1e-15);
        close(step.radiative_decay_energy_j, 0.75 * released, 4e-15);
        close(step.nonradiative_heat_j, 0.25 * released, 4e-15);
        assert!(step.residual_j.abs() < 2e-15);
    }

    #[test]
    fn g0_tiny_step_retains_pump_release_without_subtracting_large_stores() {
        let model = ExcitationDecay::try_new(1.0, 1.0).unwrap();
        let step = model.advance(0.0, 1e20, 1e-20).unwrap();
        close(step.radiative_decay_energy_j, 5e-21, 2e-15);
        assert_eq!(step.stored_energy_j, 1.0);
        // The separately reported residual includes the unrepresentable store loss.
        close(step.residual_j, -5e-21, 2e-15);
        let long = ExcitationDecay::try_new(1e308, 1.0).unwrap();
        let step = long.advance(1e308, 1e20, 1e-20).unwrap();
        close(step.radiative_decay_energy_j, 1e-20, 2e-15);
        assert_eq!(step.stored_energy_j, 1e308);
    }

    #[test]
    fn g0_prompt_limit_and_overflowing_exposure_keep_the_remaining_store() {
        let model = ExcitationDecay::try_new(1e-320, 0.8).unwrap();
        let step = model.advance(1.0, 1e300, 1.0).unwrap();
        close(step.stored_energy_j, 1e300 * 1e-320, 2e-15);
        close(step.radiative_decay_energy_j, 8e299, 2e-15);
        close(step.nonradiative_heat_j, 2e299, 2e-15);
        assert!(step.residual_j.abs() < 3e284);
    }

    #[test]
    fn g1_pump_off_releases_stored_energy_and_retains_late_tails() {
        let model = ExcitationDecay::try_new(0.5, 0.4).unwrap();
        let charged = model.advance(0.0, 2.0, 3.0).unwrap();
        let dark = model.advance(charged.stored_energy_j, 0.0, 0.5).unwrap();
        close(
            dark.stored_energy_j,
            charged.stored_energy_j * (-1.0_f64).exp(),
            2e-15,
        );
        close(
            dark.radiative_decay_energy_j,
            0.4 * charged.stored_energy_j * (1.0 - (-1.0_f64).exp()),
            2e-15,
        );
        let tail = model.advance(1e300, 0.0, 400.0).unwrap();
        close(tail.stored_energy_j, (1e300_f64.ln() - 800.0).exp(), 2e-13);
    }

    #[test]
    fn g3_constant_rate_step_splitting_preserves_all_integrated_channels() {
        for lifetime in [1e-9, 0.2, 3.0, 1e6] {
            let model = ExcitationDecay::try_new(lifetime, 0.63).unwrap();
            let whole = model.advance(4.0, 7.0, 0.8).unwrap();
            let first = model.advance(4.0, 7.0, 0.3).unwrap();
            let second = model.advance(first.stored_energy_j, 7.0, 0.5).unwrap();
            close(second.stored_energy_j, whole.stored_energy_j, 2e-14);
            close(
                first.radiative_decay_energy_j + second.radiative_decay_energy_j,
                whole.radiative_decay_energy_j,
                2e-14,
            );
            close(
                first.nonradiative_heat_j + second.nonradiative_heat_j,
                whole.nonradiative_heat_j,
                2e-14,
            );
        }
    }

    #[test]
    fn g0_quenching_changes_lifetime_and_yield_with_one_radiative_rate() {
        let original = ExcitationDecay::try_new(2.0, 0.8).unwrap();
        let quenched = original.with_nonradiative_rate(1.5).unwrap();
        assert_eq!(quenched.lifetime_s(), 0.5);
        assert_eq!(quenched.quantum_yield(), 0.2);
        close(
            quenched.quantum_yield() / quenched.lifetime_s(),
            original.quantum_yield() / original.lifetime_s(),
            1e-15,
        );
        let strong = ExcitationDecay::try_new(1e200, 0.8)
            .unwrap()
            .with_nonradiative_rate(1e100)
            .unwrap();
        close(strong.lifetime_s(), 1e-100, 2e-15);
        close(strong.quantum_yield(), 8e-301, 2e-15);
    }

    #[test]
    fn g1_arrhenius_uses_kelvin_joules_and_retains_scaled_rates() {
        let barrier = 1.380_649e-23 * 300.0;
        let law = ArrheniusQuenching::try_new(1e9, barrier).unwrap();
        close(
            law.rate_at_temperature(300.0).unwrap(),
            1e9 / 1.0_f64.exp(),
            2e-15,
        );
        assert!(law.rate_at_temperature(600.0).unwrap() > law.rate_at_temperature(300.0).unwrap());
        let large = ArrheniusQuenching::try_new(1e300, 800.0 * barrier).unwrap();
        close(
            large.rate_at_temperature(300.0).unwrap(),
            (1e300_f64.ln() - 800.0).exp(),
            3e-13,
        );
        assert_eq!(
            ArrheniusQuenching::try_new(7.0, 0.0)
                .unwrap()
                .rate_at_temperature(f64::from_bits(1))
                .unwrap(),
            7.0
        );
    }

    #[test]
    fn g0_zero_time_yield_endpoints_and_invalid_inputs() {
        for yield_ in [0.0, 1.0] {
            let model = ExcitationDecay::try_new(1.0, yield_).unwrap();
            let zero = model.advance(3.0, f64::MAX, 0.0).unwrap();
            assert_eq!(zero.stored_energy_j, 3.0);
            assert_eq!(zero.pumped_energy_j, 0.0);
            assert_eq!(zero.residual_j, 0.0);
            let step = model.advance(1.0, 0.0, 1.0).unwrap();
            assert_eq!(
                if yield_ == 0.0 {
                    step.radiative_decay_energy_j
                } else {
                    step.nonradiative_heat_j
                },
                0.0
            );
            for invalid in [-1.0, f64::INFINITY, f64::NAN] {
                assert!(model.advance(invalid, 0.0, 1.0).is_err());
                assert!(model.advance(0.0, invalid, 1.0).is_err());
                assert!(model.advance(0.0, 0.0, invalid).is_err());
                assert!(model.with_nonradiative_rate(invalid).is_err());
                assert!(ArrheniusQuenching::try_new(invalid, 0.0).is_err());
            }
            assert!(model.advance(0.0, f64::MAX, 2.0).is_err());
        }
        for lifetime in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(ExcitationDecay::try_new(lifetime, 0.5).is_err());
        }
        for yield_ in [-0.1, 1.1, f64::NAN] {
            assert!(ExcitationDecay::try_new(1.0, yield_).is_err());
        }
    }
}
