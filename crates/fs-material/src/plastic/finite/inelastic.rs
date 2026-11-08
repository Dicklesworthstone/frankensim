//! Complete finite material points for the original spatial Newton owner.
//!
//! A spatial solve must retain the actual constitutive state and kinetics, not
//! unwrap a creep law to its elastic/plastic inspection backbone. This protocol
//! preserves the complete candidate history and the original frozen interval.
//! Every stress, tangent, potential and physical heat remains the native result.

use super::combined::FiniteNortonJ2;
use super::creep::FiniteNortonCreep;
use super::creep::recovery::{FiniteRecoveringCreepState, FiniteRecoveringNortonCreep};
use super::rate::FiniteOverstressJ2;
use super::{FiniteJ2Error, FiniteJ2Material, FiniteJ2Response, FiniteJ2State};

/// An immutable finite material point with its complete, explicitly typed memory.
///
/// All energies refer to the original volume. The returned incremental potential
/// is stored energy plus the integrated flow potential; physical inelastic heat
/// is [`FiniteJ2Response::inelastic_dissipation_j_m3`] and can differ from that
/// flow potential. Consumers must preserve both accounts independently.
///
/// A spatial Newton iteration repeatedly evaluates the same accepted history.
/// It commits the complete returned `State` only after the surrounding step is
/// accepted. Neither line search nor accepted-state inspection consumes another
/// interval or recovers a second time.
pub trait FiniteInelasticMaterial {
    /// Full native history, including independent resistance when selected.
    type State: Clone + core::fmt::Debug;

    /// Virgin history in this law's native convention.
    fn initial_state(&self) -> Self::State;

    /// Frozen physical interval for rate-dependent flow, or `None` for a
    /// genuinely rate-independent law. A dynamic consumer must require this
    /// interval, when present, to equal its own physical step duration. A
    /// quasistatic consumer must not mix different kinetic clocks in one step.
    fn duration_s(&self) -> Option<f64>;

    /// Native implicit response and exact discrete tangent at frozen old history.
    fn response(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2Error>;

    /// The identical return and potential without assembling a tangent.
    fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<([f64; 9], f64), FiniteJ2Error>;

    /// Inspect energy at fixed complete history without advancing any mechanism.
    fn stored_energy(&self, f: &[f64; 9], state: &Self::State) -> Result<f64, FiniteJ2Error>;

    /// Inspect finite history and the accepted deformation. The nonnegative
    /// finite tolerance applies to equilibrium plastic yield when present;
    /// legitimate accepted creep stress or Perzyna overstress is retained.
    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error>;
}

/// Existing native and source-resolved rate-independent callers retain their
/// original owner, state, yield gate and exact arithmetic.
impl<T: FiniteJ2Material + ?Sized> FiniteInelasticMaterial for T {
    type State = FiniteJ2State;

    fn initial_state(&self) -> Self::State {
        self.finite_j2().initial_state()
    }

    fn duration_s(&self) -> Option<f64> {
        None
    }

    fn response(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2Error> {
        self.finite_j2().response(f, old)
    }

    fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<([f64; 9], f64), FiniteJ2Error> {
        self.finite_j2().stress_and_potential(f, old)
    }

    fn stored_energy(&self, f: &[f64; 9], state: &Self::State) -> Result<f64, FiniteJ2Error> {
        self.finite_j2().stored_energy(f, state)
    }

    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        self.finite_j2()
            .validate_committed(f, state, yield_tolerance_pa)
    }
}

impl FiniteInelasticMaterial for FiniteNortonCreep {
    type State = FiniteJ2State;

    fn initial_state(&self) -> Self::State {
        Self::initial_state(self)
    }

    fn duration_s(&self) -> Option<f64> {
        Some(Self::duration_s(self))
    }

    fn response(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2Error> {
        Self::response(self, f, old)
    }

    fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<([f64; 9], f64), FiniteJ2Error> {
        Self::stress_and_potential(self, f, old)
    }

    fn stored_energy(&self, f: &[f64; 9], state: &Self::State) -> Result<f64, FiniteJ2Error> {
        Self::stored_energy(self, f, state)
    }

    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        validate_tolerance(yield_tolerance_pa)?;
        Self::validate_committed(self, f, state)
    }
}

impl FiniteInelasticMaterial for FiniteRecoveringNortonCreep {
    type State = FiniteRecoveringCreepState;

    fn initial_state(&self) -> Self::State {
        Self::initial_state(self)
    }

    fn duration_s(&self) -> Option<f64> {
        Some(self.norton().duration_s())
    }

    fn response(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2Error> {
        Self::response(self, f, old)
    }

    fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<([f64; 9], f64), FiniteJ2Error> {
        Self::stress_and_potential(self, f, old)
    }

    fn stored_energy(&self, f: &[f64; 9], state: &Self::State) -> Result<f64, FiniteJ2Error> {
        Self::stored_energy(self, f, state)
    }

    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        validate_tolerance(yield_tolerance_pa)?;
        Self::validate_committed(self, f, state)
    }
}

// These existing mechanisms already expose exactly this accepted-state API.
// The macro forwards to their concrete inherent methods, retaining the
// Perzyna overstress gate and the combined law's equilibrium plastic gate.
macro_rules! impl_native_kinetic {
    ($law:ty) => {
        impl FiniteInelasticMaterial for $law {
            type State = FiniteJ2State;

            fn initial_state(&self) -> Self::State {
                <$law>::initial_state(self)
            }

            fn duration_s(&self) -> Option<f64> {
                Some(<$law>::duration_s(self))
            }

            fn response(
                &self,
                f: &[f64; 9],
                old: &Self::State,
            ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2Error> {
                <$law>::response(self, f, old)
            }

            fn stress_and_potential(
                &self,
                f: &[f64; 9],
                old: &Self::State,
            ) -> Result<([f64; 9], f64), FiniteJ2Error> {
                <$law>::stress_and_potential(self, f, old)
            }

            fn stored_energy(
                &self,
                f: &[f64; 9],
                state: &Self::State,
            ) -> Result<f64, FiniteJ2Error> {
                <$law>::stored_energy(self, f, state)
            }

            fn validate_committed(
                &self,
                f: &[f64; 9],
                state: &Self::State,
                yield_tolerance_pa: f64,
            ) -> Result<(), FiniteJ2Error> {
                <$law>::validate_committed(self, f, state, yield_tolerance_pa)
            }
        }
    };
}

impl_native_kinetic!(FiniteOverstressJ2);
impl_native_kinetic!(FiniteNortonJ2);

fn validate_tolerance(value: f64) -> Result<(), FiniteJ2Error> {
    if !value.is_finite() || value < 0.0 {
        Err(FiniteJ2Error::Parameters(
            "yield tolerance must be finite and nonnegative",
        ))
    } else {
        Ok(())
    }
}
