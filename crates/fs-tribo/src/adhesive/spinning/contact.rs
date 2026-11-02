//! Simultaneous adhesive sliding and twisting against the actual body response.
//!
//! The static set is the resultant of all local tractions with magnitude at
//! most the sourced static shear strength. Its boundary, and the kinetic
//! response, use the same uniform-disk integral as [`super::AdhesiveSlidingPatch`].
//! Neither a friction ellipse nor independent force/torque saturation is used.
//! The existing bounded circular-support projection retains all nine entries
//! of the mechanical compliance and all of its work and cancellation gates.

use super::{AdhesiveSlidingPatch, AdhesiveSlidingPatchParameters, shape_factors};
use crate::TriboError;
use crate::sliding_patch::contact::{CircularWrenchLaw, advance_circular};

/// Numerical controls shared by the circular force/moment projection.
pub use crate::sliding_patch::contact::HertzWrenchConfig as AdhesiveWrenchConfig;
/// Refusals from the common mechanical projection or original traction law.
pub use crate::sliding_patch::contact::HertzWrenchError as AdhesiveWrenchError;
/// Released translation, rotation and complete radius-scaled compliance.
pub use crate::sliding_patch::contact::HertzWrenchInput as AdhesiveWrenchInput;
/// Static, sliding or zero-strength branch of the resolved adhesive patch.
pub use crate::sliding_patch::contact::HertzWrenchMode as AdhesiveWrenchMode;
/// Physically supplied stationary or moving phase; no stopping event is inferred.
pub use crate::sliding_patch::contact::HertzWrenchPhase as AdhesiveWrenchPhase;
/// Candidate coupled force/moment, motion, local dissipation and work residuals.
pub use crate::sliding_patch::contact::HertzWrenchStep as AdhesiveWrenchStep;

/// Actual attached radius and independently supplied uniform junction strengths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveWrenchParameters {
    /// Actual solid contact radius [m], strictly positive for this attached solve.
    pub contact_radius_m: f64,
    /// Uniform local static shear bound [Pa], at least the kinetic strength.
    pub static_shear_strength_pa: f64,
    /// Uniform local gross-sliding shear strength [Pa], nonnegative.
    pub kinetic_shear_strength_pa: f64,
}

/// Frozen attached junction; normal force is not an input to its shear capacity.
/// It can therefore transmit both shear and torque at zero or tensile net load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveWrenchContact {
    parameters: AdhesiveWrenchParameters,
    static_patch: AdhesiveSlidingPatch,
    kinetic_patch: AdhesiveSlidingPatch,
}

impl AdhesiveWrenchContact {
    /// Bind real geometry and strengths without assigning material authority.
    /// A zero-area detached state belongs to the normal owner; it has no patch
    /// wrench and must not be given an invented radius for this scaled solve.
    pub fn new(parameters: AdhesiveWrenchParameters) -> Result<Self, AdhesiveWrenchError> {
        if !parameters.contact_radius_m.is_finite() || parameters.contact_radius_m <= 0.0 {
            return Err(AdhesiveWrenchError::Invalid(
                "positive attached adhesive radius",
            ));
        }
        let patch = |shear_strength_pa| {
            AdhesiveSlidingPatch::new(AdhesiveSlidingPatchParameters {
                contact_radius_m: parameters.contact_radius_m,
                shear_strength_pa,
            })
            .map_err(AdhesiveWrenchError::Law)
        };
        let static_patch = patch(parameters.static_shear_strength_pa)?;
        let kinetic_patch = patch(parameters.kinetic_shear_strength_pa)?;
        if parameters.static_shear_strength_pa < parameters.kinetic_shear_strength_pa {
            return Err(AdhesiveWrenchError::Invalid(
                "static adhesive strength below kinetic strength",
            ));
        }
        Ok(Self {
            parameters,
            static_patch,
            kinetic_patch,
        })
    }

    /// Original real contact radius and local shear strengths.
    #[must_use]
    pub const fn parameters(&self) -> AdhesiveWrenchParameters {
        self.parameters
    }

    /// Original kinetic integral for an independent actual-motion work check.
    #[must_use]
    pub const fn kinetic_patch(&self) -> &AdhesiveSlidingPatch {
        &self.kinetic_patch
    }

    /// Resolve a static-admissible held-effort sample against the full response.
    /// Normal evolution, actual stopping events, shear elasticity and any heat
    /// destination remain with their existing physical owners.
    pub fn advance(
        &self,
        input: &AdhesiveWrenchInput,
        config: AdhesiveWrenchConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdhesiveWrenchStep, AdhesiveWrenchError> {
        advance_circular(self, input, config, cancelled)
    }

    /// A moving junction cannot be captured by the larger static strength.
    /// This preserves the physical distinction between a zero mean increment
    /// and an endpoint stop, which the receiving body must resolve separately.
    pub fn advance_phase(
        &self,
        input: &AdhesiveWrenchInput,
        phase: AdhesiveWrenchPhase,
        config: AdhesiveWrenchConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdhesiveWrenchStep, AdhesiveWrenchError> {
        let law = match phase {
            AdhesiveWrenchPhase::StaticAdmissible => *self,
            AdhesiveWrenchPhase::Kinetic => Self {
                static_patch: self.kinetic_patch,
                ..*self
            },
        };
        advance_circular(&law, input, config, cancelled)
    }
}

impl CircularWrenchLaw for AdhesiveWrenchContact {
    fn radius(&self) -> f64 {
        self.parameters.contact_radius_m
    }
    fn static_capacity(&self) -> f64 {
        self.static_patch.sliding_capacity_n()
    }
    fn kinetic_capacity(&self) -> f64 {
        self.kinetic_patch.sliding_capacity_n()
    }
    fn spin_radius_fraction(&self) -> f64 {
        2.0 / 3.0
    }
    fn shape(&self, speed: f64, edge_speed: f64) -> Result<(f64, f64), TriboError> {
        shape_factors(speed, edge_speed)
    }
}
