//! Irreversible pressure-weighted sliding on an original admitted source path.
//!
//! The changing-load constitutive path runs once. Its local sliding measure
//! retains loading, unloading and separation, including sliding caused by
//! release of retained shear memory with no imposed tangential displacement.

use super::{SourceMindlinPathError, SourceMindlinPathInput, SourceMindlinPathStep};
use crate::interface_binding::mindlin::SourceMindlinError;
use crate::interface_binding::mindlin::loading::{
    SourceHertzMindlinContact, SourceHertzMindlinState,
};
use crate::normal_patch::NormalPatchError;

/// Read-only original-source path and its pressure-weighted sliding integral.
/// Consume this receipt to transfer the candidate into the accepted history.
#[derive(Debug)]
pub struct SourceMindlinPathSlidingAdvance {
    path: SourceMindlinPathStep,
    pressure_weighted_sliding_j: f64,
}

impl SourceMindlinPathSlidingAdvance {
    /// The original source-bound history, solved motion and separate work ports.
    pub const fn path_step(&self) -> &SourceMindlinPathStep {
        &self.path
    }

    /// Integral of pressure times unsigned local irreversible sliding [N m].
    /// This is independent of normal elastic storage and is not wear volume.
    pub const fn pressure_weighted_sliding_j(&self) -> f64 {
        self.pressure_weighted_sliding_j
    }

    /// Transfer the complete candidate without recomputing its constitutive path.
    pub fn into_path_step(self) -> SourceMindlinPathStep {
        self.path
    }
}

impl SourceHertzMindlinContact {
    /// Resolve the original simultaneous normal/shear path and its local sliding.
    ///
    /// For positive constant friction the native pressure/sliding integral is
    /// dissipated shear work divided by that original coefficient. At zero
    /// friction the entire footprint slides freely: the original normal source
    /// integrates the nonlinear Hertz force along affine indentation, then
    /// multiplies it by the actual solved tangential travel. Neither an endpoint
    /// load average nor offered body motion replaces this path measure.
    ///
    /// Admission, state ownership, work and cancellation are shared with
    /// [`Self::advance_coupled_path`]. No accepted state is mutated.
    pub fn advance_coupled_path_with_sliding(
        &self,
        state: &SourceHertzMindlinState,
        input: &SourceMindlinPathInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceMindlinPathSlidingAdvance, SourceMindlinPathError> {
        let (response, native_input) = self.admit_coupled_path(state, input, cancelled)?;
        let native = self.law.advance_coupled_path_with_sliding(
            &state.constitutive,
            &native_input,
            cancelled,
        )?;
        let measure = if self.law.parameters().friction_coefficient == 0.0 {
            let mean = self
                .source
                .normal_source()
                .elastic_path_average(state.indentation_m, input.indentation_m, input.dt_s)
                .map_err(SourceMindlinError::from)?;
            mean.mean_force_n * native.path_step().tangential.slip_increment_m[0].abs()
        } else {
            native.pressure_weighted_sliding_j()
        };
        if !measure.is_finite() || measure < 0.0 {
            return Err(SourceMindlinError::from(NormalPatchError::Overflow {
                field: "source pressure-weighted irreversible sliding",
            })
            .into());
        }
        let path =
            self.bind_coupled_path(state, input, response, native.into_path_step(), cancelled)?;
        Ok(SourceMindlinPathSlidingAdvance {
            path,
            pressure_weighted_sliding_j: measure,
        })
    }
}
