//! Physical wall-material sensitivity on the same immutable FE pole owner.

pub use super::super::super::frequency::coupled::poles::{
    MechanicalRelaxationParameterDirection, MechanicalResonanceParameterDerivative,
    MechanicalWallParameterDerivative, MechanicalWallParameterDirection,
};
use super::{MechanicalPoleError, MechanicalPoleLimits, PreparedMechanicalResonance};

impl PreparedMechanicalResonance<'_, '_, '_, '_> {
    /// Differentiate an isolated simple free pole and its full original mode
    /// along explicitly supplied wall-material coefficient directions.
    ///
    /// Directions are derivatives per caller parameter, not replacement values.
    /// Every selected patch supplies its original surface mass, stiffness,
    /// resistance and complete Maxwell spectrum directions. Unselected patches,
    /// bulk medium, geometry, boundary traces and pressure-release constraints
    /// remain fixed. No temperature law, validity extension or finite parameter
    /// step is inferred. Signed directions do not alter the admitted base state.
    ///
    /// The existing full bordered factor is reused at O(N^2) solve cost, with
    /// O(N+nnz) physical assembly/recovery work and O(N) additional payload.
    /// Numerical equilibration and the reported physical normalization pivot
    /// remain fixed. Original memory derivatives include the square-root
    /// stiffness coordinate's chain rule, including dark material modes.
    ///
    /// `limits` bounds the extra physical payload and one independently checked
    /// differentiated-row audit; shared preparation bounds the solve. Any source,
    /// repeated/defective-pole, residual or cancellation refusal publishes no
    /// partial result. This is a finite-model free-pole derivative, not a driven
    /// response-peak, line-width or certified material-uncertainty bound.
    pub fn wall_parameter_derivative(
        &self,
        directions: &[MechanicalWallParameterDirection<'_>],
        limits: MechanicalPoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalResonanceParameterDerivative, MechanicalPoleError> {
        self.response
            .wall_parameter_derivative(directions, limits, cancelled)
    }
}
