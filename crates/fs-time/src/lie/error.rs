//! Shared refusal type for the native SO(3) integration owners.

use fs_ga::GaError;

/// Typed failures for `SO(3)` integration steps.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LieStepError {
    /// A group/tangent value failed `fs-ga` validation.
    Group(GaError),
    /// Step size or principal inertia was non-finite, zero, or non-positive as
    /// required by the operation.
    InvalidParameter {
        /// Parameter family that was invalid.
        context: &'static str,
    },
}

impl core::fmt::Display for LieStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Group(error) => write!(f, "SO(3) group refusal: {error}"),
            Self::InvalidParameter { context } => write!(f, "invalid {context}"),
        }
    }
}

impl std::error::Error for LieStepError {}

impl From<GaError> for LieStepError {
    fn from(value: GaError) -> Self {
        Self::Group(value)
    }
}
