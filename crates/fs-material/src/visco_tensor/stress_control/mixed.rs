//! Componentwise stress/strain control on the same original tensor history.

use super::{
    PreparedTensorStressControl, Result, TensorPronyHistory, TensorStressControlBudget,
    TensorStressRampStep, Voigt,
};

/// Exact mixed loading of a full three-dimensional generalized Maxwell law.
///
/// Each `true` entry in `stress_components` prescribes stress in that component;
/// each `false` entry prescribes strain. Order is [xx, yy, zz, xy, yz, zx], with
/// tensor shear at the public boundary. The mask belongs to the prepared
/// controller, while the complete memory belongs to the original history.
/// Controllers with different masks may therefore alternate on that history
/// without translating, truncating or reinitializing its Maxwell branches.
///
/// For a plate in the xy plane, use `[false, false, true, false, true, true]`:
/// prescribe xx, yy and xy strain and set zz, yz and zx stress to zero. This
/// solves the original hereditary plane-stress constraints throughout each
/// ramp. An elastic Schur complement applied separately to each branch would
/// not, in general, give the same history.
///
/// Preparation shares the full-stress positive compliance factorization;
/// ramps use exact scalar relaxation responses and independent positive heat.
/// Both prescribed and reaction ports contribute to mechanical work. The
/// original whole-path strain bound, resource ceilings and atomic cancellation
/// apply. Nonzero Kelvin viscosity is not admitted by this controller.
#[derive(Clone, Debug)]
pub struct PreparedTensorMixedControl {
    inner: PreparedTensorStressControl,
}

impl PreparedTensorMixedControl {
    /// Bind a control partition to the immutable model of an existing history.
    /// All-stress and all-strain masks are permitted.
    ///
    /// # Errors
    /// Nonzero Kelvin viscosity, exhausted preparation budgets, unresolved
    /// numerical factors, unrepresentable source scales, or cancellation.
    pub fn new(
        history: &TensorPronyHistory,
        stress_components: [bool; 6],
        budget: TensorStressControlBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Ok(Self {
            inner: PreparedTensorStressControl::prepare(
                history,
                stress_components,
                budget,
                cancelled,
            )?,
        })
    }

    /// The fixed partition in tensor-shear Voigt order; true means stress.
    #[must_use]
    pub const fn stress_components(&self) -> [bool; 6] {
        self.inner.stress_components
    }

    /// Number of original retained memory coordinates, including null rows.
    #[must_use]
    pub fn memory_coordinates(&self) -> usize {
        self.inner.memory_coordinates()
    }

    /// Advance a linear ramp from the history's current control values.
    ///
    /// `prescribed[i]` is an endpoint stress [Pa] for a stress-controlled
    /// component, and an endpoint strain for a strain-controlled component.
    /// Shears are tensor shears in both cases, never engineering shears.
    /// Zero duration is the corresponding coupled elastic jump. The report
    /// contains all six physical stresses and strains, including reactions.
    ///
    /// # Errors
    /// Invalid data, a different material binding, unresolved source or energy
    /// balance, whole-path strain limit exhaustion, or cancellation. The
    /// history remains unchanged on every refusal.
    pub fn advance(
        &self,
        history: &mut TensorPronyHistory,
        prescribed: Voigt,
        duration_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<TensorStressRampStep> {
        self.inner
            .advance_controlled(history, prescribed, duration_s, cancelled)
    }
}
