//! Directional layer strength on the original thermal shell strain polynomial.
//!
//! The material owner recovers signed plane stresses on both faces of every
//! layer, including thermal self-stress, and rotates them into each layer's
//! elastic material axes. This adapter supplies the SAME native membrane and
//! director-bending field, including the condensed thermal bubble. No stress
//! is reconstructed from averaged section resultants.
//!
//! The six quadratic Bernstein controls enclose the represented strain field.
//! Layer stress is affine in strain and thickness. Maximum-stress utilization
//! is convex; the admitted Tsai-Wu radial utilization is a norm plus a linear
//! form and is also convex. Evaluating each control and both layer faces
//! therefore bounds the entire triangle and layer. Edge controls are not
//! physical samples: a bound reaching one does not establish actual failure.
//! [`PreparedThermalShell::laminate_strength_at`] observes a physical point.
//!
//! These are fixed-state, small-strain, in-plane elastic failure indications.
//! They do not evolve damage, remove stiffness, recover interlaminar failure,
//! or bound time intervals, model error or floating-point roundoff. The
//! material source owns all strength and interaction selections.

use super::{LaminateThermalField, PreparedThermalShell};
use crate::PlateError;
use crate::shell::mitc::strain::RealElementStrain;
use fs_material::laminate::strength::tsai_wu::{LaminateTsaiWuAssessment, PreparedLaminateTsaiWu};
use fs_material::laminate::strength::{
    LaminaFace, LaminateStrengthAssessment, LaminateStrengthError, PreparedLaminateStrength,
};
use fs_material::laminate::thermal::source::ResolvedLaminateThermal;

/// An explicit native criterion, already bound to its complete material state.
#[derive(Clone, Copy, Debug)]
pub enum LaminateStrengthCriterion<'a> {
    /// Independent tensile, compressive and in-plane shear limits.
    MaximumStress(&'a PreparedLaminateStrength),
    /// Coupled plane-stress envelope with explicitly sourced interaction.
    TsaiWu(&'a PreparedLaminateTsaiWu),
}

impl<'a> LaminateStrengthCriterion<'a> {
    fn thermal(self) -> &'a ResolvedLaminateThermal {
        match self {
            Self::MaximumStress(source) => source.thermal(),
            Self::TsaiWu(source) => source.thermal(),
        }
    }

    fn assess(
        self,
        field: &LaminateThermalField,
        membrane: [f64; 3],
        gradient: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ShellLayerStrength, ShellStrengthError> {
        let membrane = field.to_section(membrane);
        let gradient = field.to_section(gradient);
        Ok(match self {
            Self::MaximumStress(source) => {
                ShellLayerStrength::MaximumStress(source.assess(membrane, gradient, cancelled)?)
            }
            Self::TsaiWu(source) => {
                ShellLayerStrength::TsaiWu(source.assess(membrane, gradient, cancelled)?)
            }
        })
    }
}

/// Governing native layer face, retaining the selected criterion's full result.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ShellLayerStrength {
    /// All directional mode ratios and signed material-axis stress remain here.
    MaximumStress(LaminateStrengthAssessment),
    /// Coupled radial utilization and signed material-axis stress remain here.
    TsaiWu(LaminateTsaiWuAssessment),
}

impl ShellLayerStrength {
    pub fn utilization(&self) -> f64 {
        match self {
            Self::MaximumStress(value) => value.assessment.utilization(),
            Self::TsaiWu(value) => value.assessment.utilization(),
        }
    }

    pub const fn layer(&self) -> usize {
        match self {
            Self::MaximumStress(value) => value.layer,
            Self::TsaiWu(value) => value.layer,
        }
    }

    pub const fn face(&self) -> LaminaFace {
        match self {
            Self::MaximumStress(value) => value.face,
            Self::TsaiWu(value) => value.face,
        }
    }
}

/// Conservative spatial envelope at one actual nodal snapshot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateStrengthBound {
    control: usize,
    governing: ShellLayerStrength,
}

impl LaminateStrengthBound {
    /// Upper bound over the complete represented triangle and all layer depths.
    pub fn utilization_upper_bound(&self) -> f64 {
        self.governing.utilization()
    }

    /// Control 0/1/2 is a triangle vertex; 3/4/5 is the Bernstein edge control
    /// for 01/12/20. The latter is not the physical edge midpoint.
    pub const fn governing_control(&self) -> usize {
        self.control
    }

    /// Native criterion result at the governing control and layer face.
    /// Its stress is a control coefficient, not generally a physical sample.
    pub const fn governing(&self) -> &ShellLayerStrength {
        &self.governing
    }
}

/// Complete assessment refusal; no mechanics or source state is modified.
#[derive(Debug)]
pub enum ShellStrengthError {
    /// Invalid native element, displacement, strain, or cancellation.
    Native(PlateError),
    /// Original source strength or constitutive assessment refused.
    Material(LaminateStrengthError),
    /// The shell was not prepared from this exact shared material snapshot.
    SourceMismatch,
    /// Two face evaluations per layer and observation point must fit the budget.
    FaceBudget { required: usize, maximum: usize },
    /// The requested work count cannot be represented.
    BudgetOverflow,
}

impl core::fmt::Display for ShellStrengthError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-PLATE-LAMINATE-STRENGTH: {self:?}")
    }
}
impl std::error::Error for ShellStrengthError {}
impl From<PlateError> for ShellStrengthError {
    fn from(value: PlateError) -> Self {
        Self::Native(value)
    }
}
impl From<LaminateStrengthError> for ShellStrengthError {
    fn from(value: LaminateStrengthError) -> Self {
        Self::Material(value)
    }
}

impl PreparedThermalShell<'_> {
    /// Bound the native layer utilization over one whole original facet.
    ///
    /// Prepare the facet's `LaminateThermalField` with
    /// `Arc::clone(strength.thermal().prepared())` from the selected material
    /// object, then use this SAME prepared shell for mechanics and assessment.
    /// Cloning that material object preserves the shared snapshot. Independently
    /// reconstructed equal stiffness is insufficient to establish source,
    /// temperature, free-strain or strength identity and refuses.
    ///
    /// `displacement_6n` is the actual full world nodal state. The exact native
    /// strain observer includes thermally condensed director displacement.
    /// Work is 12 face assessments per layer, no quadrature peak search, heap
    /// allocation, source query or new solve. Cancellation is polled throughout.
    pub fn laminate_strength_bound(
        &self,
        element: usize,
        criterion: LaminateStrengthCriterion<'_>,
        displacement_6n: &[f64],
        maximum_face_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateStrengthBound, ShellStrengthError> {
        let (source, strain) = self.strength_field(
            element,
            criterion,
            displacement_6n,
            6,
            maximum_face_evaluations,
            &mut cancelled,
        )?;
        let mut bound: Option<LaminateStrengthBound> = None;
        for (control, point) in strain.bernstein_controls().into_iter().enumerate() {
            let governing =
                criterion.assess(source, point.membrane, point.gradient_per_m, &mut cancelled)?;
            if bound
                .as_ref()
                .is_none_or(|old| governing.utilization() > old.utilization_upper_bound())
            {
                bound = Some(LaminateStrengthBound { control, governing });
            }
        }
        if cancelled() {
            return Err(PlateError::Cancelled.into());
        }
        Ok(bound.expect("native quadratic triangle has six controls"))
    }

    /// Assess both layer faces at one physical barycentric point of the same
    /// native field. This diagnoses actual local limit crossing separately from
    /// the conservative whole-facet bound. Work is two face evaluations per
    /// layer. It does not establish safety away from this point.
    pub fn laminate_strength_at(
        &self,
        element: usize,
        criterion: LaminateStrengthCriterion<'_>,
        displacement_6n: &[f64],
        barycentric: [f64; 3],
        maximum_face_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShellLayerStrength, ShellStrengthError> {
        let (source, strain) = self.strength_field(
            element,
            criterion,
            displacement_6n,
            1,
            maximum_face_evaluations,
            &mut cancelled,
        )?;
        let point = strain.at_barycentric(barycentric)?;
        criterion.assess(source, point.membrane, point.gradient_per_m, &mut cancelled)
    }

    fn strength_field(
        &self,
        element: usize,
        criterion: LaminateStrengthCriterion<'_>,
        displacement_6n: &[f64],
        points: usize,
        maximum: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(&LaminateThermalField, RealElementStrain), ShellStrengthError> {
        if cancelled() {
            return Err(PlateError::Cancelled.into());
        }
        let facet = self.facets.get(element).ok_or(PlateError::BadSection {
            what: "laminate strength facet out of range",
        })?;
        let source = facet.laminate().ok_or(ShellStrengthError::SourceMismatch)?;
        if !core::ptr::eq(source.source(), criterion.thermal().prepared().as_ref()) {
            return Err(ShellStrengthError::SourceMismatch);
        }
        let required = criterion
            .thermal()
            .layers()
            .len()
            .checked_mul(2)
            .and_then(|faces| faces.checked_mul(points))
            .ok_or(ShellStrengthError::BudgetOverflow)?;
        if required > maximum {
            return Err(ShellStrengthError::FaceBudget { required, maximum });
        }
        if displacement_6n.len() != 6 * self.model.mesh.nodes.len() {
            return Err(PlateError::BadConstraint {
                what: "one physical displacement/rotation per original shell DOF required",
            }
            .into());
        }
        let tri = self.model.mesh.tris[element];
        let nodal = core::array::from_fn(|i| displacement_6n[6 * tri[i / 6] + i % 6]);
        let strain = self.strain_field(element, &nodal, cancelled)?;
        Ok((source, strain))
    }
}
