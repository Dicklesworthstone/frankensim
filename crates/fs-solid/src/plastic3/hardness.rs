//! Separately sourced hardness on accepted native tetrahedral J2 histories.
//!
//! The existing solid solver owns equilibrium, plastic strain and stress. This
//! observer reads each element's own accumulated strain and original fixed-state
//! material source without another return map. A complete scalar field is
//! published only after every endpoint passes native history and source-support
//! checks. Callers accept the mechanical candidate and this field together.
//!
//! Hardness is not averaged over regions or used to alter yield stress, elastic
//! stiffness, fracture resistance or heat. Endpoint observation does not certify
//! a supplied solution's origin or unobserved intermediate mechanical states.

use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_material::hardness::cyclic_history::{
    CyclicHardnessError, ResolvedCyclicHardness, SourceBoundHardnessLaw,
};
use fs_material::hardness::history_field::{AcceptedHardnessField, HardnessFieldError};
use fs_material::hardness::plastic_history::{J2HardnessError, ResolvedJ2Hardness};
use fs_material::plastic::J2State;
use fs_material::plastic::source::ResolvedJ2Plasticity;

use super::{SmallStrainTetProblem, SmallStrainTetSolution};

/// Original prepared hardness sources in the native element connectivity order.
/// Borrowing one source at multiple elements never copies its data or receipts.
#[derive(Clone, Copy, Debug)]
pub enum J2HardnessField<'a> {
    Uniform(&'a ResolvedJ2Hardness),
    PerElement(&'a [&'a ResolvedJ2Hardness]),
}

/// No partial hardness field is returned on any source, history or cancellation
/// refusal. The mechanical solution and prepared sources remain unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum J2HardnessObservationError {
    Invalid(&'static str),
    Element {
        element: usize,
        error: J2HardnessError,
    },
    Cancelled,
}

impl core::fmt::Display for J2HardnessObservationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "solid J2 hardness observation: {reason}"),
            Self::Element { element, error } => {
                write!(f, "solid J2 hardness element {element}: {error}")
            }
            Self::Cancelled => f.write_str("solid J2 hardness observation cancelled"),
        }
    }
}
impl std::error::Error for J2HardnessObservationError {}

impl SmallStrainTetProblem<'_, ResolvedJ2Plasticity> {
    /// Evaluate the original hardness response at each accepted native element.
    ///
    /// Both uniform and heterogeneous laws use the same native solution. Source
    /// identity, complete J2 history, finite small-strain bounds, instantaneous
    /// yield and the prepared alpha range are checked at every point. The yield
    /// allowance is an explicit positive absolute numerical tolerance in Pa.
    ///
    /// The returned allocation contains only one f64 per element. Source
    /// receipts and native histories are not cloned. A selected endpoint can be
    /// exported through `ResolvedJ2Hardness::observe_accepted` followed by its
    /// `resolve_hardness_state` method when a contact/wear consumer needs it.
    pub fn observe_hardness(
        &self,
        solution: &SmallStrainTetSolution<J2State>,
        sources: J2HardnessField<'_>,
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, J2HardnessObservationError> {
        let poll = || {
            cx.checkpoint()
                .map_err(|_| J2HardnessObservationError::Cancelled)
        };
        poll()?;
        let count = self.tetrahedra.len();
        if count == 0
            || count > self.budget.maximum_tetrahedra
            || self.nodes_m.len() > self.budget.maximum_nodes
            || self.laws.validate_count(count).is_err()
            || solution.element_strain.len() != count
            || solution.element_stress.len() != count
            || solution.element_states.len() != count
            || solution.displacement_m.len() != self.nodes_m.len()
            || matches!(sources, J2HardnessField::PerElement(values) if values.len() != count)
        {
            return Err(J2HardnessObservationError::Invalid(
                "hardness sources and native solution must cover the original admitted mesh",
            ));
        }
        let mut hardness_pa = Vec::with_capacity(count);
        for element in 0..count {
            if element % 256 == 0 {
                poll()?;
            }
            let source = match sources {
                J2HardnessField::Uniform(source) => source,
                J2HardnessField::PerElement(values) => values[element],
            };
            hardness_pa.push(
                source
                    .hardness_pa_at_accepted(
                        self.laws.at(element),
                        &solution.element_strain[element],
                        &solution.element_states[element],
                        yield_tolerance_pa,
                    )
                    .map_err(|error| J2HardnessObservationError::Element { element, error })?,
            );
        }
        poll()?;
        Ok(hardness_pa)
    }
}

/// Prepared original cyclic hardness in native element connectivity order.
/// The complete native source law, including every ordered backstress branch,
/// must match the actual element law retained by the solid problem.
#[derive(Clone, Debug)]
pub enum CyclicHardnessField<'a, L: SourceBoundHardnessLaw> {
    Uniform(&'a ResolvedCyclicHardness<L>),
    PerElement(&'a [&'a ResolvedCyclicHardness<L>]),
}

/// Atomic cyclic-hardness field refusal; native solution and sources are read-only.
#[derive(Clone, Debug, PartialEq)]
pub enum CyclicHardnessObservationError {
    Invalid(&'static str),
    Element {
        element: usize,
        error: CyclicHardnessError,
    },
    Cancelled,
}

impl core::fmt::Display for CyclicHardnessObservationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "solid cyclic hardness observation: {reason}"),
            Self::Element { element, error } => {
                write!(f, "solid cyclic hardness element {element}: {error}")
            }
            Self::Cancelled => f.write_str("solid cyclic hardness observation cancelled"),
        }
    }
}
impl std::error::Error for CyclicHardnessObservationError {}

impl<L: SourceBoundHardnessLaw> SmallStrainTetProblem<'_, L> {
    /// Observe complete native Voce/Chaboche endpoints at every element.
    /// Original branch order, nonlinear yield, history bounds and continuous
    /// hardness support are retained. Observation performs no constitutive
    /// update or per-element source allocation; the result has one scalar per
    /// element and is returned only when every point and cancellation gate pass.
    pub fn observe_cyclic_hardness(
        &self,
        solution: &SmallStrainTetSolution<L::State>,
        sources: CyclicHardnessField<'_, L>,
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, CyclicHardnessObservationError> {
        let poll = || {
            cx.checkpoint()
                .map_err(|_| CyclicHardnessObservationError::Cancelled)
        };
        poll()?;
        let count = self.tetrahedra.len();
        if count == 0
            || count > self.budget.maximum_tetrahedra
            || self.nodes_m.len() > self.budget.maximum_nodes
            || self.laws.validate_count(count).is_err()
            || solution.element_strain.len() != count
            || solution.element_stress.len() != count
            || solution.element_states.len() != count
            || solution.displacement_m.len() != self.nodes_m.len()
            || matches!(&sources, CyclicHardnessField::PerElement(values) if values.len() != count)
        {
            return Err(CyclicHardnessObservationError::Invalid(
                "hardness sources and native solution must cover the original admitted mesh",
            ));
        }
        let mut hardness_pa = Vec::with_capacity(count);
        for element in 0..count {
            if element % 256 == 0 {
                poll()?;
            }
            let source = match &sources {
                CyclicHardnessField::Uniform(source) => *source,
                CyclicHardnessField::PerElement(values) => values[element],
            };
            hardness_pa.push(
                source
                    .hardness_pa_at_accepted(
                        self.laws.at(element),
                        &solution.element_strain[element],
                        &solution.element_states[element],
                        yield_tolerance_pa,
                    )
                    .map_err(|error| CyclicHardnessObservationError::Element { element, error })?,
            );
        }
        poll()?;
        Ok(hardness_pa)
    }
}

/// A complete source-bound field export refused without changing the native
/// solution, its material owners or any previously exported field.
#[derive(Clone, Debug, PartialEq)]
pub enum HardnessFieldExportError {
    Invalid(&'static str),
    Field(HardnessFieldError),
    Cancelled,
}

impl core::fmt::Display for HardnessFieldExportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "solid hardness field export: {reason}"),
            Self::Field(error) => write!(f, "solid hardness field export: {error}"),
            Self::Cancelled => f.write_str("solid hardness field export cancelled"),
        }
    }
}
impl std::error::Error for HardnessFieldExportError {}
impl From<HardnessFieldError> for HardnessFieldExportError {
    fn from(error: HardnessFieldError) -> Self {
        match error {
            HardnessFieldError::Interrupted => Self::Cancelled,
            error => Self::Field(error),
        }
    }
}

fn validate_field_export<L: SmallStrainLaw>(
    problem: &SmallStrainTetProblem<'_, L>,
    solution: &SmallStrainTetSolution<L::State>,
    cx: &Cx<'_>,
) -> Result<(), HardnessFieldExportError> {
    cx.checkpoint()
        .map_err(|_| HardnessFieldExportError::Cancelled)?;
    let count = problem.tetrahedra.len();
    if count == 0
        || count > problem.budget.maximum_tetrahedra
        || problem.nodes_m.len() > problem.budget.maximum_nodes
        || problem.laws.validate_count(count).is_err()
        || solution.element_strain.len() != count
        || solution.element_stress.len() != count
        || solution.element_states.len() != count
        || solution.displacement_m.len() != problem.nodes_m.len()
    {
        return Err(HardnessFieldExportError::Invalid(
            "native solution must cover the original admitted mesh in element order",
        ));
    }
    Ok(())
}

impl SmallStrainTetProblem<'_, ResolvedJ2Plasticity> {
    /// Export one original source and compact accepted alpha/H samples in the
    /// native tetrahedron order. Every actual element law must match the fixed
    /// source; a heterogeneous material field cannot be relabeled homogeneous.
    ///
    /// This property transfer does not export residual stress, displacements or
    /// geometry. A contact consumer owns its explicit cell-to-element map and
    /// homogeneous half-space approximation, without duplicating indentation.
    pub fn export_hardness_field(
        &self,
        solution: &SmallStrainTetSolution<J2State>,
        source: &ResolvedJ2Hardness,
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<AcceptedHardnessField, HardnessFieldExportError> {
        validate_field_export(self, solution, cx)?;
        let field = AcceptedHardnessField::from_j2_checked(
            source,
            (0..self.tetrahedra.len()).map(|element| {
                (
                    self.laws.at(element),
                    &solution.element_strain[element],
                    &solution.element_states[element],
                )
            }),
            yield_tolerance_pa,
            || cx.checkpoint().map_err(|_| HardnessFieldError::Interrupted),
        )?;
        cx.checkpoint()
            .map_err(|_| HardnessFieldExportError::Cancelled)?;
        Ok(field)
    }
}

impl<L: SourceBoundHardnessLaw> SmallStrainTetProblem<'_, L> {
    /// Export a fixed-source accepted Voce/Chaboche property field. Complete
    /// actual branch histories and source roles pass the native material gate
    /// once per point. Histories and source receipts are not cloned per element.
    pub fn export_cyclic_hardness_field(
        &self,
        solution: &SmallStrainTetSolution<L::State>,
        source: &ResolvedCyclicHardness<L>,
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<AcceptedHardnessField, HardnessFieldExportError> {
        validate_field_export(self, solution, cx)?;
        let field = AcceptedHardnessField::from_cyclic_checked(
            source,
            (0..self.tetrahedra.len()).map(|element| {
                (
                    self.laws.at(element),
                    &solution.element_strain[element],
                    &solution.element_states[element],
                )
            }),
            yield_tolerance_pa,
            || cx.checkpoint().map_err(|_| HardnessFieldError::Interrupted),
        )?;
        cx.checkpoint()
            .map_err(|_| HardnessFieldExportError::Cancelled)?;
        Ok(field)
    }
}
