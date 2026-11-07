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
