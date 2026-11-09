//! Separately sourced hardness on the eight native histories of each Hex8 cell.
//!
//! Standard and B-bar cells retain the assumed strains actually used by their
//! constitutive owner. Each original J2, Voce or Chaboche history is observed
//! independently: neither strain nor hardness is averaged across Gauss points.
//! Observation performs no return map, source query or per-point source clone.
//! The caller accepts the mechanical candidate and its hardness field together.
//!
//! These are observations of caller-supplied accepted endpoints, not proof of a
//! loading trajectory or mesh origin. Hardness changes add no stress, stiffness,
//! stored energy or heat. Compact exports transfer the property only; a contact
//! consumer owns its explicit cell-to-point map and half-space approximation.

use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_material::hardness::cyclic_history::{
    CyclicHardnessError, ResolvedCyclicHardness, SourceBoundHardnessLaw,
};
use fs_material::hardness::history_field::{AcceptedHardnessField, HardnessFieldError};
use fs_material::hardness::plastic_history::{J2HardnessError, ResolvedJ2Hardness};
use fs_material::plastic::J2State;
use fs_material::plastic::source::ResolvedJ2Plasticity;

use super::{SmallStrainHexProblem, SmallStrainHexSolution};

/// Original prepared sources in native element/Gauss-point order. In each
/// eight-point group x is outermost and z innermost, as in the native solver.
#[derive(Debug)]
pub enum HexHardnessField<'a, H> {
    /// One prepared source shared by the entire mesh.
    Uniform(&'a H),
    /// One prepared source per cell, shared by its eight histories.
    PerElement(&'a [&'a H]),
    /// Eight independent original sources per cell in native quadrature order.
    PerQuadrature(&'a [[&'a H; 8]]),
}

impl<H> Copy for HexHardnessField<'_, H> {}
impl<H> Clone for HexHardnessField<'_, H> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<H> HexHardnessField<'_, H> {
    fn covers(&self, elements: usize) -> bool {
        match self {
            Self::Uniform(_) => true,
            Self::PerElement(sources) => sources.len() == elements,
            Self::PerQuadrature(sources) => sources.len() == elements,
        }
    }

    fn at(&self, element: usize, point: usize) -> &H {
        match self {
            Self::Uniform(source) => source,
            Self::PerElement(sources) => sources[element],
            Self::PerQuadrature(sources) => sources[element][point],
        }
    }
}

/// No partial field is published on source, native-history or cancellation
/// refusal. Point numbers are local to the named element, from zero to seven.
#[derive(Clone, Debug, PartialEq)]
pub enum HexHardnessError {
    /// Source shape, solution order or native allocation bounds do not match.
    Invalid(&'static str),
    /// The original J2 source or complete endpoint history refused this point.
    J2 {
        element: usize,
        point: usize,
        error: J2HardnessError,
    },
    /// The original cyclic source or complete branch history refused this point.
    Cyclic {
        element: usize,
        point: usize,
        error: CyclicHardnessError,
    },
    /// Compact export refusal, retaining the material owner's flattened index.
    Field(HardnessFieldError),
    /// Cancellation observed before publication of a complete field.
    Cancelled,
}

impl core::fmt::Display for HexHardnessError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "Hex8 hardness: {reason}"),
            Self::J2 { element, point, error } => {
                write!(f, "Hex8 J2 hardness element {element} point {point}: {error}")
            }
            Self::Cyclic { element, point, error } => {
                write!(f, "Hex8 cyclic hardness element {element} point {point}: {error}")
            }
            Self::Field(error) => write!(f, "Hex8 hardness field: {error}"),
            Self::Cancelled => f.write_str("Hex8 hardness observation cancelled"),
        }
    }
}
impl std::error::Error for HexHardnessError {}
impl From<HardnessFieldError> for HexHardnessError {
    fn from(error: HardnessFieldError) -> Self {
        match error {
            HardnessFieldError::Interrupted => Self::Cancelled,
            error => Self::Field(error),
        }
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), HexHardnessError> {
    cx.checkpoint().map_err(|_| HexHardnessError::Cancelled)
}

fn validate<L: SmallStrainLaw>(
    problem: &SmallStrainHexProblem<'_, L>,
    solution: &SmallStrainHexSolution<L::State>,
    cx: &Cx<'_>,
) -> Result<usize, HexHardnessError> {
    poll(cx)?;
    let elements = problem.hexahedra.len();
    let points = elements.checked_mul(8).ok_or(HexHardnessError::Invalid(
        "native quadrature count overflow",
    ))?;
    if elements == 0
        || elements > problem.budget.maximum_elements
        || points > problem.budget.maximum_quadrature_points
        || problem.nodes_m.is_empty()
        || problem.nodes_m.len() > problem.budget.maximum_nodes
        || problem.laws.validate_quadrature_counts(std::iter::repeat_n(8, elements)).is_err()
        || solution.displacement_m.len() != problem.nodes_m.len()
        || solution.quadrature_strain.len() != points
        || solution.quadrature_stress_pa.len() != points
        || solution.quadrature_states.len() != points
        || solution.quadrature_volume_m3.len() != points
        || solution.element_offsets.len() != elements + 1
    {
        return Err(HexHardnessError::Invalid(
            "native solution must cover the original bounded mesh with eight histories per cell",
        ));
    }
    for (element, &offset) in solution.element_offsets.iter().enumerate() {
        if element % 32 == 0 {
            poll(cx)?;
        }
        if offset != 8 * element {
            return Err(HexHardnessError::Invalid(
                "native element/Gauss-point order must remain unchanged",
            ));
        }
    }
    Ok(points)
}

impl SmallStrainHexProblem<'_, ResolvedJ2Plasticity> {
    /// Observe original H(alpha) at every accepted native J2 point. Source/law
    /// identity, complete history, instantaneous yield and the prepared alpha
    /// interval pass the original material gate. The yield allowance is in Pa.
    /// One scalar per point returns in native order only after every point passes.
    pub fn observe_hardness(
        &self,
        solution: &SmallStrainHexSolution<J2State>,
        sources: HexHardnessField<'_, ResolvedJ2Hardness>,
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HexHardnessError> {
        let points = validate(self, solution, cx)?;
        if !sources.covers(self.hexahedra.len()) {
            return Err(HexHardnessError::Invalid("hardness sources must cover every native cell"));
        }
        let mut hardness = Vec::with_capacity(points);
        for index in 0..points {
            if index % 256 == 0 {
                poll(cx)?;
            }
            let (element, point) = (index / 8, index % 8);
            hardness.push(sources.at(element, point).hardness_pa_at_accepted(
                self.laws.at_point(element, point),
                &solution.quadrature_strain[index],
                &solution.quadrature_states[index],
                yield_tolerance_pa,
            ).map_err(|error| HexHardnessError::J2 { element, point, error })?);
        }
        poll(cx)?;
        Ok(hardness)
    }

    /// Export a fixed original source and compact accepted alpha/H samples.
    /// Every actual quadrature law must match this source; a heterogeneous
    /// solid cannot be relabeled homogeneous. Contact uses the exported native
    /// point indices explicitly, without interpolation or hardness averaging.
    pub fn export_hardness_field(
        &self,
        solution: &SmallStrainHexSolution<J2State>,
        source: &ResolvedJ2Hardness,
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<AcceptedHardnessField, HexHardnessError> {
        let points = validate(self, solution, cx)?;
        let field = AcceptedHardnessField::from_j2_checked(
            source,
            (0..points).map(|index| (
                self.laws.at_point(index / 8, index % 8),
                &solution.quadrature_strain[index],
                &solution.quadrature_states[index],
            )),
            yield_tolerance_pa,
            || cx.checkpoint().map_err(|_| HardnessFieldError::Interrupted),
        )?;
        poll(cx)?;
        Ok(field)
    }
}

impl<L: SourceBoundHardnessLaw> SmallStrainHexProblem<'_, L> {
    /// Observe H(alpha) on complete native Voce/Chaboche histories, preserving
    /// each source's nonlinear yield surface and original ordered backstresses.
    /// Standard and B-bar assumed strains retain their native interpretation.
    pub fn observe_cyclic_hardness(
        &self,
        solution: &SmallStrainHexSolution<L::State>,
        sources: HexHardnessField<'_, ResolvedCyclicHardness<L>>,
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HexHardnessError> {
        let points = validate(self, solution, cx)?;
        if !sources.covers(self.hexahedra.len()) {
            return Err(HexHardnessError::Invalid("hardness sources must cover every native cell"));
        }
        let mut hardness = Vec::with_capacity(points);
        for index in 0..points {
            if index % 256 == 0 {
                poll(cx)?;
            }
            let (element, point) = (index / 8, index % 8);
            hardness.push(sources.at(element, point).hardness_pa_at_accepted(
                self.laws.at_point(element, point),
                &solution.quadrature_strain[index],
                &solution.quadrature_states[index],
                yield_tolerance_pa,
            ).map_err(|error| HexHardnessError::Cyclic { element, point, error })?);
        }
        poll(cx)?;
        Ok(hardness)
    }

    /// Compact fixed-source cyclic field in native quadrature order. Complete
    /// branch histories are admitted by their existing material owner and are
    /// not cloned into the output. The original source receipts are shared once.
    pub fn export_cyclic_hardness_field(
        &self,
        solution: &SmallStrainHexSolution<L::State>,
        source: &ResolvedCyclicHardness<L>,
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<AcceptedHardnessField, HexHardnessError> {
        let points = validate(self, solution, cx)?;
        let field = AcceptedHardnessField::from_cyclic_checked(
            source,
            (0..points).map(|index| (
                self.laws.at_point(index / 8, index % 8),
                &solution.quadrature_strain[index],
                &solution.quadrature_states[index],
            )),
            yield_tolerance_pa,
            || cx.checkpoint().map_err(|_| HardnessFieldError::Interrupted),
        )?;
        poll(cx)?;
        Ok(field)
    }
}
