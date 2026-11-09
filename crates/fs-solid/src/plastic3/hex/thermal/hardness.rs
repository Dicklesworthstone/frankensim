//! Original indentation hardness on accepted native thermal Hex8 histories.
//!
//! The eight material points retain their own caloric temperature, assumed
//! strain and complete J2/Voce/Chaboche history. One exact indentation test per
//! original cell is pinned across observations. No constitutive update, averaged
//! temperature, yield-to-hardness conversion or additional heat is introduced.
//! Full typed endpoints retain source receipts; observing them is not a
//! certification of coverage between endpoints or of a supplied trajectory.

use super::{ThermalHexProblem, ThermalHexStep};
use crate::Formulation;
use fs_exec::Cx;
use fs_matdb::PropertyKey;
use fs_material::SmallStrainLaw;
use fs_material::hardness::IndentationHardnessStatePoint;
use fs_material::hardness::thermal_history::{
    AccumulatedPlasticHistory, ThermoplasticHardnessError, ThermoplasticHardnessSource,
};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::plastic::source::ResolvedJ2Plasticity;
use fs_material::plastic::thermal::family::{ThermoplasticFamilyError, ThermoplasticSourceFamily};
use fs_material::state_point::MaterialPropertySelection;

type State<'a, F> = <<F as ThermoplasticSourceFamily<'a>>::Material as SmallStrainLaw>::State;

/// Original mesh/families and one pinned indentation test per cell. The binding
/// borrows its immutable owners; moving to a different mesh requires rebinding.
#[derive(Clone, Debug)]
pub struct ThermalHexHardnessBinding<'p, 'a, F: ThermoplasticSourceFamily<'a>> {
    nodes_m: &'p [[f64; 3]],
    hexahedra: &'p [[usize; 8]],
    families: &'p [F],
    formulation: Formulation,
    sources: Vec<ThermoplasticHardnessSource<'a, F>>,
}

/// Complete typed field in native `8 * element + local_point` order. This
/// exports hardness only, not residual stress, geometry, heat or contact work.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalHexHardnessField {
    points: Vec<IndentationHardnessStatePoint>,
}
impl ThermalHexHardnessField {
    /// Original tests, selected claims and actual T/alpha at all native points.
    pub fn points(&self) -> &[IndentationHardnessStatePoint] {
        &self.points
    }

    /// Borrow the pointwise pressure-valued hardness without another allocation.
    pub fn hardness_pa(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.points.iter().map(IndentationHardnessStatePoint::hardness_pa)
    }
}

/// Refusal publishes no partial field and changes no native state or account.
#[derive(Debug)]
pub enum ThermalHexHardnessError {
    Invalid(&'static str),
    Native {
        /// Flattened native material-point index.
        point: usize,
        error: ThermoplasticFamilyError,
    },
    Source {
        /// Flattened native material-point index.
        point: usize,
        error: ThermoplasticHardnessError,
    },
    Cancelled,
}
impl core::fmt::Display for ThermalHexHardnessError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal Hex8 hardness: {self:?}")
    }
}
impl std::error::Error for ThermalHexHardnessError {}

fn poll(cx: &Cx<'_>) -> Result<(), ThermalHexHardnessError> {
    cx.checkpoint().map_err(|_| ThermalHexHardnessError::Cancelled)
}

fn validate<'a, F: ThermoplasticSourceFamily<'a>>(
    problem: &ThermalHexProblem<'_, F>,
    step: &ThermalHexStep<State<'a, F>>,
    yield_tolerance_pa: f64,
    cx: &Cx<'_>,
) -> Result<usize, ThermalHexHardnessError> {
    poll(cx)?;
    let elements = problem.hexahedra.len();
    let count = elements.checked_mul(8).ok_or(ThermalHexHardnessError::Invalid(
        "native quadrature count overflow",
    ))?;
    let solution = &step.solution;
    if elements == 0
        || elements > problem.budget.maximum_elements
        || count > problem.budget.maximum_quadrature_points
        || problem.nodes_m.is_empty()
        || problem.nodes_m.len() > problem.budget.maximum_nodes
        || problem.families.len() != elements
        || solution.displacement_m.len() != problem.nodes_m.len()
        || solution.quadrature_strain.len() != count
        || solution.quadrature_stress_pa.len() != count
        || solution.quadrature_states.len() != count
        || solution.quadrature_volume_m3.len() != count
        || solution.element_offsets.len() != elements + 1
        || step.specific_enthalpy_j_kg.len() != count
        || step.quadrature_temperature_k.len() != count
        || step.quadrature_energy.len() != count
        || !yield_tolerance_pa.is_finite()
        || yield_tolerance_pa <= 0.0
    {
        return Err(ThermalHexHardnessError::Invalid(
            "bounded original thermal mesh, eight ordered histories per cell and positive yield tolerance required",
        ));
    }
    for (element, &offset) in solution.element_offsets.iter().enumerate() {
        if element % 32 == 0 { poll(cx)?; }
        if offset != 8 * element {
            return Err(ThermalHexHardnessError::Invalid(
                "native element/Gauss-point ordering must remain unchanged",
            ));
        }
    }
    Ok(count)
}

fn accepted_point<'a, F: ThermoplasticSourceFamily<'a>>(
    problem: &ThermalHexProblem<'_, F>,
    step: &ThermalHexStep<State<'a, F>>,
    point: usize,
    yield_tolerance_pa: f64,
) -> Result<(EquilibriumPhaseState, ResolvedJ2Plasticity), ThermalHexHardnessError> {
    let family = &problem.families[point / 8];
    let h = step.specific_enthalpy_j_kg[point];
    // The family admits its COMPLETE native state, including cyclic branches.
    // It does not advance a time interval or project onto a bare J2 history.
    family.validate_accepted(
        h,
        &step.solution.quadrature_strain[point],
        &step.solution.quadrature_states[point],
        yield_tolerance_pa,
    ).map_err(|error| ThermalHexHardnessError::Native { point, error })?;
    let (caloric, material) = family.thermal_source().resolve_enthalpy(h)
        .map_err(|error| ThermalHexHardnessError::Native { point, error: error.into() })?;
    if caloric.temperature_k() != step.quadrature_temperature_k[point] {
        return Err(ThermalHexHardnessError::Invalid(
            "hardness must observe the actual original caloric point temperature",
        ));
    }
    Ok((caloric, material))
}

impl<'p, F> ThermalHexProblem<'p, F> {
    /// Pin one exact original indentation test per native cell. The first
    /// accepted point of each cell selects the claim; `at_accepted` admits ALL
    /// points before publishing a field. Later observations keep that selection.
    pub fn bind_hardness_field<'a>(
        &self,
        step: &ThermalHexStep<State<'a, F>>,
        selections: &[(PropertyKey, MaterialPropertySelection)],
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalHexHardnessBinding<'p, 'a, F>, ThermalHexHardnessError>
    where
        F: ThermoplasticSourceFamily<'a>,
        State<'a, F>: AccumulatedPlasticHistory,
    {
        validate(self, step, yield_tolerance_pa, cx)?;
        if selections.len() != self.hexahedra.len() {
            return Err(ThermalHexHardnessError::Invalid(
                "one original hardness selection per native cell required",
            ));
        }
        let mut sources = Vec::with_capacity(selections.len());
        for (element, (key, selection)) in selections.iter().enumerate() {
            poll(cx)?;
            let point = 8 * element;
            let (caloric, material) = accepted_point(self, step, point, yield_tolerance_pa)?;
            sources.push(ThermoplasticHardnessSource::<'a, F>::try_from_native_state(
                self.families[element].thermal_source(),
                &material,
                caloric,
                &step.solution.quadrature_states[point],
                key.clone(),
                selection.clone(),
            ).map_err(|error| ThermalHexHardnessError::Source { point, error })?);
        }
        poll(cx)?;
        Ok(ThermalHexHardnessBinding {
            nodes_m: self.nodes_m,
            hexahedra: self.hexahedra,
            families: self.families,
            formulation: self.formulation,
            sources,
        })
    }
}

impl<'a, F: ThermoplasticSourceFamily<'a>> ThermalHexHardnessBinding<'_, 'a, F>
where
    State<'a, F>: AccumulatedPlasticHistory,
{
    /// Observe every actual accepted point under the original mesh, source
    /// families, quadrature order and formulation. Loading/support values may
    /// change between accepted steps. Source resolution and native history
    /// validation use their existing owners without another constitutive solve.
    ///
    /// The caller accepts the native thermal/mechanical candidate and this field
    /// together. Full typed endpoints allocate their source receipts; no compact
    /// fixed-temperature half-space transfer or coverage between times is claimed.
    pub fn at_accepted(
        &self,
        problem: &ThermalHexProblem<'_, F>,
        step: &ThermalHexStep<State<'a, F>>,
        yield_tolerance_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalHexHardnessField, ThermalHexHardnessError> {
        let count = validate(problem, step, yield_tolerance_pa, cx)?;
        if !core::ptr::eq(self.nodes_m, problem.nodes_m)
            || !core::ptr::eq(self.hexahedra, problem.hexahedra)
            || !core::ptr::eq(self.families, problem.families)
            || self.formulation != problem.formulation
        {
            return Err(ThermalHexHardnessError::Invalid(
                "hardness binding requires the original mesh, ordered families and formulation",
            ));
        }
        let mut points = Vec::with_capacity(count);
        for point in 0..count {
            poll(cx)?;
            let (caloric, material) = accepted_point(problem, step, point, yield_tolerance_pa)?;
            points.push(self.sources[point / 8].at_native_state(
                self.families[point / 8].thermal_source(),
                &material,
                caloric,
                &step.solution.quadrature_states[point],
            ).map_err(|error| ThermalHexHardnessError::Source { point, error })?);
        }
        poll(cx)?;
        Ok(ThermalHexHardnessField { points })
    }
}
