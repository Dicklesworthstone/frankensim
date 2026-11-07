//! Indentation hardness from actual spatial thermoplastic temperature/history.
//!
//! The original tetrahedral owner supplies complete J2/Chaboche histories and
//! its original enthalpy cells. Every element keeps its independently selected
//! indentation test; no yield-to-hardness relation, averaged temperature or
//! second constitutive solve is introduced. A late source refusal publishes no
//! partial field or coupled candidate. This is an observation, not a new contact
//! law: projecting element values to a surface, assigning an indentation scale,
//! and avoiding duplicate resolved/reduced indentation remain consumer duties.

use super::*;
use fs_matdb::PropertyKey;
use fs_material::hardness::IndentationHardnessStatePoint;
use fs_material::hardness::thermal_history::{
    AccumulatedPlasticHistory, ThermoplasticHardnessError, ThermoplasticHardnessSource,
};
use fs_material::state_point::MaterialPropertySelection;

/// An original source test for each tetrahedron, retained across accepted steps.
#[derive(Clone, Debug)]
pub struct ThermoplasticHardnessFieldBinding<
    'a,
    S: ThermoplasticSourceFamily<'a> = ThermoplasticSource<'a>,
> {
    geometry: Arc<TetLoadGeometry3<'a>>,
    cells: Vec<CellId>,
    sources: Vec<ThermoplasticHardnessSource<'a, S>>,
}

/// One exact source endpoint at an original spatial material point.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermoplasticElementHardness {
    pub element: usize,
    pub cell: CellId,
    /// Contains actual caloric T, native accumulated plastic strain, selected
    /// claim, evidence and the exact original indentation-test convention.
    pub hardness: IndentationHardnessStatePoint,
}

/// Complete field in original tetrahedron order, at one accepted physical time.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermoplasticHardnessField {
    time_s: f64,
    points: Vec<ThermoplasticElementHardness>,
}
impl ThermoplasticHardnessField {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn points(&self) -> &[ThermoplasticElementHardness] {
        &self.points
    }
    /// Borrow the actual per-element pressure-valued response without allocating
    /// another field or discarding its independently retained source receipts.
    pub fn hardness_pa(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.points.iter().map(|point| point.hardness.hardness_pa())
    }
}

/// Original coupled candidate and its complete hardness field, adopted together.
pub struct ThermoplasticHardnessAccepted<C> {
    pub mechanics: C,
    pub hardness: ThermoplasticHardnessField,
}

#[derive(Debug)]
pub enum ThermoplasticHardnessFieldError {
    Native(ThermoplasticError),
    Invalid(&'static str),
    Source {
        element: usize,
        cell: CellId,
        error: ThermoplasticHardnessError,
    },
}
impl From<ThermoplasticError> for ThermoplasticHardnessFieldError {
    fn from(error: ThermoplasticError) -> Self {
        Self::Native(error)
    }
}
impl core::fmt::Display for ThermoplasticHardnessFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial thermoplastic hardness: {self:?}")
    }
}
impl std::error::Error for ThermoplasticHardnessFieldError {}

impl<'a, S: ThermoplasticSourceFamily<'a>> ThermoplasticSolid<'a, S>
where
    <S::Material as SmallStrainLaw>::State: AccumulatedPlasticHistory,
{
    /// Pin one independent indentation test at each accepted original element.
    ///
    /// `maximum_points` bounds source bindings and later output fields. Existing
    /// solid/thermal geometry, clocks, inventory and caloric-law identity checks
    /// admit the entire field before any candidate can be published. The binding
    /// shares the original prepared geometry and owns no mechanical history.
    pub fn bind_hardness_field(
        &self,
        thermal: &DepositionSimulation<'_>,
        selections: &[(PropertyKey, MaterialPropertySelection)],
        maximum_points: usize,
        cx: &Cx<'_>,
    ) -> Result<ThermoplasticHardnessFieldBinding<'a, S>, ThermoplasticHardnessFieldError> {
        poll(cx)?;
        if maximum_points == 0
            || self.cells.len() > maximum_points
            || selections.len() != self.cells.len()
            || thermal.time_s() != self.time_s
        {
            return Err(ThermoplasticHardnessFieldError::Invalid(
                "matching accepted clocks and one bounded hardness selection per original tetrahedron required",
            ));
        }
        let mut cells = Vec::with_capacity(self.cells.len());
        let mut sources = Vec::with_capacity(self.cells.len());
        for (element, (bound, (key, selection))) in self.cells.iter().zip(selections).enumerate() {
            poll(cx)?;
            self.checked_thermal_cell::<Infallible>(bound, thermal, true)?;
            let original = bound.source.thermal_source();
            let caloric = original
                .caloric_source()
                .curve()
                .state_at_specific_enthalpy(bound.accepted.specific_h_j_kg)
                .map_err(|error| ThermoplasticError::Material(FeedbackError::Caloric(error)))?;
            sources.push(
                ThermoplasticHardnessSource::<'a, S>::try_from_native_state(
                    original,
                    &self.materials[element],
                    caloric,
                    &self.solution.element_states[element],
                    key.clone(),
                    selection.clone(),
                )
                .map_err(|error| ThermoplasticHardnessFieldError::Source {
                    element,
                    cell: bound.accepted.id,
                    error,
                })?,
            );
            cells.push(bound.accepted.id);
        }
        poll(cx)?;
        Ok(ThermoplasticHardnessFieldBinding {
            geometry: Arc::clone(&self.geometry),
            cells,
            sources,
        })
    }
}

impl<'a, S: ThermoplasticSourceFamily<'a>> ThermoplasticHardnessFieldBinding<'a, S>
where
    <S::Material as SmallStrainLaw>::State: AccumulatedPlasticHistory,
{
    /// Observe the exact accepted original mesh, enthalpy inventory and histories.
    /// Source support is checked at each actual endpoint; interpolation across
    /// elements or coverage between physical times is not asserted.
    pub fn at_accepted(
        &self,
        solid: &ThermoplasticSolid<'a, S>,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> Result<ThermoplasticHardnessField, ThermoplasticHardnessFieldError> {
        poll(cx)?;
        if !Arc::ptr_eq(&self.geometry, &solid.geometry)
            || self.cells.len() != solid.cells.len()
            || thermal.time_s() != solid.time_s
            || self
                .cells
                .iter()
                .zip(&solid.cells)
                .any(|(id, bound)| *id != bound.accepted.id)
        {
            return Err(ThermoplasticHardnessFieldError::Invalid(
                "hardness field requires the original mesh, ordered material cells and accepted clock",
            ));
        }
        let mut points = Vec::with_capacity(self.cells.len());
        for (element, bound) in solid.cells.iter().enumerate() {
            poll(cx)?;
            solid.checked_thermal_cell::<Infallible>(bound, thermal, true)?;
            let original = bound.source.thermal_source();
            let caloric = original
                .caloric_source()
                .curve()
                .state_at_specific_enthalpy(bound.accepted.specific_h_j_kg)
                .map_err(|error| ThermoplasticError::Material(FeedbackError::Caloric(error)))?;
            let hardness = self.sources[element]
                .at_native_state(
                    original,
                    &solid.materials[element],
                    caloric,
                    &solid.solution.element_states[element],
                )
                .map_err(|error| ThermoplasticHardnessFieldError::Source {
                    element,
                    cell: bound.accepted.id,
                    error,
                })?;
            points.push(ThermoplasticElementHardness {
                element,
                cell: bound.accepted.id,
                hardness,
            });
        }
        poll(cx)?;
        Ok(ThermoplasticHardnessField {
            time_s: solid.time_s,
            points,
        })
    }

    /// Gate an existing original coupled candidate without rerunning mechanics,
    /// heat or observers. A missing hardness endpoint consumes the private
    /// proposal and returns no partial field or advanced physical owners.
    pub fn check_candidate<'t, O>(
        &self,
        candidate: ThermoplasticCandidate<'a, 't, O, S>,
        cx: &Cx<'_>,
    ) -> Result<
        ThermoplasticHardnessAccepted<ThermoplasticCandidate<'a, 't, O, S>>,
        ThermoplasticHardnessFieldError,
    > {
        let hardness = self.at_accepted(candidate.solid(), candidate.thermal(), cx)?;
        Ok(ThermoplasticHardnessAccepted {
            mechanics: candidate,
            hardness,
        })
    }
}
