//! Hardness admitted together with the original coupled Hex8 candidate.
//!
//! The existing owner retains all eight native histories, exact point enthalpies,
//! reference inventories and the converged heat/mechanical accounts. This gate
//! observes those private accepted states through their original material source;
//! it neither reconstructs a thermal step nor repeats mechanics or transport.
//! One exact indentation test per element remains pinned across continuation.
//! Coverage is checked at actual endpoints, not between physical times.

use core::convert::Infallible;
use std::sync::Arc;

use super::feedback::ThermalBinding;
use super::{
    HexThermoplasticCandidate, HexThermoplasticError, PointCaloric, ThermoplasticHexSolid, poll,
    source_error,
};
use fs_exec::Cx;
use fs_matdb::PropertyKey;
use fs_material::SmallStrainLaw;
use fs_material::hardness::IndentationHardnessStatePoint;
use fs_material::hardness::thermal_history::{
    AccumulatedPlasticHistory, ThermoplasticHardnessError, ThermoplasticHardnessSource,
};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::plastic::section::thermal::source::ThermalSectionSourceError;
use fs_material::plastic::source::ResolvedJ2Plasticity;
use fs_material::plastic::thermal::family::ThermoplasticSourceFamily;
use fs_material::state_point::MaterialPropertySelection;
use fs_solid::Formulation;
use fs_solid::hyper3::plastic::FinitePlasticLawField;
use fs_thermal::enthalpy::{CellId, DepositionSimulation};

/// Original source/mesh owners and one exact indentation test per element.
/// The shared family allocation survives adoption of successive coupled states.
#[derive(Clone, Debug)]
pub struct HexThermoplasticHardnessBinding<'a, F: ThermoplasticSourceFamily<'a>> {
    nodes_m: &'a [[f64; 3]],
    hexahedra: &'a [[usize; 8]],
    families: Arc<[F]>,
    formulation: Formulation,
    cells: Vec<CellId>,
    sources: Vec<ThermoplasticHardnessSource<'a, F>>,
}

/// One original material point and its complete indentation source receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct HexThermoplasticPointHardness {
    pub element: usize,
    pub local_point: usize,
    pub cell: CellId,
    /// Actual caloric temperature and native accumulated plastic strain.
    pub hardness: IndentationHardnessStatePoint,
}

/// Complete field in original element-major, eight-point Gauss order.
#[derive(Clone, Debug, PartialEq)]
pub struct HexThermoplasticHardnessField {
    time_s: f64,
    points: Vec<HexThermoplasticPointHardness>,
}
impl HexThermoplasticHardnessField {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn points(&self) -> &[HexThermoplasticPointHardness] {
        &self.points
    }
    pub fn hardness_pa(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.points.iter().map(|point| point.hardness.hardness_pa())
    }
}

/// Adopt the unchanged coupled candidate and its hardness field together.
pub struct HexThermoplasticHardnessAccepted<C> {
    pub mechanics: C,
    pub hardness: HexThermoplasticHardnessField,
}

#[derive(Debug)]
pub enum HexThermoplasticHardnessError {
    Native(HexThermoplasticError),
    Invalid(&'static str),
    Source {
        point: usize,
        cell: CellId,
        error: ThermoplasticHardnessError,
    },
}
impl From<HexThermoplasticError> for HexThermoplasticHardnessError {
    fn from(error: HexThermoplasticError) -> Self {
        Self::Native(error)
    }
}
impl core::fmt::Display for HexThermoplasticHardnessError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled Hex8 hardness: {self:?}")
    }
}
impl std::error::Error for HexThermoplasticHardnessError {}

// Reuse the same identity, reference volume/mass, original caloric-object and
// exact accepted-enthalpy admission as the existing conservative heat driver.
fn checked_enthalpies<'a, F: ThermoplasticSourceFamily<'a>>(
    solid: &ThermoplasticHexSolid<'a, F>,
    thermal: &DepositionSimulation<'_>,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, HexThermoplasticHardnessError> {
    poll::<Infallible>(cx)?;
    if thermal.time_s() != solid.time_s {
        return Err(HexThermoplasticHardnessError::Invalid(
            "hardness requires matching accepted solid and caloric clocks",
        ));
    }
    let calorics: Vec<_> = (0..solid.cells.len())
        .map(|point| PointCaloric(solid.families[point / 8].thermal_source().caloric_source()))
        .collect();
    Ok(ThermalBinding {
        cells: &solid.cells,
        materials: FinitePlasticLawField::PerElement(&calorics),
    }
    .enthalpies::<Infallible>(thermal, true, cx)?)
}

fn resolved_point<'a, F: ThermoplasticSourceFamily<'a>>(
    solid: &ThermoplasticHexSolid<'a, F>,
    point: usize,
    h: f64,
) -> Result<(EquilibriumPhaseState, ResolvedJ2Plasticity), HexThermoplasticHardnessError> {
    Ok(solid.families[point / 8]
        .thermal_source()
        .resolve_enthalpy(h)
        .map_err(|error| {
            source_error::<Infallible>(point, ThermalSectionSourceError::Plastic(error.into()))
        })?)
}

impl<'a, F: ThermoplasticSourceFamily<'a>> ThermoplasticHexSolid<'a, F>
where
    <<F as ThermoplasticSourceFamily<'a>>::Material as SmallStrainLaw>::State:
        AccumulatedPlasticHistory,
{
    /// Pin an exact indentation test per original element and admit every point.
    /// Existing Hex8 budgets bound the entire field. The private native owner has
    /// already admitted the complete J2/Voce/Chaboche history; its source family
    /// and caloric inventory remain authoritative throughout continuation.
    pub fn bind_hardness_field(
        &self,
        thermal: &DepositionSimulation<'_>,
        selections: &[(PropertyKey, MaterialPropertySelection)],
        cx: &Cx<'_>,
    ) -> Result<HexThermoplasticHardnessBinding<'a, F>, HexThermoplasticHardnessError> {
        let h = checked_enthalpies(self, thermal, cx)?;
        if selections.len() != self.families.len() {
            return Err(HexThermoplasticHardnessError::Invalid(
                "one exact hardness selection per original Hex8 element required",
            ));
        }
        let mut sources = Vec::with_capacity(selections.len());
        for (element, (key, selection)) in selections.iter().enumerate() {
            poll::<Infallible>(cx)?;
            let point = 8 * element;
            let original = self.families[element].thermal_source();
            let (caloric, material) = resolved_point(self, point, h[point])?;
            let source = ThermoplasticHardnessSource::<'a, F>::try_from_native_state(
                original,
                &material,
                caloric,
                &self.solution.quadrature_states[point],
                key.clone(),
                selection.clone(),
            )
            .map_err(|error| HexThermoplasticHardnessError::Source {
                point,
                cell: self.cells[point].id,
                error,
            })?;
            // The remaining points must support this same exact selected claim.
            for point in (point + 1)..(point + 8) {
                poll::<Infallible>(cx)?;
                let (caloric, material) = resolved_point(self, point, h[point])?;
                source
                    .at_native_state(
                        original,
                        &material,
                        caloric,
                        &self.solution.quadrature_states[point],
                    )
                    .map_err(|error| HexThermoplasticHardnessError::Source {
                        point,
                        cell: self.cells[point].id,
                        error,
                    })?;
            }
            sources.push(source);
        }
        poll::<Infallible>(cx)?;
        Ok(HexThermoplasticHardnessBinding {
            nodes_m: self.nodes_m,
            hexahedra: self.hexahedra,
            families: Arc::clone(&self.families),
            formulation: self.discretization.formulation,
            cells: self.cells.iter().map(|cell| cell.id).collect(),
            sources,
        })
    }
}

impl<'a, F: ThermoplasticSourceFamily<'a>> HexThermoplasticHardnessBinding<'a, F>
where
    <<F as ThermoplasticSourceFamily<'a>>::Material as SmallStrainLaw>::State:
        AccumulatedPlasticHistory,
{
    /// Observe the original accepted point histories and exact caloric inventory.
    /// No native state, work, heat, velocity or observer is copied or advanced.
    pub fn at_accepted(
        &self,
        solid: &ThermoplasticHexSolid<'a, F>,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> Result<HexThermoplasticHardnessField, HexThermoplasticHardnessError> {
        poll::<Infallible>(cx)?;
        if !Arc::ptr_eq(&self.families, &solid.families)
            || !core::ptr::eq(self.nodes_m, solid.nodes_m)
            || !core::ptr::eq(self.hexahedra, solid.hexahedra)
            || self.formulation != solid.discretization.formulation
            || self.cells.len() != solid.cells.len()
            || self
                .cells
                .iter()
                .zip(&solid.cells)
                .any(|(id, cell)| *id != cell.id)
        {
            return Err(HexThermoplasticHardnessError::Invalid(
                "hardness requires the original mesh, ordered families, cells and formulation",
            ));
        }
        let h = checked_enthalpies(solid, thermal, cx)?;
        let mut points = Vec::with_capacity(self.cells.len());
        for (point, &h) in h.iter().enumerate() {
            poll::<Infallible>(cx)?;
            let (caloric, material) = resolved_point(solid, point, h)?;
            let hardness = self.sources[point / 8]
                .at_native_state(
                    solid.families[point / 8].thermal_source(),
                    &material,
                    caloric,
                    &solid.solution.quadrature_states[point],
                )
                .map_err(|error| HexThermoplasticHardnessError::Source {
                    point,
                    cell: self.cells[point],
                    error,
                })?;
            points.push(HexThermoplasticPointHardness {
                element: point / 8,
                local_point: point % 8,
                cell: self.cells[point],
                hardness,
            });
        }
        poll::<Infallible>(cx)?;
        Ok(HexThermoplasticHardnessField {
            time_s: solid.time_s,
            points,
        })
    }

    /// Consume an already converged original proposal and publish it only when
    /// every hardness endpoint is admitted. Late refusal leaves all accepted
    /// input owners unchanged and returns no advanced thermal/mechanical state.
    pub fn check_candidate<'t, O>(
        &self,
        candidate: HexThermoplasticCandidate<'a, 't, F, O>,
        cx: &Cx<'_>,
    ) -> Result<
        HexThermoplasticHardnessAccepted<HexThermoplasticCandidate<'a, 't, F, O>>,
        HexThermoplasticHardnessError,
    > {
        let hardness = self.at_accepted(candidate.solid(), candidate.thermal(), cx)?;
        Ok(HexThermoplasticHardnessAccepted {
            mechanics: candidate,
            hardness,
        })
    }
}
