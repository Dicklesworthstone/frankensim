//! The same source geometry supplies electrical section and thermal surface.
use super::*;
use crate::circuit_enthalpy::{BoundaryCircuitThermalCoefficients, entropy::CircuitEntropyAudit};
use fs_material::thermal::transport::expansion::surface::{
    ExpansionSurfaceSample, PreparedExpansionSurface,
};
use fs_thermal::enthalpy::backward_euler::BoundaryConductances;

/// Disjoint reference surface partitions on one original expanding resistor.
/// Face-normal areas are explicit geometry; electrical length/section cannot
/// determine them. Reservoir coverage must match the native boundary exactly.
#[derive(Clone, Debug)]
pub struct ExpandingBoundary {
    pub branch: usize,
    pub exposed: PreparedExpansionSurface,
    pub reservoir: Option<PreparedExpansionSurface>,
}

/// Actual source areas and signed native boundary heat at the thermal endpoint.
/// Watts are positive into the cell; they are OBSERVATIONS of the replacement
/// conductances, never another heat source deposited alongside that boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandingBoundarySample {
    pub branch: usize,
    pub cell: CellId,
    pub temperature_k: f64,
    pub exposed: ExpansionSurfaceSample,
    pub reservoir: Option<ExpansionSurfaceSample>,
    pub ambient_conductance_w_k: f64,
    pub reservoir_conductance_w_k: f64,
    pub ambient_heat_into_cell_w: f64,
    pub reservoir_heat_into_cell_w: f64,
}

/// Complete original transport and surface geometry at the accepted endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandingNetworkSample {
    pub material: MaterialNetworkSample,
    /// Stable branch ordering. Entries exist only for explicitly bound surfaces.
    pub boundaries: Vec<ExpandingBoundarySample>,
}

fn same_area(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a >= 0.0
        && b >= 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.abs().max(b.abs())
}
fn product(a: f64, b: f64) -> Result<f64, MaterialNetworkError> {
    let value = a * b;
    if !value.is_finite() || (a != 0.0 && b != 0.0 && value == 0.0) {
        return Err(MaterialNetworkError::Invalid(
            "expanded boundary conductance or signed heat is unrepresentable",
        ));
    }
    Ok(value)
}

impl ExpandingBoundary {
    pub(super) fn admit(
        &self,
        entry: &ExpandingBranch,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), MaterialNetworkError> {
        let invalid = MaterialNetworkError::Invalid;
        let original = entry
            .law
            .free_expansion()
            .expect("prepared expanding branch")
            .0;
        if self.exposed.source().identity() != original.identity()
            || self
                .reservoir
                .as_ref()
                .is_some_and(|s| s.source().identity() != original.identity())
        {
            return Err(invalid(
                "electrical and thermal surfaces must retain the identical expansion source and frame",
            ));
        }
        let cell = &thermal.cells()[&entry.cell];
        let reservoir_area = match (&self.reservoir, cell.boundary.reservoir) {
            (None, None) => 0.0,
            (Some(surface), Some(boundary)) => {
                let area = surface.reference_total_area_m2();
                if !same_area(area, boundary.area_m2)
                    || (area == 0.0 && boundary.conductance_w_k != 0.0)
                {
                    return Err(invalid(
                        "source reservoir faces must retain the native reference coverage and finite areal law",
                    ));
                }
                area
            }
            _ => {
                return Err(invalid(
                    "reservoir surface partition and native boundary presence disagree",
                ));
            }
        };
        let exposed = self.exposed.reference_total_area_m2();
        if !same_area(exposed + reservoir_area, entry.geometry.surface_area_m2)
            || !same_area(exposed, thermal.exposed_area_m2(entry.cell)?)
        {
            return Err(invalid(
                "disjoint source faces must exactly partition the original exposed and reservoir surface",
            ));
        }
        Ok(())
    }

    fn sample(
        &self,
        entry: &ExpandingBranch,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ExpandingBoundarySample, MaterialNetworkError> {
        // resolve_bulk has admitted this exact immutable state and partition.
        let cell = &thermal.cells()[&entry.cell];
        let temperature_k = thermal.temperature_k(entry.cell)?;
        let geometry = entry
            .law
            .free_expansion()
            .expect("prepared expanding branch")
            .0
            .at_temperature_with_tangent(temperature_k)
            .map_err(MaterialNetworkError::Expansion)?;
        let exposed = self
            .exposed
            .from_geometry(geometry)
            .map_err(MaterialNetworkError::Expansion)?;
        let ambient_conductance_w_k = product(cell.boundary.htc_w_m2_k, exposed.total_area_m2())?;
        let (reservoir, reservoir_conductance_w_k, reservoir_heat_into_cell_w) =
            match (&self.reservoir, cell.boundary.reservoir) {
                (Some(surface), Some(boundary)) => {
                    let sample = surface
                        .from_geometry(geometry)
                        .map_err(MaterialNetworkError::Expansion)?;
                    let reference = surface.reference_total_area_m2();
                    let ratio = if reference == 0.0 {
                        1.0
                    } else {
                        sample.total_area_m2() / reference
                    };
                    if !ratio.is_finite() || ratio <= 0.0 {
                        return Err(MaterialNetworkError::Invalid(
                            "expanded reservoir area ratio is unrepresentable",
                        ));
                    }
                    let g = product(boundary.conductance_w_k, ratio)?;
                    (
                        Some(sample),
                        g,
                        product(g, boundary.temperature_k - temperature_k)?,
                    )
                }
                _ => (None, 0.0, 0.0),
            };
        Ok(ExpandingBoundarySample {
            branch: self.branch,
            cell: entry.cell,
            temperature_k,
            exposed,
            reservoir,
            ambient_conductance_w_k,
            reservoir_conductance_w_k,
            ambient_heat_into_cell_w: product(
                ambient_conductance_w_k,
                cell.boundary.ambient_k - temperature_k,
            )?,
            reservoir_heat_into_cell_w,
        })
    }
}

impl ExpandingMaterialNetwork<'_, '_> {
    /// Complete native endpoint coefficients, including current surface areas.
    /// Reference h and total reservoir G retain their authored event semantics;
    /// the explicit surface partition supplies current exposed/covered areas.
    /// G_res=G_ref*A_res(T)/A_res,ref is an explicitly uniform areal law, not
    /// changing-contact mechanics or a bulk conduction-path approximation.
    /// No boundary heat is inserted into the additional-watts map.
    pub fn resolve_with_boundary_conductances(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BoundaryCircuitThermalCoefficients<ExpandingNetworkSample>, MaterialNetworkError>
    {
        let original = self.resolve_bulk(thermal, &mut cancelled)?;
        let mut boundaries = BoundaryConductances::default();
        let mut samples = Vec::new();
        for entry in self.branches.values() {
            if cancelled() {
                return Err(DepositionError::Cancelled.into());
            }
            if let Some(boundary) = &entry.boundary {
                let sample = boundary.sample(entry, thermal)?;
                boundaries
                    .ambient_w_k
                    .insert(entry.cell, sample.ambient_conductance_w_k);
                if sample.reservoir.is_some() {
                    boundaries
                        .reservoir_w_k
                        .insert(entry.cell, sample.reservoir_conductance_w_k);
                }
                samples.push(sample);
            }
        }
        if cancelled() {
            return Err(DepositionError::Cancelled.into());
        }
        Ok(BoundaryCircuitThermalCoefficients {
            bulk: CircuitThermalCoefficients {
                resistors: original.resistors,
                sources_w: original.sources_w,
                conductances_w_k: original.conductances_w_k,
                source: ExpandingNetworkSample {
                    material: original.source,
                    boundaries: samples,
                },
            },
            boundaries,
        })
    }

    /// Same native current and finite heat step, with the original reservoir
    /// entropy account evaluated only during final acceptance. Expanded area
    /// changes the actual reservoir heat/entropy; it is never unclassified work.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate_with_entropy<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        (
            CircuitEnthalpyCandidate<'m, ExpandingNetworkSample>,
            CircuitEntropyAudit,
        ),
        CircuitEnthalpyError<MaterialNetworkError>,
    > {
        let cancel = RefCell::new(cancelled);
        self.admit_circuit(thermal, circuit, || (*cancel.borrow_mut())())
            .map_err(coupled_error)?;
        let destinations = self.resistor_destinations();
        map_cancel(
            crate::circuit_enthalpy::entropy::step_candidate_with_boundary_conductances(
                thermal,
                target_s,
                CircuitEnthalpyInput {
                    circuit,
                    ports,
                    destinations: &destinations,
                },
                controls,
                guess_j,
                |state| self.resolve_with_boundary_conductances(state, || (*cancel.borrow_mut())()),
                || (*cancel.borrow_mut())(),
            ),
        )
    }
}
