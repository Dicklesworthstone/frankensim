//! Source-defined free geometry on the original finite circuit/enthalpy solve.
//!
//! Resistor and explicitly bound surface geometry are lowered here. Material
//! queries, caloric history, Kirchhoff evolution, field motion, Joule deposition
//! and boundary heat keep their owners.
//! The thermal cells retain their reference geometry and delivered mass.

use super::*;
use crate::circuit_enthalpy::motion::{FieldMotionEnthalpyCandidate, FieldMotionEnthalpyInput};
use crate::material_circuit::motion::MaterialFieldMotion;
use fs_material::conductor::{ResolvedConductorState, graph::MaterialOhmicNode};
use fs_material::thermal::transport::expansion::PreparedFreeExpansion;
use fs_thermal::enthalpy::CellGeometry;

mod boundary;
pub use boundary::{ExpandingBoundary, ExpandingBoundarySample, ExpandingNetworkSample};

/// One original resistor's measured engineering-strain path and current axis.
/// The fixed principal frame belongs to `source`; no geometry is inferred from
/// a material name, density change or thermal-expansion coefficient surrogate.
#[derive(Clone, Debug)]
pub struct ConductorFreeExpansion {
    pub branch: usize,
    pub source: PreparedFreeExpansion,
    pub longitudinal_axis: usize,
}

struct ExpandingBranch {
    cell: CellId,
    mass_kg: f64,
    geometry: CellGeometry,
    law: MaterialOhmicNode,
    boundary: Option<ExpandingBoundary>,
}

/// Prepared transport with sourced, unconstrained solid resistor expansion.
///
/// The SAME `MaterialOhmicNode` owns source admission and `L(T)/A(T)` lowering
/// for this native circuit and the constitutive graph. Its prepared geometry
/// reads no material database during a nonlinear trial. The original network
/// resolves every electrical/thermal property once at that trial's endpoint;
/// its complete bundle then moves through the common geometry law.
///
/// This is fixed-mass, quasistatic, stress-free expansion, without phase change
/// or mechanical loading. Robin/reservoir exchange requires explicit source-
/// bound surface partitions; occupied interbody contacts still need their own
/// geometric heat-transfer owner. Independent cell watts remain admissible.
/// Other cells retain the network's original heat-transfer laws. No mechanical
/// work, density-based mass update or geometry change to an inductor/capacitor
/// is inferred.
pub struct ExpandingMaterialNetwork<'n, 'a> {
    prepared: PreparedMaterialNetwork<'n, 'a>,
    branches: BTreeMap<usize, ExpandingBranch>,
}

impl<'n, 'a> PreparedMaterialNetwork<'n, 'a> {
    /// Attach one or more principal expansion paths to original whole-cell
    /// resistors. Their original lengths/areas describe the supplied reference
    /// state; every path must have that exact card, query and temperature axis.
    /// The complete caloric chart must remain solid. Electrical and geometric
    /// claims stay pinned across subsequent steps and refusals.
    pub fn with_free_expansion(
        self,
        thermal: &DepositionSimulation<'_>,
        expansions: &[ConductorFreeExpansion],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ExpandingMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        self.with_free_expansion_and_boundaries(thermal, expansions, &[], cancelled)
    }

    /// The same original free expansion with explicit disjoint ambient and
    /// reservoir surface patches. Their principal directions and strain path
    /// must match the electrical geometry; their reference areas must cover the
    /// actual original thermal surface exactly. Robin h remains per current
    /// area. Reservoir G is explicitly interpreted as a uniform areal law at
    /// its declared reference coverage; native boundary events may change h,
    /// G and reservoir temperatures while preserving this geometric partition.
    pub fn with_free_expansion_and_boundaries(
        self,
        thermal: &DepositionSimulation<'_>,
        expansions: &[ConductorFreeExpansion],
        boundaries: &[ExpandingBoundary],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExpandingMaterialNetwork<'n, 'a>, MaterialNetworkError> {
        let invalid = MaterialNetworkError::Invalid;
        if cancelled() {
            return Err(DepositionError::Cancelled.into());
        }
        if expansions.is_empty() || expansions.len() > self.network.resistors.len() {
            return Err(invalid(
                "bounded nonempty original resistor expansion bindings required",
            ));
        }
        if boundaries.len() > expansions.len()
            || boundaries.iter().enumerate().any(|(i, boundary)| {
                !expansions
                    .iter()
                    .any(|expansion| expansion.branch == boundary.branch)
                    || boundaries[..i]
                        .iter()
                        .any(|other| other.branch == boundary.branch)
            })
        {
            return Err(invalid(
                "boundary surfaces need unique expanded resistor branches",
            ));
        }
        let original = self.resolve(thermal, &mut cancelled)?;
        let mut branches = BTreeMap::new();
        for expansion in expansions {
            if cancelled() {
                return Err(DepositionError::Cancelled.into());
            }
            let source = self
                .network
                .resistors
                .iter()
                .find(|r| r.destination.branch == expansion.branch)
                .ok_or(invalid("expansion names no original resistor branch"))?;
            let binding = self
                .network
                .cells
                .iter()
                .find(|c| c.cell == source.destination.cell)
                .ok_or(invalid("expanding resistor has no caloric binding"))?;
            let cell = &thermal.cells()[&binding.cell];
            let phase = binding
                .curve
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(MaterialNetworkError::Phase)?;
            let reference_mass = phase.bulk_density_kg_m3() * cell.geometry.volume_m3;
            if binding.point.liquid_fraction_axis.is_some()
                || binding
                    .curve
                    .knots()
                    .iter()
                    .any(|k| k.liquid_mass_fraction != 0.0)
                || !reference_mass.is_finite()
                || (reference_mass / cell.mass_kg - 1.0).abs() > 64.0 * f64::EPSILON
            {
                return Err(invalid(
                    "free resistor expansion requires a wholly solid, matching fixed reference mass",
                ));
            }
            let conductor = &original
                .resistors
                .iter()
                .find(|r| r.branch == expansion.branch)
                .ok_or(invalid("expanding resistor has no original conductor"))?
                .conductor;
            let state = ResolvedConductorState::try_from_material_with_property(
                conductor.material().clone(),
                conductor.source_basis(),
                conductor.source_key(),
            )
            .map_err(|e| MaterialNetworkError::Conductor(e.into()))?;
            let axis = binding.point.temperature_axis;
            let endpoints = expansion.source.query_endpoints();
            let range = [endpoints[0], endpoints[1]].map(|p| p.axes().get(axis).copied());
            let [Some(low), Some(high)] = range else {
                return Err(invalid(
                    "free expansion must use the caloric temperature axis",
                ));
            };
            let knots = binding.curve.knots();
            if low < knots[0].temperature_k || high > knots[knots.len() - 1].temperature_k {
                return Err(invalid(
                    "free expansion interval exceeds the original solid caloric chart",
                ));
            }
            let law = MaterialOhmicNode::try_new(
                binding.point.material,
                state,
                axis,
                [low, high],
                source.length_m,
                source.area_m2,
            )
            .and_then(|law| {
                law.with_free_expansion(expansion.source.clone(), expansion.longitudinal_axis)
            })
            .map_err(|e| MaterialNetworkError::Conductor(e.into()))?;
            let entry = ExpandingBranch {
                cell: binding.cell,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
                law,
                boundary: boundaries
                    .iter()
                    .find(|b| b.branch == expansion.branch)
                    .cloned(),
            };
            entry.admit(thermal)?;
            if branches.insert(expansion.branch, entry).is_some() {
                return Err(invalid("duplicate original resistor expansion binding"));
            }
        }
        if cancelled() {
            return Err(DepositionError::Cancelled.into());
        }
        Ok(ExpandingMaterialNetwork {
            prepared: self,
            branches,
        })
    }
}

impl ExpandingBranch {
    fn admit(&self, thermal: &DepositionSimulation<'_>) -> Result<(), MaterialNetworkError> {
        let invalid = MaterialNetworkError::Invalid;
        let cell = thermal
            .cells()
            .get(&self.cell)
            .ok_or(invalid("expanding resistor thermal cell is absent"))?;
        if cell.mass_kg.to_bits() != self.mass_kg.to_bits() || cell.geometry != self.geometry {
            return Err(invalid(
                "free resistor expansion must retain delivered mass and reference thermal geometry",
            ));
        }
        if thermal.links().values().any(|link| {
            if link.a != self.cell && link.b != self.cell {
                return false;
            }
            let side = usize::from(link.b == self.cell);
            thermal
                .contact_patch(link.id)
                .is_some_and(|patch| patch.covered_area_m2[side] != 0.0)
                || link.area_m2 != 0.0
                || link.conductance_w_k != 0.0
        }) {
            return Err(invalid(
                "expanding resistor contact geometry needs its own heat-transfer owner",
            ));
        }
        if let Some(boundary) = &self.boundary {
            boundary.admit(self, thermal)?;
        } else if cell.boundary.htc_w_m2_k != 0.0 || cell.boundary.reservoir.is_some() {
            return Err(invalid(
                "expanding resistor boundary needs its explicit source surface partition",
            ));
        }
        Ok(())
    }
}

impl ExpandingMaterialNetwork<'_, '_> {
    /// The original source/geometry law for an expanded native branch.
    pub fn expansion(&self, branch: usize) -> Option<&MaterialOhmicNode> {
        self.branches.get(&branch).map(|entry| &entry.law)
    }

    /// Same original heat destinations, retaining reference cell identity.
    pub fn resistor_destinations(&self) -> Vec<ResistorCell> {
        self.prepared.resistor_destinations()
    }

    /// Resolve original transport ONCE and lower its exact endpoint bundle
    /// through prepared source geometry. No alternate R(T) or heat law exists.
    pub fn resolve(
        &self,
        thermal: &DepositionSimulation<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CircuitThermalCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        if self.branches.values().any(|entry| entry.boundary.is_some()) {
            return Err(MaterialNetworkError::Invalid(
                "surface-bound expansion requires the complete boundary-aware resolver",
            ));
        }
        self.resolve_bulk(thermal, cancelled)
    }

    fn resolve_bulk(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CircuitThermalCoefficients<MaterialNetworkSample>, MaterialNetworkError> {
        let mut result = self.prepared.resolve(thermal, &mut cancelled)?;
        for resistor in &mut result.resistors {
            if cancelled() {
                return Err(DepositionError::Cancelled.into());
            }
            if let Some(entry) = self.branches.get(&resistor.branch) {
                entry.admit(thermal)?;
                let state = ResolvedConductorState::try_from_material_with_property(
                    resistor.conductor.material().clone(),
                    resistor.conductor.source_basis(),
                    resistor.conductor.source_key(),
                )
                .map_err(|e| MaterialNetworkError::Conductor(e.into()))?;
                resistor.conductor = entry
                    .law
                    .resolve_from_state(state)
                    .map_err(|e| MaterialNetworkError::Conductor(e.into()))?;
            }
        }
        if cancelled() {
            return Err(DepositionError::Cancelled.into());
        }
        Ok(result)
    }

    fn admit_circuit(
        &self,
        thermal: &DepositionSimulation<'_>,
        circuit: &MaterialCircuit,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), MaterialNetworkError> {
        let current = self
            .resolve_with_boundary_conductances(thermal, cancelled)?
            .bulk;
        if current.resistors.len() != circuit.resistors().len()
            || current.resistors.iter().any(|source| {
                !circuit
                    .resistors()
                    .iter()
                    .any(|accepted| accepted == source)
            })
        {
            return Err(MaterialNetworkError::Invalid(
                "accepted circuit must retain the actual thermal source and expanded geometry",
            ));
        }
        Ok(())
    }

    /// Continue the original charge/flux and finite-enthalpy histories together.
    /// Source-driven geometry changes native current and its ONCE-deposited
    /// Joule heat in the same endpoint residual. This is the existing first-order
    /// thermal coupling, not a new solver or an accuracy/passivity certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        CircuitEnthalpyCandidate<'m, ExpandingNetworkSample>,
        CircuitEnthalpyError<MaterialNetworkError>,
    > {
        let cancel = RefCell::new(cancelled);
        self.admit_circuit(thermal, circuit, || (*cancel.borrow_mut())())
            .map_err(coupled_error)?;
        let destinations = self.resistor_destinations();
        map_cancel(
            crate::circuit_enthalpy::step_candidate_with_boundary_conductances(
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

    /// The same prepared resistor geometry drives the existing reciprocal
    /// electric/magnetic gap actuator. The resistor bodies expand freely;
    /// actuator geometry, force and mechanical work retain their own owner.
    #[allow(clippy::too_many_arguments)]
    pub fn field_motion_candidate<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        motion: &MaterialFieldMotion,
        ports: &[f64],
        external_force_n: f64,
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        FieldMotionEnthalpyCandidate<'m, ExpandingNetworkSample>,
        CircuitEnthalpyError<MaterialNetworkError>,
    > {
        let cancel = RefCell::new(cancelled);
        self.admit_circuit(thermal, motion.circuit(), || (*cancel.borrow_mut())())
            .map_err(coupled_error)?;
        let destinations = self.resistor_destinations();
        map_cancel(
            crate::circuit_enthalpy::motion::step_candidate_with_boundary_conductances(
                thermal,
                target_s,
                FieldMotionEnthalpyInput {
                    motion,
                    ports,
                    external_force_n,
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

fn coupled_error(error: MaterialNetworkError) -> CircuitEnthalpyError<MaterialNetworkError> {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(other),
    }
}

fn map_cancel<T>(
    result: Result<T, CircuitEnthalpyError<MaterialNetworkError>>,
) -> Result<T, CircuitEnthalpyError<MaterialNetworkError>> {
    result.map_err(|error| match error {
        CircuitEnthalpyError::Source(source) => coupled_error(source),
        other => other,
    })
}
