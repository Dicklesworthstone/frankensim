//! Spatial winding Joule heat on the original finite enthalpy/Newton owner.
//!
//! The magnetic basis and native Kirchhoff history remain with fs-em/fs-phs.
//! This adapter reuses the existing MaterialCell and prepared thermal-transport
//! bindings: caloric storage, conductor resistivity and thermal conductivity
//! must come from the SAME cards and actual endpoint temperatures. Mapped
//! magnetic triangles cover whole fixed thermal volumes; no guessed lumped resistor,
//! fill fraction, copper coefficient or alternate heat integrator is introduced.

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ResistorCell, ThermalDriver, cell_driver_candidate, finite,
    material::{
        ConductionLink, MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
    },
};
use fs_em::winding::ohmic::{
    OhmicWindingCandidate, OhmicWindingCircuit, OhmicWindingControls, OhmicWindingError,
    OhmicWindingStep, WindingConductor,
};
use fs_material::{
    conductor::{ConductorError, ResolvedConductorState},
    thermal::source::ThermalMaterialError,
};
use fs_phs::circuit::Branch;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// Complete correspondence for one ORIGINAL magnetic triangle. Conducting
/// triangles must use a thermal Cell; their loss is deposited there. Fixed
/// explicitly holds an unheated magnetic region at its prepared material state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MagneticRegion {
    /// Existing thermal control volume owns this region's temperature.
    Cell { triangle: usize, cell: CellId },
    /// Explicit isothermal exterior/core region outside the thermal model.
    Fixed { triangle: usize },
}
impl MagneticRegion {
    fn triangle(self) -> usize {
        match self {
            Self::Cell { triangle, .. } | Self::Fixed { triangle } => triangle,
        }
    }
}

/// The complete exact magnetic source used/admitted at the thermal endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticThermalSample {
    pub triangle: usize,
    /// None retains the explicitly fixed prepared source; no heat is sent there.
    pub cell: Option<CellId>,
    pub material: fs_material::state_point::ResolvedMaterialStatePoint,
}

/// Original caloric/conductive samples and every magnetic region's source.
#[derive(Clone, Debug, PartialEq)]
pub struct WindingThermalSample {
    pub network: MaterialNetworkSample,
    pub magnetic: Vec<MagneticThermalSample>,
}

/// Data-only material, geometry and complete heat routing for stranded windings.
/// Thermal cells outside the magnetic geometry keep their original latent and
/// boundary behavior. Mapped magnetic cells admit fixed solid geometry/density;
/// thermal expansion, changing permeability/inductance and moving conductors
/// require their own coupled work/geometry owner.
pub struct WindingThermalNetwork<'a> {
    /// Same caloric/material bindings as the existing thermal network, complete
    /// for all active cells including passive graph-resistor heat destinations.
    pub cells: &'a [MaterialCell<'a>],
    /// Complete existing thermal-link material/conduction bindings.
    pub links: &'a [ConductionLink],
    /// Exactly one explicit thermal or fixed-state correspondence per magnetic
    /// triangle, including nonconducting cores and surrounding field regions.
    pub magnetic_regions: &'a [MagneticRegion],
    /// Exactly one destination for EVERY original numeric graph resistor.
    /// These retain their authored resistance; only intrinsic winding loss
    /// derives from the represented source material and current profile.
    pub graph_resistors: &'a [ResistorCell],
    /// Absolute agreement allowance [m3] for covered cell volume and mass/rho.
    /// This is geometry admission, not a license to infer unresolved fill.
    pub volume_tolerance_m3: f64,
}

/// Source and topology refusals remain distinct from nonlinear/energy failures.
#[derive(Debug)]
pub enum WindingThermalError {
    Material(MaterialNetworkError),
    Conductor(ConductorError),
    Winding(OhmicWindingError),
    Magnetic(fs_material::field_storage::FieldStorageError),
    Invalid(&'static str),
}
impl core::fmt::Display for WindingThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material winding/enthalpy: {self:?}")
    }
}
impl std::error::Error for WindingThermalError {}
type Error = CircuitEnthalpyError<WindingThermalError>;
fn invalid(what: &'static str) -> Error {
    CircuitEnthalpyError::Source(WindingThermalError::Invalid(what))
}
fn material(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(WindingThermalError::Material(other)),
    }
}
fn electrical(error: OhmicWindingError) -> Error {
    match error {
        OhmicWindingError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(WindingThermalError::Winding(other)),
    }
}
fn check(cancelled: &RefCell<impl FnMut() -> bool>) -> Result<(), Error> {
    if (*cancelled.borrow_mut())() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}

/// Both accepted candidates belong to the same interval and include all heat.
/// Original inputs remain unchanged through every Newton/Jacobian trial.
pub struct WindingEnthalpyCandidate<'t, 'c, 'm> {
    thermal: DepositionSimulation<'t>,
    circuit: OhmicWindingCircuit<'c, 'm>,
    report: CircuitEnthalpyReport<WindingThermalSample, OhmicWindingStep>,
}
impl<'t, 'c, 'm> WindingEnthalpyCandidate<'t, 'c, 'm> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn circuit(&self) -> &OhmicWindingCircuit<'c, 'm> {
        &self.circuit
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<WindingThermalSample, OhmicWindingStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        OhmicWindingCircuit<'c, 'm>,
        CircuitEnthalpyReport<WindingThermalSample, OhmicWindingStep>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}
impl<'c, 'm> ThermalDriver for OhmicWindingCandidate<'c, 'm> {
    type Model = OhmicWindingCircuit<'c, 'm>;
    type Step = OhmicWindingStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            stored_change_j: step.stored_after_j - step.stored_before_j,
            supplied_j: step.native.supplied,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        OhmicWindingCandidate::into_parts(self)
    }
}

struct Routes {
    magnetic: BTreeMap<usize, Option<CellId>>,
    windings: BTreeMap<usize, CellId>,
    graph: BTreeMap<usize, CellId>,
}
impl WindingThermalNetwork<'_> {
    /// Solve temperature-dependent winding current and ALL circuit Joule heat
    /// on the SAME native backward-Euler cell residual. Coefficients use each
    /// actual proposed endpoint; native electrical currents use its midpoint.
    /// This is generally first-order coupled accuracy, not a stability claim.
    ///
    /// Existing caloric law instances/card identities, fixed coordinates,
    /// complete link sources and accepted electrical source choices remain pinned.
    /// Thermal transport choices are pinned across this interval's trials.
    /// The initial electrical snapshots must match the accepted thermal state
    /// and clock. Only same-card temperature refresh is permitted by fs-em;
    /// changed mu or temperature-dependent magnetic storage refuses, including
    /// nonconducting cores. Every other magnetic triangle is explicitly Fixed
    /// at its prepared material state; omission never implies isothermal support.
    ///
    /// Thermally mapped magnetic cells require a wholly solid constant-density caloric chart
    /// on this fixed-geometry lane. Numeric graph-resistor heat may feed any
    /// explicitly bound active cell. Process start events first and split at
    /// interior source/topology events; target events remain pending.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, 'c, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &OhmicWindingCircuit<'c, 'm>,
        ports: &[f64],
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<WindingEnthalpyCandidate<'t, 'c, 'm>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        if circuit.time_s().to_bits() != thermal.time_s().to_bits()
            || !self.volume_tolerance_m3.is_finite()
            || self.volume_tolerance_m3 < 0.0
            || self.cells.len() != thermal.cells().len()
            || self.cells.len() > controls.maximum_cells
        {
            return Err(invalid(
                "matching electrical/thermal clocks, complete bounded cells and finite volume tolerance required",
            ));
        }
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let routes = self.admit_routes(thermal, circuit, &cancel)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let (initial, _) =
            self.materials_at(thermal, circuit, &accepted.source, &routes, &cancel)?;
        if initial.len() != circuit.conductors().len()
            || initial
                .iter()
                .zip(circuit.conductors())
                .any(|(a, b)| !same_conductor_state(a, b))
        {
            return Err(invalid(
                "accepted winding source snapshots must match the actual thermal endpoint",
            ));
        }
        let ohmic = OhmicWindingControls {
            power_tolerance_w: controls.circuit.power_tolerance_w,
            energy_tolerance_j: controls.circuit.energy_tolerance_j,
        };
        let (thermal, circuit, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let (sources, magnetic) =
                    self.materials_at(state, circuit, &coefficients.source, &routes, &cancel)?;
                let circuit = circuit
                    .step_to_candidate(target_s, sources, ports, ohmic, &mut || {
                        (*cancel.borrow_mut())()
                    })
                    .map_err(electrical)?;
                Ok(DrivenCoefficients {
                    circuit,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: WindingThermalSample {
                        network: coefficients.source,
                        magnetic,
                    },
                    circuit_trials: 1,
                })
            },
            |candidate, dt, checkpoint| routes.deposit(candidate, dt, checkpoint),
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(WindingEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }

    fn admit_routes(
        &self,
        thermal: &DepositionSimulation<'_>,
        circuit: &OhmicWindingCircuit<'_, '_>,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Routes, Error> {
        let problem = circuit.magnetic().basis().problem();
        if self.magnetic_regions.len() != problem.mesh().faces().len()
            || self.graph_resistors.len()
                != circuit
                    .magnetic()
                    .graph()
                    .branches
                    .iter()
                    .filter(|(_, _, b)| matches!(b, Branch::Resistor { .. }))
                    .count()
        {
            return Err(invalid(
                "every magnetic triangle needs an explicit thermal/fixed region, and every graph resistor needs one heat destination",
            ));
        }
        let mut magnetic = BTreeMap::new();
        let mut graph = BTreeMap::new();
        let mut volumes = BTreeMap::new();
        for &region in self.magnetic_regions {
            check(cancelled)?;
            let triangle = region.triangle();
            let cell = match region {
                MagneticRegion::Cell { cell, .. } => Some(cell),
                MagneticRegion::Fixed { .. } => None,
            };
            if triangle >= problem.mesh().faces().len()
                || cell.is_some_and(|id| !thermal.cells().contains_key(&id))
                || magnetic.insert(triangle, cell).is_some()
            {
                return Err(invalid(
                    "duplicate, foreign or inactive magnetic region correspondence",
                ));
            }
            if let Some(cell) = cell {
                let volume = problem
                    .mesh()
                    .face_measure(triangle)
                    .ok_or_else(|| invalid("missing original magnetic triangle geometry"))?;
                let sum = volumes.entry(cell).or_insert(0.0);
                *sum = finite(*sum + volume)?;
            }
        }
        let mut windings = BTreeMap::new();
        for source in circuit.conductors() {
            check(cancelled)?;
            let cell = magnetic
                .get(&source.triangle)
                .copied()
                .flatten()
                .ok_or_else(|| {
                    invalid(
                        "every conducting triangle must have an actual thermal heat destination",
                    )
                })?;
            windings.insert(source.triangle, cell);
        }
        for (id, volume) in volumes {
            check(cancelled)?;
            let cell = &thermal.cells()[&id];
            if (volume - cell.geometry.volume_m3).abs() > self.volume_tolerance_m3 {
                return Err(invalid(
                    "mapped magnetic triangle volumes must cover the whole receiving thermal cell",
                ));
            }
            let binding = self
                .cells
                .iter()
                .find(|s| s.cell == id)
                .ok_or_else(|| invalid("magnetic thermal cell lacks caloric binding"))?;
            for knot in binding.curve.knots() {
                check(cancelled)?;
                let density_volume = finite(cell.mass_kg / knot.bulk_density_kg_m3)?;
                if knot.liquid_mass_fraction != 0.0
                    || (density_volume - cell.geometry.volume_m3).abs() > self.volume_tolerance_m3
                {
                    return Err(invalid(
                        "fixed magnetic geometry requires a fully solid caloric chart with matching constant density",
                    ));
                }
            }
        }
        for route in self.graph_resistors {
            check(cancelled)?;
            if !thermal.cells().contains_key(&route.cell)
                || graph.insert(route.branch, route.cell).is_some()
                || !matches!(
                    circuit.magnetic().graph().branches.get(route.branch),
                    Some((_, _, Branch::Resistor { .. }))
                )
            {
                return Err(invalid(
                    "duplicate, foreign or inactive graph-resistor heat destination",
                ));
            }
        }
        check(cancelled)?;
        Ok(Routes {
            magnetic,
            windings,
            graph,
        })
    }

    fn materials_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        circuit: &OhmicWindingCircuit<'_, '_>,
        sample: &MaterialNetworkSample,
        routes: &Routes,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(Vec<WindingConductor>, Vec<MagneticThermalSample>), Error> {
        use fs_material::field_storage::{
            MAGNETIC_PERMEABILITY_PROPERTY, linear_magnetic_permeability,
        };
        let bindings = self
            .cells
            .iter()
            .map(|c| (c.cell, c))
            .collect::<BTreeMap<_, _>>();
        let phases = sample
            .cells
            .iter()
            .map(|c| (c.cell, c.phase))
            .collect::<BTreeMap<_, _>>();
        let conductors = circuit
            .conductors()
            .iter()
            .map(|c| (c.triangle, c))
            .collect::<BTreeMap<_, _>>();
        let problem = circuit.magnetic().basis().problem();
        let mut result = Vec::with_capacity(conductors.len());
        let mut magnetic = Vec::with_capacity(routes.magnetic.len());
        for (&triangle, &cell) in &routes.magnetic {
            check(cancelled)?;
            let original = &problem.materials()[triangle];
            let source = conductors.get(&triangle);
            let snapshot = if let Some(id) = cell {
                let binding = bindings
                    .get(&id)
                    .ok_or_else(|| invalid("magnetic cell lacks original material binding"))?;
                let phase = *phases
                    .get(&id)
                    .ok_or_else(|| invalid("magnetic cell lacks actual endpoint phase"))?;
                if !phase.is_fully_solid()
                    || (finite(thermal.cells()[&id].mass_kg / phase.bulk_density_kg_m3())?
                        - thermal.cells()[&id].geometry.volume_m3)
                        .abs()
                        > self.volume_tolerance_m3
                {
                    return Err(invalid("fixed solid magnetic density/geometry changed"));
                }
                let query = binding
                    .point
                    .at(phase)
                    .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                let base = source.map_or(original, |s| s.conductor.material());
                let snapshot = base
                    .at_query_point(binding.point.material, &query)
                    .map_err(|e| {
                        material(MaterialNetworkError::Material(
                            ThermalMaterialError::Material(e),
                        ))
                    })?;
                let mut changed_temperature = false;
                for ((name, value), (original_name, original_value)) in
                    snapshot.query_point().iter().zip(original.query_point())
                {
                    if name != original_name {
                        return Err(invalid("magnetic coordinate correspondence changed"));
                    }
                    if value.to_bits() != original_value.to_bits() {
                        if name != binding.point.temperature_axis {
                            return Err(invalid(
                                "only the owned magnetic-region temperature may change",
                            ));
                        }
                        changed_temperature = true;
                    }
                }
                let mu = linear_magnetic_permeability(&snapshot)
                    .map_err(|e| CircuitEnthalpyError::Source(WindingThermalError::Magnetic(e)))?;
                let property = snapshot
                    .property(MAGNETIC_PERMEABILITY_PROPERTY)
                    .expect("admitted permeability");
                binding
                    .point
                    .admit_static(property)
                    .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                if mu.to_bits() != problem.permeability_h_m()[triangle].to_bits()
                    || changed_temperature
                        && property.source_axes().iter().any(|(axis, _)| {
                            axis == binding.point.temperature_axis
                        })
                {
                    return Err(invalid(
                        "heated magnetic region changes prepared storage; new extraction and magnetic parameter work required",
                    ));
                }
                snapshot
            } else {
                original.clone()
            };
            if let Some(source) = source {
                let conductor = ResolvedConductorState::try_from_material_with_basis(
                    snapshot.clone(),
                    source.conductor.source_basis(),
                )
                .map_err(|e| CircuitEnthalpyError::Source(WindingThermalError::Conductor(e)))?;
                result.push(WindingConductor {
                    triangle,
                    conductor,
                });
            }
            magnetic.push(MagneticThermalSample {
                triangle,
                cell,
                material: snapshot,
            });
        }
        check(cancelled)?;
        Ok((result, magnetic))
    }
}

// Refresh pins the exact selected claims. Its honest usage receipt therefore
// changes policy even at an unchanged physical state; receipt identity is not
// an equality test for accepted temperature. Keep both full receipts intact.
fn same_conductor_state(a: &WindingConductor, b: &WindingConductor) -> bool {
    let x = a.conductor.material();
    let y = b.conductor.material();
    a.triangle == b.triangle
        && a.conductor.source_basis() == b.conductor.source_basis()
        && a.conductor.resistivity_ohm_m().to_bits() == b.conductor.resistivity_ohm_m().to_bits()
        && x.card_identity() == y.card_identity()
        && x.material() == y.material()
        && x.axis_quantities() == y.axis_quantities()
        && x.query_point().len() == y.query_point().len()
        && x.query_point().iter().zip(y.query_point()).all(
            |((name, value), (other_name, other_value))| {
                name == other_name && value.to_bits() == other_value.to_bits()
            },
        )
        && x.properties().len() == y.properties().len()
        && x.properties().iter().zip(y.properties()).all(|(p, q)| {
            p.requirement() == q.requirement()
                && p.source_axes() == q.source_axes()
                && p.answer().receipt.selected == q.answer().receipt.selected
                && p.value_si().to_bits() == q.value_si().to_bits()
        })
}

impl Routes {
    fn deposit(
        &self,
        candidate: &OhmicWindingCandidate<'_, '_>,
        dt: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Error> {
        let step = candidate.step();
        if step.cells.len() != self.windings.len() || step.graph_resistors.len() != self.graph.len()
        {
            return Err(invalid(
                "native winding/graph heat receipts changed correspondence",
            ));
        }
        let mut powers = BTreeMap::new();
        let mut add = |cell, heat: f64| -> Result<(), Error> {
            checkpoint()?;
            let power = finite(heat / dt)?;
            if heat < 0.0 || (heat != 0.0 && power == 0.0) {
                return Err(invalid(
                    "represented winding heat disappeared in power conversion",
                ));
            }
            let sum = powers.entry(cell).or_insert(0.0);
            *sum = finite(*sum + power)?;
            Ok(())
        };
        for loss in &step.cells {
            let cell = *self
                .windings
                .get(&loss.source.triangle)
                .ok_or_else(|| invalid("native conductor heat has no destination"))?;
            add(cell, loss.heat_j)?;
        }
        for loss in &step.graph_resistors {
            let cell = *self
                .graph
                .get(&loss.branch)
                .ok_or_else(|| invalid("native graph-resistor heat has no destination"))?;
            add(cell, loss.heat_j)?;
        }
        Ok(powers)
    }
}
