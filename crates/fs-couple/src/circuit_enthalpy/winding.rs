//! Spatial winding Joule heat on the original finite enthalpy/Newton owner.
//!
//! The magnetic field/basis and native Kirchhoff history remain with fs-em/fs-phs.
//! This adapter reuses the existing MaterialCell and prepared thermal-transport
//! bindings: caloric storage, conductor resistivity and thermal conductivity
//! must come from the SAME cards and actual endpoint temperatures. Mapped
//! magnetic triangles cover whole fixed thermal volumes; no guessed lumped resistor,
//! fill fraction, copper coefficient or alternate heat integrator is introduced.

/// Circuit-driven nonlinear fields with the same material and heat receivers.
pub mod nonlinear;

use super::magnetic::{
    MagneticRoutes, MagneticThermalNetwork, check, invalid, material, same_material,
};
use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ResistorCell, ThermalDriver, cell_driver_candidate, finite,
    material::{ConductionLink, MaterialCell, MaterialNetwork, MaterialNetworkSample},
};
use fs_em::winding::ohmic::{
    OhmicWindingCandidate, OhmicWindingCircuit, OhmicWindingControls, OhmicWindingError,
    OhmicWindingStep, WindingConductor,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet},
};

pub use super::magnetic::{
    MagneticRegion, MagneticThermalError as WindingThermalError, MagneticThermalSample,
};

/// Original caloric/conductive samples and every magnetic region's source.
#[derive(Clone, Debug, PartialEq)]
pub struct WindingThermalSample {
    pub network: MaterialNetworkSample,
    pub magnetic: Vec<MagneticThermalSample>,
}

/// Data-only material, geometry and complete heat routing for stranded windings.
/// Thermal cells outside the magnetic geometry keep their original latent and
/// boundary behavior. Mapped magnetic cells admit fixed solid geometry/density;
/// thermal expansion, changing magnetic constitutive data and moving conductors
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

type Error = CircuitEnthalpyError<WindingThermalError>;
fn electrical(error: OhmicWindingError) -> Error {
    match error {
        OhmicWindingError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(WindingThermalError::Winding(other)),
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
            numerical_dissipation_j: 0.0,
            stored_change_j: step.stored_after_j - step.stored_before_j,
            supplied_j: step.native.supplied,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        OhmicWindingCandidate::into_parts(self)
    }
}

struct Routes {
    magnetic: MagneticRoutes,
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

    fn magnetic_layout(&self) -> MagneticThermalNetwork<'_> {
        MagneticThermalNetwork {
            cells: self.cells,
            links: self.links,
            magnetic_regions: self.magnetic_regions,
            volume_tolerance_m3: self.volume_tolerance_m3,
        }
    }

    fn admit_routes(
        &self,
        thermal: &DepositionSimulation<'_>,
        circuit: &OhmicWindingCircuit<'_, '_>,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Routes, Error> {
        let triangles: Vec<_> = circuit.conductors().iter().map(|s| s.triangle).collect();
        self.admit_routes_for_mesh(
            thermal,
            circuit.magnetic().basis().problem().mesh(),
            &triangles,
            circuit.magnetic().graph(),
            cancelled,
        )
    }

    // Linear and nonlinear circuit fields share the original volumes and heat routes.
    fn admit_routes_for_mesh(
        &self,
        thermal: &DepositionSimulation<'_>,
        mesh: &fs_rep_mesh::TriComplex2,
        conducting_triangles: &[usize],
        circuit: &CircuitGraph,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Routes, Error> {
        let magnetic = self.magnetic_layout().routes_for_mesh(
            thermal,
            mesh,
            conducting_triangles,
            cancelled,
        )?;
        if self.graph_resistors.len()
            != circuit
                .branches
                .iter()
                .filter(|(_, _, b)| matches!(b, Branch::Resistor { .. }))
                .count()
        {
            return Err(invalid("every graph resistor needs one heat destination"));
        }
        let mut graph = BTreeMap::new();
        for route in self.graph_resistors {
            check(cancelled)?;
            if !thermal.cells().contains_key(&route.cell)
                || graph.insert(route.branch, route.cell).is_some()
                || !matches!(
                    circuit.branches.get(route.branch),
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
            windings: magnetic.conductors.clone(),
            magnetic,
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
        let conductors = circuit
            .conductors()
            .iter()
            .map(|s| (s.triangle, &s.conductor))
            .collect();
        let (sources, samples) = self.magnetic_layout().materials_at(
            thermal,
            circuit.magnetic().basis().problem(),
            &conductors,
            sample,
            &routes.magnetic,
            cancelled,
        )?;
        Ok((
            sources
                .into_iter()
                .map(|(triangle, conductor)| WindingConductor {
                    triangle,
                    conductor,
                })
                .collect(),
            samples,
        ))
    }
}

// Refresh pins the exact selected claims. Its honest usage receipt therefore
// changes policy even at an unchanged physical state; receipt identity is not
// an equality test for accepted temperature. Keep both full receipts intact.
fn same_conductor_state(a: &WindingConductor, b: &WindingConductor) -> bool {
    a.triangle == b.triangle
        && a.conductor.source_basis() == b.conductor.source_basis()
        && a.conductor.resistivity_ohm_m().to_bits() == b.conductor.resistivity_ohm_m().to_bits()
        && same_material(a.conductor.material(), b.conductor.material())
}

impl Routes {
    fn deposit(
        &self,
        candidate: &OhmicWindingCandidate<'_, '_>,
        dt: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Error> {
        let step = candidate.step();
        self.deposit_losses(&step.cells, &[], &step.graph_resistors, dt, checkpoint)
    }

    // Only physical winding/core/graph losses enter the thermal residual.
    fn deposit_losses(
        &self,
        winding: &[fs_em::winding::ohmic::WindingCellHeat],
        core: &[fs_em::planar::eddy::EddyCellHeat],
        graph: &[fs_em::winding::ohmic::WindingGraphResistorHeat],
        dt: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Error> {
        if winding.len().checked_add(core.len()) != Some(self.windings.len())
            || graph.len() != self.graph.len()
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
        let mut triangles = BTreeSet::new();
        for (triangle, heat) in winding
            .iter()
            .map(|loss| (loss.source.triangle, loss.heat_j))
            .chain(core.iter().map(|loss| (loss.source.triangle, loss.heat_j)))
        {
            if !triangles.insert(triangle) {
                return Err(invalid(
                    "winding and core heat must have disjoint original triangles",
                ));
            }
            let cell = *self
                .windings
                .get(&triangle)
                .ok_or_else(|| invalid("native conductor heat has no destination"))?;
            add(cell, heat)?;
        }
        let mut branches = BTreeSet::new();
        for loss in graph {
            if !branches.insert(loss.branch) {
                return Err(invalid(
                    "native graph-resistor heat destination is duplicated",
                ));
            }
            let cell = *self
                .graph
                .get(&loss.branch)
                .ok_or_else(|| invalid("native graph-resistor heat has no destination"))?;
            add(cell, loss.heat_j)?;
        }
        Ok(powers)
    }
}
