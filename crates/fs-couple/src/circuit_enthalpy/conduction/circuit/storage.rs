//! Material field storage on the original spatial conductor/circuit/heat path.
//!
//! Dielectric permittivity and core permeability use the SAME source reductions,
//! graph stamping and fixed-coordinate work account as the lumped circuit.
//! Their material states and geometry are held during electrothermal steps;
//! the spatial conductor and other material resistors still evolve with heat.
//! A storage-only reconfiguration records required external parameter work.
//! It does not invent the thermal entropy or mechanical transfer of that change.

use super::*;
use crate::material_circuit::storage::{StorageParameterWork, lower, parameter_work};

/// Proposed source/geometry change with unchanged electrical history and time.
/// Publish it together with the external owner that pays or receives its work.
pub struct ConductionStorageCandidate<C = ResolvedConductorState> {
    circuit: ConductionCircuit<C>,
    work: Vec<StorageParameterWork>,
    parameter_work_j: f64,
}

impl<C: PlanarConductorMaterial> ConductionStorageCandidate<C> {
    pub const fn circuit(&self) -> &ConductionCircuit<C> {
        &self.circuit
    }

    /// Original branch order, using the actual source-derived energy on each side.
    pub fn work(&self) -> &[StorageParameterWork] {
        &self.work
    }

    /// Positive work enters field storage; negative work is returned externally.
    pub const fn parameter_work_j(&self) -> f64 {
        self.parameter_work_j
    }

    pub fn into_parts(self) -> (ConductionCircuit<C>, Vec<StorageParameterWork>, f64) {
        (self.circuit, self.work, self.parameter_work_j)
    }
}

impl<C: PlanarConductorMaterial> ConductionCircuit<C> {
    /// Construct a spatial conductor circuit with source-bound R, C and L.
    ///
    /// The spatial branch is stamped from its actual field extraction; every
    /// other resistor comes from `material_resistors`. Every capacitor and
    /// inductor uses the existing [`MaterialStorage`] geometry/source law.
    /// Numeric graph coefficients are placeholders, never fallbacks. Complete
    /// source records survive all ordinary coupled steps and their continuation.
    ///
    /// This is fixed-state, linear, nondispersive dielectric/core storage. Its
    /// original material law admits the actual context and source applicability;
    /// this path does not freeze unsupported field/frequency-dependent data.
    /// Dielectric/core thermal state is not bound to the conductor's heat cells.
    /// Source-backed C/L stay held throughout each electrothermal interval;
    /// evolving their temperature requires reversible field entropy accounting,
    /// not a change silently charged to the conductor's thermal store.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_material_storage(
        problem: PlanarConduction<C>,
        graph: CircuitGraph,
        terminal: ConductionCircuitTerminal,
        material_resistors: Vec<MaterialResistor>,
        mut storage: Vec<MaterialStorage>,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        field_controls: ConductionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        // Lower storage first so authored placeholder coefficients need not be
        // valid physical values. The original DAE still admits all topology and
        // initial source constraints after every material has stamped its slot.
        let graph = lower(&graph, &mut storage)?;
        let mut circuit = Self::new_with_material_resistors(
            problem,
            graph,
            terminal,
            material_resistors,
            time_s,
            fluxes_wb,
            charges_c,
            ports,
            field_controls,
            &mut cancelled,
        )?;
        circuit.storage = storage;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(circuit)
    }

    /// Full original dielectric/core snapshots and geometry in branch order.
    /// Empty for explicitly authored ideal storage on the original constructors.
    pub fn storage(&self) -> &[MaterialStorage] {
        &self.storage
    }

    /// Propose a source/geometry change at unchanged charge and flux linkage.
    ///
    /// Uses the original circuit's field-energy evaluator and constraint solve.
    /// An ideal-voltage-clamped capacitor cannot generally keep both its charge
    /// and voltage when C changes; that case refuses without an invented impulse.
    /// The spatial mesh, complete conductivity operator, source records, time and
    /// material resistor states remain unchanged, without another field solve.
    ///
    /// The result reports W = delta(q²/(2C) + lambda²/(2L)) as required external
    /// parameter work. It is not Joule heat and is not already deposited into any
    /// thermal cell. This instantaneous control does not implement a mechanical
    /// actuator or temperature-dependent field thermodynamics; those owners must
    /// account for the actual work/entropy when assembling such a coupled model.
    /// A refusal or cancellation leaves all accepted inputs unchanged.
    pub fn rebind_storage_candidate(
        &self,
        mut storage: Vec<MaterialStorage>,
        ports: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionStorageCandidate<C>, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        let count = self
            .graph
            .branches
            .iter()
            .filter(|(_, _, branch)| {
                matches!(branch, Branch::Capacitor { .. } | Branch::Inductor { .. })
            })
            .count();
        if self.storage.len() != count {
            return Err(invalid(
                "live storage changes require source-bound original capacitors and inductors",
            ));
        }
        let graph = lower(&self.graph, &mut storage)?;
        let dae = assemble_circuit(&graph).map_err(electrical)?;
        consistent(&dae, &self.fluxes_wb, &self.charges_c, ports)?;
        let (work, parameter_work_j) = match parameter_work(
            &self.graph,
            &self.storage,
            &storage,
            &self.fluxes_wb,
            &self.charges_c,
            &mut cancelled,
        ) {
            Err(MaterialCircuitError::Cancelled) => return Err(CircuitEnthalpyError::Cancelled),
            other => other?,
        };
        let circuit = Self {
            problem: self.problem.clone(),
            graph,
            terminal: self.terminal.clone(),
            fluxes_wb: self.fluxes_wb.clone(),
            charges_c: self.charges_c.clone(),
            time_s: self.time_s,
            material_resistors: self.material_resistors.clone(),
            storage,
        };
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(ConductionStorageCandidate {
            circuit,
            work,
            parameter_work_j,
        })
    }
}
