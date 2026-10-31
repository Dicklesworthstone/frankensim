//! A source-bound piezoelectric rod inside the original Kirchhoff descriptor.
//!
//! The selected capacitor branch is the rod's electrode pair, oriented from
//! its initial (+3) electrode to its final electrode. Its clamped capacitance
//! is lowered from the rod's original material. Its quadratic charge energy is
//! REPLACED by the original charge/strain energy, not added to it. All other
//! circuit stamps, fluxes, charges and source ports retain their native owner.
//! Ordinary stepping holds the rod source context fixed. The explicit thermal
//! entry rebinds that same rod with frozen conserved coordinates and records
//! parameter work for its caloric owner. Dielectric loss, depoling and nonlinear
//! strain are not added.

use super::{MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve,
    MaterialResistor, ResistorHeat, account_step, bind, extended, poll, positive};
use crate::piezoelectric::{PiezoelectricCouplingError, PiezoelectricRod};
use fs_material::piezoelectric::thermal::ThermalLongitudinalPiezoelectric;
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, StepRecord, Storage, step_descriptor};
use fs_phs::circuit::{Branch, CircuitDae, CircuitError, CircuitGraph};

impl From<PiezoelectricCouplingError> for MaterialCircuitError {
    fn from(error: PiezoelectricCouplingError) -> Self { Self::Piezoelectric(error) }
}
fn invalid(message: &'static str) -> MaterialCircuitError { MaterialCircuitError::Invalid(message) }
fn native(error: fs_phs::PhsError) -> MaterialCircuitError { CircuitError::from(error).into() }

/// Accepted circuit storage and source-derived mechanical state.
#[derive(Clone, Debug, PartialEq)]
pub struct PiezoelectricCircuit {
    circuit: MaterialCircuit,
    rod: PiezoelectricRod,
    branch: usize,
    charge_slot: usize,
    extension_m: f64,
    momentum_kg_m_s: f64,
}

/// Complete work account; piezoelectric energy includes both its electrical and
/// mechanical storage, whose reversible exchange is internal to this component.
pub struct PiezoelectricCircuitStep {
    pub native: StepRecord,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub piezoelectric_energy_j: [f64; 2],
    pub other_electrical_energy_j: [f64; 2],
    pub electrical_supply_j: f64,
    pub mechanical_supply_j: f64,
    /// Frozen-coordinate ceramic Helmholtz change; internal reversible transfer
    /// for a thermal owner, not electrical supply or irreversible heat.
    pub parameter_work_j: f64,
    pub resistors: Vec<ResistorHeat>,
    pub heat_j: f64,
    pub energy_residual_j: f64,
}

/// One proposed circuit/rod history; construction never mutates accepted inputs.
pub struct PiezoelectricCircuitCandidate {
    model: PiezoelectricCircuit,
    step: PiezoelectricCircuitStep,
}
impl PiezoelectricCircuitCandidate {
    pub const fn model(&self) -> &PiezoelectricCircuit { &self.model }
    pub const fn step(&self) -> &PiezoelectricCircuitStep { &self.step }
    pub fn into_parts(self) -> (PiezoelectricCircuit, PiezoelectricCircuitStep) {
        (self.model, self.step)
    }
}

impl PiezoelectricCircuit {
    /// Bind the selected native capacitor branch to the original piezo material.
    ///
    /// Numeric capacitance on that branch is a placeholder and is replaced.
    /// Charges/fluxes remain in their original respective branch orders. All
    /// resistor values come from their material bindings. Initial algebraic
    /// admission uses the TRUE coupled energy, never an uncoupled Q/C voltage;
    /// incompatible ideal clamps refuse instead of projecting charge or strain.
    /// The rod's +3 direction and electrode order follow the selected branch.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mut graph: CircuitGraph, sources: Vec<MaterialResistor>,
        fluxes_wb: Vec<f64>, charges_c: Vec<f64>, branch: usize,
        rod: PiezoelectricRod, extension_m: f64, velocity_m_s: f64,
        initial_ports: &[f64], mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialCircuitError> {
        poll(&mut cancelled)?;
        if graph.node_count > super::MAX_MATERIAL_CIRCUIT_NODES
            || graph.branches.len() > super::MAX_MATERIAL_CIRCUIT_BRANCHES {
            return Err(invalid("piezoelectric graph exceeds bounded circuit capacity"));
        }
        let entry = graph.branches.get_mut(branch)
            .ok_or_else(|| invalid("piezoelectric branch is outside the original graph"))?;
        if !matches!(entry.2, Branch::Capacitor { .. }) {
            return Err(invalid("piezoelectric electrode pair must replace a capacitor branch"));
        }
        entry.2 = Branch::Capacitor { farads: rod.axially_clamped_capacitance_f() };
        let charge_slot = graph.branches[..branch].iter()
            .filter(|(_, _, element)| matches!(element, Branch::Capacitor { .. })).count();
        let (graph, resistors, dae) = bind(&graph, sources, &[])?;
        let momentum_kg_m_s = rod.modal_mass_kg() * velocity_m_s;
        if !velocity_m_s.is_finite() || !momentum_kg_m_s.is_finite()
            || velocity_m_s != 0.0 && momentum_kg_m_s == 0.0 {
            return Err(invalid("piezoelectric momentum is not representable"));
        }
        // The native parent owns this exact storage representation. Only the
        // composite descriptor below may admit/advance the active capacitor.
        let circuit = MaterialCircuit {
            graph, mutual: Vec::new(), magnetic_couplings: Vec::new(),
            capacitor_curves: Vec::new().into(),
            resistors, contacts: Vec::new(), fluxes_wb, charges_c, storage: Vec::new(),
        };
        let model = Self { circuit, rod, branch, charge_slot, extension_m, momentum_kg_m_s };
        let (coupled, _) = model.assemble(&dae)?;
        let x = model.extended(&dae, &coupled, initial_ports)?;
        let mut ports = initial_ports.to_vec();
        ports.push(0.0);
        model.admit(&coupled, &x, &ports)?;
        poll(&mut cancelled)?;
        Ok(model)
    }

    /// Routing/source view only. Advance this composite, not the isolated view:
    /// the selected capacitor's voltage also depends on the rod's strain.
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    pub const fn rod(&self) -> &PiezoelectricRod { &self.rod }
    pub const fn branch(&self) -> usize { self.branch }
    pub const fn extension_m(&self) -> f64 { self.extension_m }
    pub const fn momentum_kg_m_s(&self) -> f64 { self.momentum_kg_m_s }
    pub fn velocity_m_s(&self) -> f64 { self.momentum_kg_m_s / self.rod.modal_mass_kg() }
    pub fn charge_c(&self) -> f64 { self.circuit.charges_c[self.charge_slot] }
    pub fn voltage_v(&self) -> Result<f64, MaterialCircuitError> {
        Ok(self.rod.length_m() * self.rod.response(self.extension_m, self.charge_c())?.electric_field_v_m)
    }
    /// Original rod material/kinetic storage, counted exactly once.
    pub fn piezoelectric_energy_j(&self) -> Result<f64, MaterialCircuitError> {
        Ok(self.rod.energy(&[self.extension_m, self.momentum_kg_m_s, self.charge_c()])?)
    }

    /// Refresh the SAME ceramic from its pinned affine thermal law while
    /// retaining original extension, momentum, charge, flux and circuit topology.
    /// Then run the original descriptor interval at the proposed endpoint T.
    /// No switching impulse is inferred: incompatible ideal-source constraints
    /// at the frozen-coordinate rebind still refuse through native admission.
    /// The caloric caller must pay -P-delta(TS) and close physical internal energy.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_thermal_material_candidate(
        &self, law: &ThermalLongitudinalPiezoelectric, temperature_k: f64,
        sources: Vec<MaterialResistor>, electrical_ports: &[f64], force_n: f64,
        duration_s: f64, controls: MaterialCircuitSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PiezoelectricCircuitCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        law.at_material(self.rod.material()).map_err(PiezoelectricCouplingError::from)?;
        let at = law.at_temperature(temperature_k).map_err(PiezoelectricCouplingError::from)?;
        let rod = PiezoelectricRod::new(at.material().clone(), self.rod.length_m(),
            self.rod.area_m2(), self.rod.limits())?;
        if rod.modal_mass_kg().to_bits() != self.rod.modal_mass_kg().to_bits() {
            return Err(invalid("thermal piezoelectric update cannot change reference inertia"));
        }
        let volume = self.rod.length_m()*self.rod.area_m2();
        let before = volume*self.rod.response(self.extension_m,self.charge_c())?.stored_energy_j_m3;
        let frozen = volume*rod.response(self.extension_m,self.charge_c())?.stored_energy_j_m3;
        let parameter_work_j = frozen-before;
        if !before.is_finite() || !frozen.is_finite() || !parameter_work_j.is_finite() {
            return Err(invalid("piezoelectric thermal parameter work is unrepresentable"));
        }
        let mut rebound = self.clone();
        rebound.circuit.graph.branches[self.branch].2 = Branch::Capacitor {
            farads: rod.axially_clamped_capacitance_f(),
        };
        rebound.rod = rod;
        let mut candidate = rebound.step_candidate(sources,electrical_ports,force_n,
            duration_s,controls,&mut cancelled)?;
        candidate.step.piezoelectric_energy_j[0] = self.piezoelectric_energy_j()?;
        candidate.step.stored_before_j = candidate.step.piezoelectric_energy_j[0]
            + candidate.step.other_electrical_energy_j[0];
        candidate.step.parameter_work_j = parameter_work_j;
        candidate.step.energy_residual_j = candidate.step.stored_after_j-candidate.step.stored_before_j
            + candidate.step.heat_j-candidate.step.electrical_supply_j-candidate.step.mechanical_supply_j
            - parameter_work_j;
        if !candidate.step.energy_residual_j.is_finite()
            || candidate.step.energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(MaterialCircuitError::Balance("piezoelectric original storage and thermal parameter work do not close"));
        }
        poll(&mut cancelled)?;
        Ok(candidate)
    }

    /// Solve the material circuit and rod together for one held-load interval.
    ///
    /// Resistors are refreshed from their sources for this trial. Every trial
    /// starts from the same accepted charge, flux, extension and momentum.
    /// Actual native algebraic midpoint voltages determine Joule heat. The
    /// signed force-displacement work plus native electrical source work close
    /// against combined storage and that heat before any endpoint is returned.
    pub fn step_candidate(
        &self, sources: Vec<MaterialResistor>, electrical_ports: &[f64],
        force_n: f64, duration_s: f64, controls: MaterialCircuitSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PiezoelectricCircuitCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        if !force_n.is_finite() || !positive(duration_s)
            || !positive(controls.energy_tolerance_j) || !positive(controls.power_tolerance_w) {
            return Err(invalid("finite force and positive duration/work/power tolerances required"));
        }
        let (graph, resistors, dae) = bind(&self.circuit.graph, sources, &[])?;
        let (coupled, law) = self.assemble(&dae)?;
        let x0 = self.extended(&dae, &coupled, electrical_ports)?;
        let mut ports = electrical_ports.to_vec();
        ports.push(force_n);
        self.admit(&coupled, &x0, &ports)?;
        poll(&mut cancelled)?;
        let record = step_descriptor(&coupled.system, &x0, &ports, duration_s).map_err(native)?;
        let old_d = dae.system.differential_dim();
        let piezoelectric_energy_j = [law.component_energy(&x0), law.component_energy(&record.x)];
        let other_electrical_energy_j = [law.other.hamiltonian(&x0), law.other.hamiltonian(&record.x)];
        let stored_before_j = piezoelectric_energy_j[0] + other_electrical_energy_j[0];
        let stored_after_j = piezoelectric_energy_j[1] + other_electrical_energy_j[1];
        let (heat, contacts, heat_j, residual) = account_step(
            &graph, &resistors, &[], &coupled.node_potential_index, &x0, &record,
            stored_before_j, stored_after_j, duration_s, controls, &mut cancelled,
        )?;
        if !contacts.is_empty() { return Err(invalid("unexpected piezoelectric contact heat")); }
        let electrical_supply_j = duration_s * electrical_ports.iter().zip(&record.y)
            .map(|(input, output)| input * output).sum::<f64>();
        let mechanical_supply_j = force_n * (record.x[old_d] - self.extension_m);
        let native_mechanical = duration_s * force_n * record.y[dae.ports];
        let energy_residual_j =
            stored_after_j - stored_before_j + heat_j - electrical_supply_j - mechanical_supply_j;
        if [electrical_supply_j, mechanical_supply_j, native_mechanical, energy_residual_j]
            .iter().any(|value| !value.is_finite())
            || (mechanical_supply_j-native_mechanical).abs() > controls.energy_tolerance_j
            || (electrical_supply_j+mechanical_supply_j-record.supplied).abs() > controls.energy_tolerance_j
            || (energy_residual_j-residual).abs() > controls.energy_tolerance_j
            || energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(MaterialCircuitError::Balance("piezoelectric load work, circuit supply and heat disagree"));
        }
        // Reject incompatible endpoint constraints without changing conserved
        // coordinates. Native rank-aware algebraic admission remains the owner.
        self.admit(&coupled, &record.x, &ports)?;
        let mut circuit = self.circuit.clone();
        circuit.graph = graph;
        circuit.resistors = resistors;
        circuit.fluxes_wb = dae.flux_index.iter().map(|&index| record.x[index]).collect();
        circuit.charges_c = dae.charge_index.iter().map(|&index| record.x[index]).collect();
        let model = Self {
            circuit, rod: self.rod.clone(), branch: self.branch, charge_slot: self.charge_slot,
            extension_m: record.x[old_d], momentum_kg_m_s: record.x[old_d+1],
        };
        model.piezoelectric_energy_j()?;
        poll(&mut cancelled)?;
        Ok(PiezoelectricCircuitCandidate { model, step: PiezoelectricCircuitStep {
            native: record, stored_before_j, stored_after_j, piezoelectric_energy_j,
            other_electrical_energy_j, electrical_supply_j, mechanical_supply_j, parameter_work_j: 0.0,
            resistors: heat, heat_j, energy_residual_j,
        } })
    }

    fn admit(&self, dae: &CircuitDae, x: &[f64], ports: &[f64])
        -> Result<(), MaterialCircuitError> {
        let admitted = dae.consistent_initial_state(x, ports)?;
        for index in 0..dae.system.differential_dim() {
            if !admitted[index].is_finite()
                || (admitted[index]-x[index]).abs() > 64.0*f64::EPSILON*x[index].abs() {
                return Err(invalid("piezoelectric algebraic admission changed conserved state"));
            }
        }
        if !dae.system.hamiltonian(x).is_finite() {
            return Err(invalid("piezoelectric state left its original source envelope"));
        }
        Ok(())
    }

    fn extended(&self, original: &CircuitDae, coupled: &CircuitDae, ports: &[f64])
        -> Result<Vec<f64>, MaterialCircuitError> {
        let old = extended(original, &self.circuit.fluxes_wb, &self.circuit.charges_c, ports)?;
        let d = original.system.differential_dim();
        let mut state = vec![0.0; coupled.system.state_dim()];
        for (index, value) in old.into_iter().enumerate() {
            state[if index < d { index } else { index+2 }] = value;
        }
        state[d] = self.extension_m;
        state[d+1] = self.momentum_kg_m_s;
        Ok(state)
    }

    fn assemble(&self, original: &CircuitDae)
        -> Result<(CircuitDae, PiezoelectricCircuitStorage), MaterialCircuitError> {
        let old_n = original.system.state_dim();
        let old_d = original.system.differential_dim();
        let n = old_n+2;
        let d = old_d+2;
        let m = original.ports+1;
        if n > super::MAX_STATES || m > super::MAX_PORTS {
            return Err(invalid("piezoelectric circuit exceeds bounded native descriptor capacity"));
        }
        let charge = *original.charge_index.get(self.charge_slot)
            .ok_or_else(|| invalid("piezoelectric charge coordinate is absent"))?;
        let original_q = original.system.storage().exact_quadratic_matrix()
            .ok_or_else(|| invalid("original circuit quadratic storage is required"))?;
        let inverse_c = 1.0/self.rod.axially_clamped_capacitance_f();
        if original_q.len() != old_n*old_n || original_q[charge*old_n+charge] != inverse_c
            || (0..old_n).any(|index| index != charge
                && (original_q[charge*old_n+index] != 0.0 || original_q[index*old_n+charge] != 0.0)) {
            return Err(invalid("piezoelectric capacitance does not match the selected native stamp"));
        }
        let mut q = vec![0.0; d*d];
        for row in 0..old_d {
            for column in 0..old_d { q[row*d+column] = original_q[row*old_n+column]; }
        }
        q[charge*d+charge] = 0.0;
        let law = PiezoelectricCircuitStorage {
            other: QuadraticStorage::new(q, d).map_err(native)?,
            rod: self.rod.clone(), charge, extension: old_d, momentum: old_d+1,
        };
        let (old_j, old_r, old_g) = original.system.structure();
        let (mut j, mut r, mut g) = (vec![0.0; n*n], vec![0.0; n*n], vec![0.0; n*m]);
        let map = |index: usize| if index < old_d { index } else { index+2 };
        for row in 0..old_n {
            for column in 0..old_n {
                j[map(row)*n+map(column)] = old_j[row*old_n+column];
                r[map(row)*n+map(column)] = old_r[row*old_n+column];
            }
            for port in 0..original.ports { g[map(row)*m+port] = old_g[row*original.ports+port]; }
        }
        j[old_d*n+old_d+1] = 1.0;
        j[(old_d+1)*n+old_d] = -1.0;
        g[(old_d+1)*m+original.ports] = 1.0;
        let system = DescriptorPortHamiltonian::new(n, d, m, j, r, g, Box::new(law.clone())).map_err(native)?;
        Ok((CircuitDae {
            system, flux_index: original.flux_index.clone(), charge_index: original.charge_index.clone(),
            node_potential_index: original.node_potential_index.iter().copied().map(map).collect(),
            ports: m,
        }, law))
    }
}

#[derive(Clone)]
struct PiezoelectricCircuitStorage {
    other: QuadraticStorage,
    rod: PiezoelectricRod,
    charge: usize,
    extension: usize,
    momentum: usize,
}
impl PiezoelectricCircuitStorage {
    fn component_energy(&self, x: &[f64]) -> f64 {
        self.rod.energy(&[x[self.extension], x[self.momentum], x[self.charge]]).unwrap_or(f64::NAN)
    }
}
impl Storage for PiezoelectricCircuitStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.other.hamiltonian(x) + self.component_energy(x)
    }
    fn gradient(&self, x: &[f64], output: &mut [f64]) {
        self.other.gradient(x, output);
        match self.rod.response(x[self.extension], x[self.charge]) {
            Ok(response) => {
                output[self.charge] = self.rod.length_m() * response.electric_field_v_m;
                output[self.extension] = self.rod.area_m2() * response.stress_pa;
                output[self.momentum] = x[self.momentum]/self.rod.modal_mass_kg();
            }
            Err(_) => output.fill(f64::NAN),
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], midpoint: &[f64]) -> Option<f64> {
        Some(if [a,b,midpoint].iter().all(|state| self.hamiltonian(state).is_finite()) {
            0.0
        } else { f64::NAN })
    }
}
