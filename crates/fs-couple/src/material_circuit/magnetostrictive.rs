//! Reciprocal magnetostrictive motion inside the original material RLC graph.
//!
//! One inductor is the ideal coil of an existing MagnetostrictiveRod. Its
//! uncoupled lambda^2/(2 L) storage is REPLACED by the rod's original coupled
//! strain/flux and kinetic energy. The existing Kirchhoff descriptor advances
//! that same energy with all other charges, fluxes, sources and resistors.
//! Source-derived stress drives physical axial motion; motion changes current
//! through the SAME energy derivative, without an extra actuator/sensor gain.
//!
//! This entire network uses PERTURBATIONS about the rod's held DC-biased
//! reference. Sources and branch voltages/currents are increments. Resistance
//! dissipates R*(Delta i)^2 in this incremental account. That is NOT generally
//! the total R*(I_bias+Delta i)^2 heat of a biased winding, or its change from
//! the DC heating baseline. Bias preparation, bias supply and their cross terms
//! remain outside this incremental owner. The `thermal` child adds actual
//! total winding heat, terminal work and finite winding enthalpy.
//! Fixed source context, ideal uniform linkage, small axial deformation and the
//! original strain/H envelope apply. Hysteresis, eddy loss, magnetic-path motion
//! and acoustic radiation remain separate physical owners.

pub mod thermal;

use super::{
    MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor, ResistorHeat,
    account_step, bind, extended, poll, positive,
};
use crate::magnetostrictive::{MagnetostrictiveCouplingError, MagnetostrictiveRod};
use fs_phs::circuit::{Branch, CircuitDae, CircuitError, CircuitGraph};
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, StepRecord, Storage, step_descriptor};

/// Original source/rod or native circuit refusal; no accepted history is changed.
#[derive(Debug)]
pub enum MagnetostrictiveCircuitError {
    Circuit(MaterialCircuitError),
    Magnetostrictive(MagnetostrictiveCouplingError),
}
impl From<MaterialCircuitError> for MagnetostrictiveCircuitError {
    fn from(error: MaterialCircuitError) -> Self {
        Self::Circuit(error)
    }
}
impl From<CircuitError> for MagnetostrictiveCircuitError {
    fn from(error: CircuitError) -> Self {
        Self::Circuit(error.into())
    }
}
impl From<MagnetostrictiveCouplingError> for MagnetostrictiveCircuitError {
    fn from(error: MagnetostrictiveCouplingError) -> Self {
        Self::Magnetostrictive(error)
    }
}
impl core::fmt::Display for MagnetostrictiveCircuitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "magnetostrictive material circuit: {self:?}")
    }
}
impl std::error::Error for MagnetostrictiveCircuitError {}
type Result<T> = core::result::Result<T, MagnetostrictiveCircuitError>;
fn invalid(message: &'static str) -> MagnetostrictiveCircuitError {
    MaterialCircuitError::Invalid(message).into()
}
fn native(error: fs_phs::PhsError) -> MagnetostrictiveCircuitError {
    CircuitError::from(error).into()
}

/// Joint accepted incremental electrical and physical axial mechanical history.
#[derive(Clone, Debug, PartialEq)]
pub struct MagnetostrictiveCircuit {
    circuit: MaterialCircuit,
    rod: MagnetostrictiveRod,
    branch: usize,
    flux_slot: usize,
    extension_m: f64,
    momentum_kg_m_s: f64,
}

/// One incremental work account, excluding the fixed bias and its preparation.
pub struct MagnetostrictiveCircuitStep {
    /// Original descriptor result; electrical ports followed by the force port.
    pub native: StepRecord,
    /// Combined rod/mechanical and other RLC perturbation storage [J].
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Original coupled rod energy, including kinetic energy, counted once.
    pub rod_incremental_energy_j: [f64; 2],
    /// Native C/L energy EXCLUDING the selected coil's replaced storage.
    pub other_incremental_electrical_energy_j: [f64; 2],
    pub incremental_electrical_supply_j: f64,
    pub incremental_mechanical_supply_j: f64,
    /// Original source-bearing branch records in the INCREMENTAL convention.
    /// Their voltage/current are increments and their heat_j is R*(Delta i)^2 dt;
    /// never deposit these as total DC-biased winding heat without its own owner.
    pub incremental_resistors: Vec<ResistorHeat>,
    /// Sum of those incremental resistance losses, independently checked.
    pub incremental_resistor_dissipation_j: f64,
    pub incremental_energy_residual_j: f64,
}

/// Inspectable proposal; joint publication consumes the candidate.
pub struct MagnetostrictiveCircuitCandidate {
    model: MagnetostrictiveCircuit,
    step: MagnetostrictiveCircuitStep,
}
impl MagnetostrictiveCircuitCandidate {
    pub const fn model(&self) -> &MagnetostrictiveCircuit {
        &self.model
    }
    pub const fn step(&self) -> &MagnetostrictiveCircuitStep {
        &self.step
    }
    pub fn into_parts(self) -> (MagnetostrictiveCircuit, MagnetostrictiveCircuitStep) {
        (self.model, self.step)
    }
}

impl MagnetostrictiveCircuit {
    /// Bind one original inductor to the existing source-derived coil/core rod.
    ///
    /// Its numeric inductance is a placeholder, replaced by the material's
    /// clamped inductance. Positive first-to-second branch current produces +3 H.
    /// Fluxes and charges retain their original respective branch orders; the
    /// selected flux is the rod's linked perturbation. Source-port gaps retain
    /// their original indices. The true coupled energy admits initial ideal
    /// constraints, never the uncoupled lambda/L approximation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mut graph: CircuitGraph,
        sources: Vec<MaterialResistor>,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        branch: usize,
        rod: MagnetostrictiveRod,
        extension_m: f64,
        velocity_m_s: f64,
        initial_ports: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if graph.node_count > super::MAX_MATERIAL_CIRCUIT_NODES
            || graph.branches.len() > super::MAX_MATERIAL_CIRCUIT_BRANCHES
        {
            return Err(invalid(
                "magnetostrictive graph exceeds bounded native circuit capacity",
            ));
        }
        let entry = graph
            .branches
            .get_mut(branch)
            .ok_or_else(|| invalid("magnetostrictive coil branch is outside the original graph"))?;
        if !matches!(entry.2, Branch::Inductor { .. }) {
            return Err(invalid(
                "magnetostrictive coil must replace an inductor branch",
            ));
        }
        entry.2 = Branch::Inductor {
            henries: rod.clamped_inductance_h(),
        };
        let flux_slot = graph.branches[..branch]
            .iter()
            .filter(|(_, _, element)| matches!(element, Branch::Inductor { .. }))
            .count();
        let (graph, resistors, dae) = bind(&graph, sources, &[])?;
        let momentum_kg_m_s = rod.modal_mass_kg() * velocity_m_s;
        if !velocity_m_s.is_finite()
            || !momentum_kg_m_s.is_finite()
            || velocity_m_s != 0.0 && momentum_kg_m_s == 0.0
        {
            return Err(invalid("magnetostrictive momentum is unrepresentable"));
        }
        let circuit = MaterialCircuit {
            graph,
            resistors,
            contacts: Vec::new(),
            fluxes_wb,
            charges_c,
            storage: Vec::new(),
        };
        let model = Self {
            circuit,
            rod,
            branch,
            flux_slot,
            extension_m,
            momentum_kg_m_s,
        };
        let (coupled, _) = model.assemble(&dae)?;
        let x = model.extended(&dae, &coupled, initial_ports)?;
        let mut ports = initial_ports.to_vec();
        ports.push(0.0);
        model.admit(&coupled, &x, &ports)?;
        poll(&mut cancelled)?;
        Ok(model)
    }

    /// Source/routing view only. Advance this composite, not the isolated view:
    /// the active inductor's current depends on actual rod strain.
    pub const fn circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    pub const fn rod(&self) -> &MagnetostrictiveRod {
        &self.rod
    }
    pub const fn branch(&self) -> usize {
        self.branch
    }
    pub const fn extension_m(&self) -> f64 {
        self.extension_m
    }
    pub const fn momentum_kg_m_s(&self) -> f64 {
        self.momentum_kg_m_s
    }
    pub fn velocity_m_s(&self) -> f64 {
        self.momentum_kg_m_s / self.rod.modal_mass_kg()
    }
    pub fn flux_linkage_increment_wb(&self) -> f64 {
        self.circuit.fluxes_wb[self.flux_slot]
    }
    pub fn current_increment_a(&self) -> Result<f64> {
        Ok(self.rod.length_m()
            * self
                .rod
                .response(self.extension_m, self.flux_linkage_increment_wb())?
                .field_increment_a_m
            / f64::from(self.rod.turns()))
    }
    pub fn rod_incremental_energy_j(&self) -> Result<f64> {
        Ok(self.rod.energy(&[
            self.extension_m,
            self.momentum_kg_m_s,
            self.flux_linkage_increment_wb(),
        ])?)
    }

    /// Advance circuit and physical rod together with the original descriptor.
    ///
    /// Source, geometry and bias of the rod remain fixed. Electrical and force
    /// inputs are held perturbations. Native midpoint efforts determine circuit
    /// resistance loss and reciprocal force/current. Fresh original source-
    /// envelope, algebraic and independent work admission precedes publication.
    /// The midpoint reduction is second order for this fixed linear tangent and
    /// held inputs; cut source discontinuities and refine the time step to resolve
    /// motion. An energy residual alone is not a trajectory-accuracy certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(
        &self,
        sources: Vec<MaterialResistor>,
        electrical_ports: &[f64],
        force_increment_n: f64,
        duration_s: f64,
        controls: MaterialCircuitSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MagnetostrictiveCircuitCandidate> {
        poll(&mut cancelled)?;
        if !force_increment_n.is_finite()
            || !positive(duration_s)
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.power_tolerance_w)
        {
            return Err(invalid(
                "finite force and positive duration/work/power tolerances required",
            ));
        }
        let (graph, resistors, dae) = bind(&self.circuit.graph, sources, &[])?;
        let (coupled, law) = self.assemble(&dae)?;
        let x0 = self.extended(&dae, &coupled, electrical_ports)?;
        let mut ports = electrical_ports.to_vec();
        ports.push(force_increment_n);
        self.admit(&coupled, &x0, &ports)?;
        poll(&mut cancelled)?;
        let record = step_descriptor(&coupled.system, &x0, &ports, duration_s).map_err(native)?;
        let old_d = dae.system.differential_dim();
        let rod_incremental_energy_j = [law.component_energy(&x0), law.component_energy(&record.x)];
        let other_incremental_electrical_energy_j =
            [law.other.hamiltonian(&x0), law.other.hamiltonian(&record.x)];
        let stored_before_j =
            rod_incremental_energy_j[0] + other_incremental_electrical_energy_j[0];
        let stored_after_j = rod_incremental_energy_j[1] + other_incremental_electrical_energy_j[1];
        let (incremental_resistors, contacts, incremental_resistor_dissipation_j, native_residual) =
            account_step(
                &graph,
                &resistors,
                &[],
                &coupled.node_potential_index,
                &x0,
                &record,
                stored_before_j,
                stored_after_j,
                duration_s,
                controls,
                &mut cancelled,
            )?;
        if !contacts.is_empty() {
            return Err(invalid("unexpected magnetostrictive contact loss"));
        }
        let incremental_electrical_supply_j = duration_s
            * electrical_ports
                .iter()
                .zip(&record.y)
                .map(|(input, output)| input * output)
                .sum::<f64>();
        let incremental_mechanical_supply_j =
            force_increment_n * (record.x[old_d] - self.extension_m);
        let native_mechanical = duration_s * force_increment_n * record.y[dae.ports];
        let incremental_energy_residual_j = stored_after_j - stored_before_j
            + incremental_resistor_dissipation_j
            - incremental_electrical_supply_j
            - incremental_mechanical_supply_j;
        if [
            incremental_electrical_supply_j,
            incremental_mechanical_supply_j,
            native_mechanical,
            incremental_energy_residual_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || (incremental_mechanical_supply_j - native_mechanical).abs()
                > controls.energy_tolerance_j
            || (incremental_electrical_supply_j + incremental_mechanical_supply_j - record.supplied)
                .abs()
                > controls.energy_tolerance_j
            || (incremental_energy_residual_j - native_residual).abs() > controls.energy_tolerance_j
            || incremental_energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(MaterialCircuitError::Balance(
                "magnetostrictive motion, circuit work and incremental resistor loss do not close",
            )
            .into());
        }
        self.admit(&coupled, &record.x, &ports)?;
        let mut circuit = self.circuit.clone();
        circuit.graph = graph;
        circuit.resistors = resistors;
        circuit.fluxes_wb = dae
            .flux_index
            .iter()
            .map(|&index| record.x[index])
            .collect();
        circuit.charges_c = dae
            .charge_index
            .iter()
            .map(|&index| record.x[index])
            .collect();
        let model = Self {
            circuit,
            rod: self.rod.clone(),
            branch: self.branch,
            flux_slot: self.flux_slot,
            extension_m: record.x[old_d],
            momentum_kg_m_s: record.x[old_d + 1],
        };
        model.rod_incremental_energy_j()?;
        poll(&mut cancelled)?;
        Ok(MagnetostrictiveCircuitCandidate {
            model,
            step: MagnetostrictiveCircuitStep {
                native: record,
                stored_before_j,
                stored_after_j,
                rod_incremental_energy_j,
                other_incremental_electrical_energy_j,
                incremental_electrical_supply_j,
                incremental_mechanical_supply_j,
                incremental_resistors,
                incremental_resistor_dissipation_j,
                incremental_energy_residual_j,
            },
        })
    }

    fn admit(&self, dae: &CircuitDae, x: &[f64], ports: &[f64]) -> Result<()> {
        let admitted = dae.consistent_initial_state(x, ports)?;
        for index in 0..dae.system.differential_dim() {
            if !admitted[index].is_finite()
                || (admitted[index] - x[index]).abs() > 64.0 * f64::EPSILON * x[index].abs()
            {
                return Err(invalid(
                    "magnetostrictive algebraic admission changed conserved state",
                ));
            }
        }
        if !dae.system.hamiltonian(x).is_finite() {
            return Err(invalid(
                "magnetostrictive state left its original source envelope",
            ));
        }
        Ok(())
    }

    fn extended(
        &self,
        original: &CircuitDae,
        coupled: &CircuitDae,
        ports: &[f64],
    ) -> Result<Vec<f64>> {
        let old = extended(
            original,
            &self.circuit.fluxes_wb,
            &self.circuit.charges_c,
            ports,
        )?;
        let d = original.system.differential_dim();
        let mut x = vec![0.0; coupled.system.state_dim()];
        for (index, value) in old.into_iter().enumerate() {
            x[if index < d { index } else { index + 2 }] = value;
        }
        x[d] = self.extension_m;
        x[d + 1] = self.momentum_kg_m_s;
        Ok(x)
    }

    fn assemble(&self, original: &CircuitDae) -> Result<(CircuitDae, CoupledStorage)> {
        let old_n = original.system.state_dim();
        let old_d = original.system.differential_dim();
        let (n, d, m) = (old_n + 2, old_d + 2, original.ports + 1);
        if n > super::MAX_STATES || m > super::MAX_PORTS {
            return Err(invalid(
                "magnetostrictive circuit exceeds bounded descriptor capacity",
            ));
        }
        let flux = *original
            .flux_index
            .get(self.flux_slot)
            .ok_or_else(|| invalid("selected native coil flux coordinate is absent"))?;
        let original_q = original
            .system
            .storage()
            .exact_quadratic_matrix()
            .ok_or_else(|| invalid("original circuit quadratic storage is required"))?;
        let inverse_l = 1.0 / self.rod.clamped_inductance_h();
        if original_q.len() != old_n * old_n
            || original_q[flux * old_n + flux] != inverse_l
            || (0..old_n).any(|index| {
                index != flux
                    && (original_q[flux * old_n + index] != 0.0
                        || original_q[index * old_n + flux] != 0.0)
            })
        {
            return Err(invalid(
                "selected native inductance does not match the rod's clamped coil",
            ));
        }
        let mut q = vec![0.0; d * d];
        for row in 0..old_d {
            for column in 0..old_d {
                q[row * d + column] = original_q[row * old_n + column];
            }
        }
        q[flux * d + flux] = 0.0;
        let other = QuadraticStorage::new(q.clone(), d).map_err(native)?;
        // Exact Hessian for the native allocation-free quadratic roundoff bound;
        // actual energy/efforts still come from the original material/rod owner.
        q[flux * d + flux] = inverse_l;
        q[old_d * d + old_d] = self.rod.fixed_flux_stiffness_n_m();
        q[(old_d + 1) * d + old_d + 1] = 1.0 / self.rod.modal_mass_kg();
        q[old_d * d + flux] = -self.rod.coupling_wb_m() * inverse_l;
        q[flux * d + old_d] = q[old_d * d + flux];
        let law = CoupledStorage {
            other,
            quadratic: QuadraticStorage::new(q, d).map_err(native)?,
            rod: self.rod.clone(),
            flux,
            extension: old_d,
            momentum: old_d + 1,
        };
        let (old_j, old_r, old_g) = original.system.structure();
        let (mut j, mut r, mut g) = (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * m]);
        let map = |index: usize| if index < old_d { index } else { index + 2 };
        for row in 0..old_n {
            for column in 0..old_n {
                j[map(row) * n + map(column)] = old_j[row * old_n + column];
                r[map(row) * n + map(column)] = old_r[row * old_n + column];
            }
            for port in 0..original.ports {
                g[map(row) * m + port] = old_g[row * original.ports + port];
            }
        }
        j[old_d * n + old_d + 1] = 1.0;
        j[(old_d + 1) * n + old_d] = -1.0;
        g[(old_d + 1) * m + original.ports] = 1.0;
        let system = DescriptorPortHamiltonian::new(n, d, m, j, r, g, Box::new(law.clone()))
            .map_err(native)?;
        Ok((
            CircuitDae {
                system,
                flux_index: original.flux_index.clone(),
                charge_index: original.charge_index.clone(),
                node_potential_index: original
                    .node_potential_index
                    .iter()
                    .copied()
                    .map(map)
                    .collect(),
                ports: m,
            },
            law,
        ))
    }
}

#[derive(Clone)]
struct CoupledStorage {
    other: QuadraticStorage,
    quadratic: QuadraticStorage,
    rod: MagnetostrictiveRod,
    flux: usize,
    extension: usize,
    momentum: usize,
}
impl CoupledStorage {
    fn component_energy(&self, x: &[f64]) -> f64 {
        self.rod
            .energy(&[x[self.extension], x[self.momentum], x[self.flux]])
            .unwrap_or(f64::NAN)
    }
}
impl Storage for CoupledStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.other.hamiltonian(x) + self.component_energy(x)
    }
    fn gradient(&self, x: &[f64], output: &mut [f64]) {
        self.other.gradient(x, output);
        match self.rod.response(x[self.extension], x[self.flux]) {
            Ok(response) => {
                output[self.flux] = self.rod.length_m() * response.field_increment_a_m
                    / f64::from(self.rod.turns());
                output[self.extension] = self.rod.area_m2() * response.stress_increment_pa;
                output[self.momentum] = x[self.momentum] / self.rod.modal_mass_kg();
            }
            Err(_) => output.fill(f64::NAN),
        }
    }
    fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
        self.quadratic.exact_quadratic_matrix()
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        Some(
            if [a, b, mid]
                .iter()
                .all(|state| self.hamiltonian(state).is_finite())
            {
                0.0
            } else {
                f64::NAN
            },
        )
    }
}

#[cfg(test)]
mod tests;
