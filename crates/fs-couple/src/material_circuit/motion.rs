//! Reciprocal electric-field and rigid-body motion on the existing circuit DAE.
//!
//! A uniform, constant-area dielectric gap stores `q^2 d / (2 epsilon A)`.
//! Voltage and attractive force are derivatives of that ONE energy with respect
//! to charge and gap. The same implicit descriptor step therefore evolves the
//! circuit and the moving electrode; motion is not an externally imposed C update.
//! Body mass comes from the existing regional material/inertia compiler.

use fs_matdb::EvaluationDecision;
use fs_material::field_storage::{ELECTRIC_PERMITTIVITY_PROPERTY, ResolvedCapacitor};
use fs_mbd_material::ResolvedBody;
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, StepRecord, Storage, step_descriptor};
use fs_phs::circuit::{Branch, CircuitDae, CircuitError};

use super::storage::MaterialStorageElement;
use super::contacts::{self, ContactHeat, SeriesElectricalContact};
use super::{MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
    ResistorHeat, account_step, bind, extended, poll, positive};

/// A source-bound capacitor and ONE rigid, translationally guided moving body.
/// Positive gap velocity separates the electrodes; the other electrode is fixed.
/// Regional density/geometry, dielectric source, other circuit storage and graph
/// topology remain attached to this accepted state.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialCapacitorMotion {
    circuit: MaterialCircuit,
    body: ResolvedBody,
    storage_slot: usize,
    charge_slot: usize,
    gap_limits_m: [f64; 2],
    momentum_kg_m_s: f64,
}

/// Held external inputs over one common electrical/mechanical interval.
#[derive(Clone, Copy, Debug)]
pub struct CapacitorMotionInput<'a> {
    /// Original circuit port order, unchanged by the added mechanical coordinate.
    pub electrical_ports: &'a [f64],
    /// Force applied by an external owner, positive toward increasing gap [N].
    pub external_force_n: f64,
    /// Common physical duration [s].
    pub duration_s: f64,
}

/// Joint field, mechanical, source and heat accounting from one descriptor step.
pub struct CapacitorMotionStep {
    /// Native state order: original fluxes/charges, gap, momentum, then original
    /// algebraic coordinates. Outputs preserve circuit ports and append the
    /// discrete-gradient velocity, which is conjugate to external force.
    pub native: StepRecord,
    /// Total capacitor/inductor energy before and after the coupled step [J].
    pub electrical_energy_j: [f64; 2],
    /// Actual moving body's kinetic energy before and after the step [J].
    pub kinetic_energy_j: [f64; 2],
    /// Work supplied through original electrical ports [J].
    pub electrical_supply_j: f64,
    /// External mechanical work, force times actual gap change [J].
    pub mechanical_supply_j: f64,
    /// Work transferred from electrical storage to the body [J]. This is internal
    /// exchange, already present in both storage changes; it is not extra heat.
    pub field_work_on_body_j: f64,
    /// Existing material conductor's branch-local heat at native terminal efforts.
    pub resistors: Vec<ResistorHeat>,
    /// Original series-interface losses, separate from bulk Joule heat.
    pub contacts: Vec<ContactHeat>,
    /// Sum of bulk and contact heat [J], checked against native dissipation.
    pub heat_j: f64,
    /// Combined storage change plus heat minus both external work supplies [J].
    pub energy_residual_j: f64,
}

/// A proposed circuit AND moving body. A failure never changes accepted history.
pub struct CapacitorMotionCandidate {
    model: MaterialCapacitorMotion,
    step: CapacitorMotionStep,
}
impl CapacitorMotionCandidate {
    /// Joint endpoint; electrical charge cannot be accepted without the same gap.
    pub const fn model(&self) -> &MaterialCapacitorMotion { &self.model }
    /// Original solver result and independently partitioned physical work.
    pub const fn step(&self) -> &CapacitorMotionStep { &self.step }
    /// Publish the joint model and its work/heat record together.
    pub fn into_parts(self) -> (MaterialCapacitorMotion, CapacitorMotionStep) {
        (self.model, self.step)
    }
}

fn invalid(what: &'static str) -> MaterialCircuitError { MaterialCircuitError::Invalid(what) }
fn phs(error: fs_phs::PhsError) -> MaterialCircuitError { CircuitError::from(error).into() }
fn within(gap: f64, limits: [f64; 2]) -> bool {
    gap.is_finite() && gap > limits[0] && gap < limits[1]
}
fn epsilon_area(capacitor: &ResolvedCapacitor) -> Result<f64, MaterialCircuitError> {
    let property = capacitor.material().property(ELECTRIC_PERMITTIVITY_PROPERTY)
        .ok_or_else(|| invalid("moving capacitor lost its dielectric source"))?;
    // A temperature-dependent free-energy coefficient does not by itself give
    // internal energy: the missing field entropy would need a thermal owner.
    // Do not label that missing reversible heat as conserved mechanical energy.
    if property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity {
        return Err(invalid("conservative gap motion requires a validity-wide constant dielectric coefficient"));
    }
    let value = property.value_si() * capacitor.area_m2();
    if !positive(value) || !positive(1.0 / value) {
        return Err(invalid("dielectric area product is unrepresentable"));
    }
    Ok(value)
}

impl MaterialCapacitorMotion {
    /// Bind an existing material RLC circuit to a source-resolved rigid body.
    /// `branch` is the ORIGINAL graph index of its moving capacitor. Its existing
    /// gap is the initial position; mass is ONLY `body.inertia().mass()`.
    ///
    /// The caller supplies two positive, finite, ordered OPEN gap limits from
    /// the admitted geometry model and an initial velocity [m/s]. Area stays
    /// fixed. Body regions translate together; rotation/deformation are guided
    /// out. Source and geometry swaps before construction change the actual
    /// dynamics without changing this code. Live body/material replacement needs
    /// its own mass, momentum and thermodynamic transaction.
    ///
    /// This reduction assumes ideal conducting electrodes, a uniform linear
    /// nondispersive dielectric, negligible fringing, and a fixed opposing plate.
    /// It supplies no electrode contact, pull-in stop, breakdown, leakage, electrostriction
    /// or electrocaloric law. Dielectric temperature curves are refused because
    /// their reversible heat requires field entropy data. Plate rigidity and
    /// agreement of its geometry with the exposed capacitor area remain explicit
    /// assembly assumptions, as for the existing material piston.
    pub fn new(circuit: MaterialCircuit, branch: usize, body: ResolvedBody,
        velocity_m_s: f64, gap_limits_m: [f64; 2]) -> Result<Self, MaterialCircuitError>
    {
        let slot = circuit.storage.iter().position(|source| source.branch == branch)
            .ok_or_else(|| invalid("moving capacitor requires an original material storage binding"))?;
        let MaterialStorageElement::Capacitor(capacitor) = &circuit.storage[slot].element else {
            return Err(invalid("moving branch must be a capacitor"));
        };
        if gap_limits_m.iter().any(|&v| !positive(v))
            || !within(capacitor.gap_m(), gap_limits_m) || !velocity_m_s.is_finite() {
            return Err(invalid("finite velocity and positive open gap limits containing the capacitor required"));
        }
        let eps_area = epsilon_area(capacitor)?;
        if gap_limits_m.iter().any(|gap| !positive(eps_area / gap) || !positive(gap / eps_area)) {
            return Err(invalid("capacitance is unrepresentable within the declared gap domain"));
        }
        let mass = body.inertia().mass();
        let momentum = mass * velocity_m_s;
        if !positive(mass) || !positive(1.0 / mass) || !momentum.is_finite()
            || (velocity_m_s != 0.0 && momentum == 0.0) {
            return Err(invalid("source-derived plate mass or momentum is unrepresentable"));
        }
        let charge_slot = circuit.graph.branches[..branch].iter()
            .filter(|(_, _, b)| matches!(b, Branch::Capacitor { .. })).count();
        let result = Self { circuit, body, storage_slot: slot, charge_slot,
            gap_limits_m, momentum_kg_m_s: momentum };
        result.capacitor().stored_energy_j(result.charge_c())?;
        result.kinetic_energy()?;
        result.electrical_force_n()?;
        Ok(result)
    }

    /// Existing circuit at the SAME accepted plate gap and electrical coordinates.
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    /// Canonical regional geometry, material receipts and spatial inertia.
    pub const fn body(&self) -> &ResolvedBody { &self.body }
    /// Original dielectric snapshot with its CURRENT physical gap.
    pub fn capacitor(&self) -> &ResolvedCapacitor {
        match &self.circuit.storage[self.storage_slot].element {
            MaterialStorageElement::Capacitor(capacitor) => capacitor,
            MaterialStorageElement::Inductor(_) => unreachable!("private capacitor binding invariant"),
        }
    }
    /// Positive electrode separation [m]. There is no competing position cache.
    pub fn gap_m(&self) -> f64 { self.capacitor().gap_m() }
    /// Canonical translational momentum [kg m/s].
    pub const fn momentum_kg_m_s(&self) -> f64 { self.momentum_kg_m_s }
    /// Physical endpoint velocity [m/s], positive toward increasing separation.
    pub fn velocity_m_s(&self) -> f64 { self.momentum_kg_m_s / self.body.inertia().mass() }
    /// Instantaneous fixed-charge force `-dH/dd = -q^2/(2 epsilon A)` [N].
    /// Its discrete-gradient counterpart, not this endpoint sample, drives a step.
    pub fn electrical_force_n(&self) -> Result<f64, MaterialCircuitError> {
        let q = self.charge_c();
        let force = -0.5 * q * (q / epsilon_area(self.capacitor())?);
        if !force.is_finite() || (q != 0.0 && force == 0.0) {
            return Err(invalid("electric plate force is unrepresentable"));
        }
        Ok(force)
    }
    fn charge_c(&self) -> f64 { self.circuit.charges_c[self.charge_slot] }
    fn kinetic_energy(&self) -> Result<f64, MaterialCircuitError> {
        let energy = 0.5 * self.momentum_kg_m_s * self.velocity_m_s();
        if !energy.is_finite() || (self.momentum_kg_m_s != 0.0 && energy == 0.0) {
            return Err(invalid("plate kinetic energy is unrepresentable"));
        }
        Ok(energy)
    }

    /// Solve electrical storage, plate gap/momentum and KCL/KVL together using
    /// the original descriptor integrator. Resistor material bindings are held
    /// for this trial, as on `MaterialCircuit::step_candidate`; their heat uses
    /// the same native terminal efforts and conductor law.
    ///
    /// The moving capacitor's old quadratic term is REPLACED by `q^2 d/(2 eps A)`.
    /// The original circuit's J/R/G stamps, transformer constraints and remaining
    /// storage are retained. No force gain, separate circuit integrator, spring
    /// constant or manually overridden inertial mass is introduced.
    ///
    /// The endpoint must also admit its unchanged proposed charge/flux under the
    /// original algebraic constraints at the NEW geometry. Incompatible ideal
    /// clamps refuse; no projection, switching impulse or lost charge is hidden.
    /// Out-of-domain Newton trials, gap limits, imbalance and cancellation leave
    /// the whole model unchanged. Work tolerances do not certify time accuracy.
    pub fn step_candidate(&self, sources: Vec<MaterialResistor>, input: CapacitorMotionInput<'_>,
        controls: MaterialCircuitSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<CapacitorMotionCandidate, MaterialCircuitError>
    {
        poll(&mut cancelled)?;
        if !self.circuit.contacts.is_empty() {
            return Err(invalid("contact-bound motion requires explicit refreshed contact sources"));
        }
        self.step_with_contacts_candidate(sources, Vec::new(), input, controls, cancelled)
    }

    /// Advance the same reciprocal field/body transaction with explicit refreshed
    /// series contacts. Each original contact ID and branch assignment is retained.
    /// Bulk and interface heat remain separate, and both enter total dissipation.
    /// Contact source/area changes follow the existing zero-storage interface
    /// model; electrode gap contact, arc energy and interface memory are not added.
    pub fn step_with_contacts_candidate(&self, sources: Vec<MaterialResistor>,
        mut contacts: Vec<SeriesElectricalContact>, input: CapacitorMotionInput<'_>,
        controls: MaterialCircuitSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<CapacitorMotionCandidate, MaterialCircuitError>
    {
        poll(&mut cancelled)?;
        contacts::admit(&mut contacts)?;
        if contacts.len() != self.circuit.contacts.len()
            || contacts.iter().zip(&self.circuit.contacts)
                .any(|(a,b)| (a.id(),a.branch()) != (b.id(),b.branch())) {
            return Err(invalid("electrical contact identities changed inside a physical motion step"));
        }
        let dt = input.duration_s;
        if !positive(dt) || !input.external_force_n.is_finite()
            || !positive(controls.energy_tolerance_j) || !positive(controls.power_tolerance_w) {
            return Err(invalid("finite held force and positive duration/work/power tolerances required"));
        }
        let (graph, resistors, dae) = bind(&self.circuit.graph, sources, &contacts)?;
        let electrical_x0 = extended(&dae, &self.circuit.fluxes_wb, &self.circuit.charges_c,
            input.electrical_ports)?;
        let (system, law) = self.assemble(&dae)?;
        let nd = dae.system.differential_dim();
        let map = |i: usize| if i < nd { i } else { i + 2 };
        let mut x0 = vec![0.0; system.state_dim()];
        for (i, &v) in electrical_x0.iter().enumerate() { x0[map(i)] = v; }
        x0[nd] = self.gap_m();
        x0[nd+1] = self.momentum_kg_m_s;
        let mut ports = input.electrical_ports.to_vec();
        ports.push(input.external_force_n);
        poll(&mut cancelled)?;
        let record = step_descriptor(&system, &x0, &ports, dt).map_err(phs)?;
        let gap = record.x[nd];
        if !within(gap, self.gap_limits_m) {
            return Err(invalid("coupled endpoint left the admitted open gap domain"));
        }
        let electrical = [law.electrical_energy(&x0), law.electrical_energy(&record.x)];
        let kinetic = [law.kinetic_energy(&x0), law.kinetic_energy(&record.x)];
        let nodes: Vec<_> = dae.node_potential_index.iter().map(|&i| map(i)).collect();
        let (heat, contact_heat, heat_j, residual) = account_step(&graph, &resistors, &contacts, &nodes, &x0, &record,
            electrical[0]+kinetic[0], electrical[1]+kinetic[1], dt, controls, &mut cancelled)?;
        let electrical_supply = dt * input.electrical_ports.iter().zip(&record.y)
            .map(|(&u, &y)| u*y).sum::<f64>();
        let mechanical_supply = input.external_force_n * (gap-self.gap_m());
        let native_mechanical_supply = dt * input.external_force_n * record.y[dae.ports];
        let field_work = kinetic[1]-kinetic[0]-mechanical_supply;
        let electrical_work = electrical_supply-heat_j-(electrical[1]-electrical[0]);
        let total_residual = (electrical[1]-electrical[0])+(kinetic[1]-kinetic[0])
            +heat_j-electrical_supply-mechanical_supply;
        if [electrical_supply, mechanical_supply, native_mechanical_supply, field_work,
            electrical_work, total_residual].iter().any(|v| !v.is_finite())
            || (mechanical_supply-native_mechanical_supply).abs() > controls.energy_tolerance_j
            || (field_work-electrical_work).abs() > controls.energy_tolerance_j
            || (total_residual-residual).abs() > controls.energy_tolerance_j
            || total_residual.abs() > controls.energy_tolerance_j {
            return Err(MaterialCircuitError::Balance("reciprocal field work, displacement work and branch heat disagree"));
        }
        let mut storage = self.circuit.storage.clone();
        storage[self.storage_slot].element = MaterialStorageElement::Capacitor(
            ResolvedCapacitor::try_from_material(self.capacitor().material().clone(),
                self.capacitor().area_m2(), gap)?);
        poll(&mut cancelled)?;
        let graph = super::storage::lower(&graph, &mut storage)?;
        let mut circuit = MaterialCircuit::new_with_contacts(graph, resistors, contacts,
            dae.flux_index.iter().map(|&i| record.x[i]).collect(),
            dae.charge_index.iter().map(|&i| record.x[i]).collect(), input.electrical_ports)?;
        circuit.storage = storage;
        let model = Self { circuit, body: self.body.clone(), storage_slot: self.storage_slot,
            charge_slot: self.charge_slot, gap_limits_m: self.gap_limits_m,
            momentum_kg_m_s: record.x[nd+1] };
        model.capacitor().stored_energy_j(model.charge_c())?;
        model.kinetic_energy()?;
        model.electrical_force_n()?;
        poll(&mut cancelled)?;
        Ok(CapacitorMotionCandidate { model, step: CapacitorMotionStep { native: record,
            electrical_energy_j: electrical, kinetic_energy_j: kinetic,
            electrical_supply_j: electrical_supply, mechanical_supply_j: mechanical_supply,
            field_work_on_body_j: field_work, resistors: heat, contacts: contact_heat, heat_j,
            energy_residual_j: total_residual } })
    }

    fn assemble(&self, dae: &CircuitDae)
        -> Result<(DescriptorPortHamiltonian, GapStorage), MaterialCircuitError>
    {
        let old_n = dae.system.state_dim();
        let old_d = dae.system.differential_dim();
        let n = old_n+2;
        let nd = old_d+2;
        let m = dae.ports+1;
        if n > super::MAX_STATES || m > super::MAX_PORTS {
            return Err(invalid("moving capacitor exceeds bounded descriptor capacity"));
        }
        let q_index = dae.charge_index[self.charge_slot];
        let old_q = dae.system.storage().exact_quadratic_matrix()
            .ok_or_else(|| invalid("gap coupling requires the original quadratic circuit storage"))?;
        if old_q.len() != old_n*old_n
            || old_q[q_index*old_n+q_index] != 1.0/self.capacitor().capacitance_f()
            || (0..old_n).any(|i| i != q_index &&
                (old_q[q_index*old_n+i] != 0.0 || old_q[i*old_n+q_index] != 0.0)) {
            return Err(invalid("active capacitor does not match the original circuit storage stamp"));
        }
        let mut q = vec![0.0; nd*nd];
        for i in 0..old_d {
            for k in 0..old_d { q[i*nd+k] = old_q[i*old_n+k]; }
        }
        // Remove only the active capacitor's original energy; all other C/L
        // terms remain exactly as stamped by the existing Kirchhoff owner.
        q[q_index*nd+q_index] = 0.0;
        let law = GapStorage { other: QuadraticStorage::new(q, nd).map_err(phs)?,
            q_index, gap_index: old_d, momentum_index: old_d+1,
            mass: self.body.inertia().mass(), eps_area: epsilon_area(self.capacitor())?,
            limits: self.gap_limits_m };
        let (old_j, old_r, old_g) = dae.system.structure();
        let (mut j, mut r, mut g) = (vec![0.0; n*n], vec![0.0; n*n], vec![0.0; n*m]);
        let map = |i: usize| if i < old_d { i } else { i+2 };
        for i in 0..old_n {
            for k in 0..old_n {
                j[map(i)*n+map(k)] = old_j[i*old_n+k];
                r[map(i)*n+map(k)] = old_r[i*old_n+k];
            }
            for p in 0..dae.ports { g[map(i)*m+p] = old_g[i*dae.ports+p]; }
        }
        j[old_d*n+old_d+1] = 1.0;
        j[(old_d+1)*n+old_d] = -1.0;
        g[(old_d+1)*m+dae.ports] = 1.0;
        let system = DescriptorPortHamiltonian::new(n, nd, m, j, r, g,
            Box::new(law.clone())).map_err(phs)?;
        Ok((system, law))
    }
}

#[derive(Clone)]
struct GapStorage {
    other: QuadraticStorage,
    q_index: usize,
    gap_index: usize,
    momentum_index: usize,
    mass: f64,
    eps_area: f64,
    limits: [f64; 2],
}
impl GapStorage {
    fn electrical_energy(&self, x: &[f64]) -> f64 {
        let charge = x[self.q_index];
        self.other.hamiltonian(x) + 0.5*charge*(charge/(self.eps_area/x[self.gap_index]))
    }
    fn kinetic_energy(&self, x: &[f64]) -> f64 {
        0.5*x[self.momentum_index]*(x[self.momentum_index]/self.mass)
    }
}
impl Storage for GapStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if !within(x[self.gap_index], self.limits) { return f64::NAN; }
        self.electrical_energy(x)+self.kinetic_energy(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if !within(x[self.gap_index], self.limits) { out.fill(f64::NAN); return; }
        self.other.gradient(x, out);
        let charge = x[self.q_index];
        out[self.q_index] += charge/(self.eps_area/x[self.gap_index]);
        out[self.gap_index] += 0.5*charge*(charge/self.eps_area);
        out[self.momentum_index] += x[self.momentum_index]/self.mass;
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], _: &[f64]) -> Option<f64> {
        if !within(a[self.gap_index], self.limits) || !within(b[self.gap_index], self.limits) {
            return Some(f64::NAN);
        }
        let dq = b[self.q_index]-a[self.q_index];
        let dd = b[self.gap_index]-a[self.gap_index];
        // Exact cubic midpoint remainder. The retained circuit and kinetic
        // terms are quadratic and contribute zero; avoid subtracting large Hs.
        Some(0.125*dq*(dq/self.eps_area)*dd)
    }
}
