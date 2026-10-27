//! Source-bound resistor heat on the existing Kirchhoff storage/port solver.
//!
//! Capacitors and inductors keep their original charges and fluxes while each
//! resistor is re-resolved by its material owner. No RC/RLC time integrator,
//! circuit topology solver, or second resistance law is introduced here.
//! Only differential storage is history. Node potentials and constraint currents
//! are algebraic efforts; they must not become independent energy stores.

/// Joint finite-supply, material-enthalpy, gas and piston evolution.
pub mod thermal;
/// Material/geometry-derived C and L, with explicit work on live reconfiguration.
pub mod storage;
/// Reciprocal electric/magnetic gap motion driven by the same field energy.
pub mod motion;
/// Self-consistent conductivity at solved electric field and current density.
pub mod conductivity;
/// Source-bound electrical contacts in series with original bulk branches.
pub mod contacts;
/// Source-derived piezoelectric motion on the same Kirchhoff storage solver.
pub mod piezoelectric;
use contacts::{SeriesElectricalContact, ContactHeat};

use fs_material::conductor::{ConductorError, OhmicDrive, ResolvedConductor};
use fs_phs::circuit::{Branch, CircuitDae, CircuitError, CircuitGraph, assemble_circuit};
use fs_phs::StepRecord;

/// Bounds on the existing allocating dense descriptor solve, not a sparse solver.
pub const MAX_MATERIAL_CIRCUIT_NODES: usize = 32;
/// Maximum branches, including source and storage branches.
pub const MAX_MATERIAL_CIRCUIT_BRANCHES: usize = 64;
const MAX_STATES: usize = 128;
const MAX_PORTS: usize = 32;
const MAX_TRANSFORMERS: usize = 16;

/// A complete source-bound resistor at one explicitly selected physical state.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialResistor {
    /// Index in the ORIGINAL `CircuitGraph::branches`, not resistor-only order.
    pub branch: usize,
    /// Original material bundle, current geometry and resistance owner.
    pub conductor: ResolvedConductor,
}

/// Absolute consistency tolerances; these do not bound time-discretization error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialCircuitSolve {
    /// Circuit storage/supply and branch-heat partition allowance [J].
    pub energy_tolerance_j: f64,
    /// Per-resistor voltage/current versus constitutive power allowance [W].
    pub power_tolerance_w: f64,
}

/// A refused proposal never changes the accepted circuit or any thermal store.
#[derive(Debug)]
pub enum MaterialCircuitError {
    /// Incomplete binding, invalid geometry/topology/state, or bounded work exceeded.
    Invalid(&'static str),
    /// Original Kirchhoff/descriptor owner refusal.
    Circuit(CircuitError),
    /// Original material-conductor evaluation refusal.
    Conductor(ConductorError),
    /// Source-derived electric/magnetic storage or field-energy refusal.
    Storage(fs_material::field_storage::FieldStorageError),
    /// Original source-bound mechanical section or layer-strain refusal.
    Material(fs_material::state_point::MaterialStatePointError),
    /// Original piezoelectric material, geometry or operating-domain refusal.
    Piezoelectric(crate::piezoelectric::PiezoelectricCouplingError),
    /// Independent branch heat or total electrical balance did not close.
    Balance(&'static str),
    /// The material response did not agree with the solved branch fields.
    ConductivityConvergence {
        /// Actual native circuit trials performed.
        iterations: usize,
        /// Largest absolute applied/evaluated resistance difference [ohm].
        resistance_defect_ohm: f64,
        /// Largest absolute constitutive branch-power difference [W].
        power_defect_w: f64,
    },
    /// Cancellation before publication.
    Cancelled,
}
impl From<CircuitError> for MaterialCircuitError {
    fn from(e: CircuitError) -> Self { Self::Circuit(e) }
}
impl From<ConductorError> for MaterialCircuitError {
    fn from(e: ConductorError) -> Self { Self::Conductor(e) }
}
impl From<fs_material::field_storage::FieldStorageError> for MaterialCircuitError {
    fn from(e: fs_material::field_storage::FieldStorageError) -> Self { Self::Storage(e) }
}
impl From<fs_material::state_point::MaterialStatePointError> for MaterialCircuitError {
    fn from(e: fs_material::state_point::MaterialStatePointError) -> Self { Self::Material(e) }
}
impl core::fmt::Display for MaterialCircuitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material circuit: {self:?}")
    }
}
impl std::error::Error for MaterialCircuitError {}

/// One accepted graph and its conserved electrical coordinates. C/L are held
/// during an ordinary circuit step; source-bound changes between steps use
/// [explicit parameter work](Self::rebind_storage_candidate). Source-port
/// mapping and transformer topology stay fixed. [`motion::MaterialFieldMotion`]
/// evolves a source-bound capacitor's or inductor's gap with its moving body.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialCircuit {
    graph: CircuitGraph,
    resistors: Vec<MaterialResistor>,
    contacts: Vec<SeriesElectricalContact>,
    fluxes_wb: Vec<f64>,
    charges_c: Vec<f64>,
    storage: Vec<storage::MaterialStorage>,
}

/// Branch quantities use the SAME algebraic midpoint as the descriptor's
/// dissipation. Averaging separately reconstructed endpoint voltages would
/// not implement this contract when sources or resistances change.
#[derive(Clone, Debug, PartialEq)]
pub struct ResistorHeat {
    /// Original source binding, retained without dropping any bundle properties.
    pub material: MaterialResistor,
    /// Signed voltage across the BULK conductor [V], excluding series contacts.
    /// With no contacts this is the original terminal-to-terminal voltage.
    pub voltage_v: f64,
    /// Signed first-to-second-terminal current [A].
    pub current_a: f64,
    /// Original conductor's nonnegative Joule power [W].
    pub power_w: f64,
    /// Heat delivered to this resistor during the physical step [J].
    pub heat_j: f64,
}

/// Circuit owner results plus an independently reconstructed material heat split.
pub struct MaterialCircuitStep {
    /// Original step, including external port work and solver residual.
    pub native: StepRecord,
    /// Original circuit Hamiltonian at accepted start [J].
    pub stored_before_j: f64,
    /// Same Hamiltonian at the proposed endpoint [J].
    pub stored_after_j: f64,
    /// Heat by original branch index, once per resistor (never by node degree).
    pub resistors: Vec<ResistorHeat>,
    /// Heat in each original electrical interface, in stable contact-ID order.
    pub contacts: Vec<ContactHeat>,
    /// Sum of bulk AND contact heat, checked against the original native loss.
    pub heat_j: f64,
    /// Actual storage change plus branch heat minus source work [J].
    pub energy_residual_j: f64,
}
/// Electrical endpoint and its original per-branch source/heat records. Consumers
/// must publish this with ALL receiving thermal candidates, not debit it in a
/// trial callback. Neither constructing nor inspecting a candidate mutates input.
pub struct MaterialCircuitCandidate {
    circuit: MaterialCircuit,
    step: MaterialCircuitStep,
}
impl MaterialCircuitCandidate {
    /// Proposed electrical history; no thermal store has been updated.
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    /// Actual source work, electrical storage and branch-local heat.
    pub const fn step(&self) -> &MaterialCircuitStep { &self.step }
    /// Move the electrical endpoint and accounting into a coupled transaction.
    pub fn into_parts(self) -> (MaterialCircuit, MaterialCircuitStep) { (self.circuit, self.step) }
}

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialCircuitError> {
    if cancelled() { Err(MaterialCircuitError::Cancelled) } else { Ok(()) }
}
/// Shared bounded graph/storage admission for lumped and spatial conductors.
/// Returns resistor count; the native circuit owner still admits all topology.
pub(crate) fn admit_graph(graph: &CircuitGraph) -> Result<usize, MaterialCircuitError> {
    use MaterialCircuitError::Invalid;
    if graph.node_count > MAX_MATERIAL_CIRCUIT_NODES
        || graph.branches.len() > MAX_MATERIAL_CIRCUIT_BRANCHES
        || graph.transformers.len() > MAX_TRANSFORMERS {
        return Err(Invalid("circuit exceeds bounded dense graph capacity"));
    }
    let mut n_l = 0; let mut n_c = 0; let mut n_v = 0;
    let mut n_r = 0;
    for &(a, b, branch) in &graph.branches {
        // A self-loop has no physical terminal voltage; do not stamp a fictitious
        // shunt loss (including the legacy owner's resistor self-loop case).
        if a == b { return Err(Invalid("circuit branches require distinct terminals")); }
        match branch {
            Branch::Inductor { henries } => {
                if !positive(henries) || !positive(1.0/henries) { return Err(Invalid("unrepresentable inductive storage")); }
                n_l += 1;
            }
            Branch::Capacitor { farads } => {
                if !positive(farads) || !positive(1.0/farads) { return Err(Invalid("unrepresentable capacitive storage")); }
                n_c += 1;
            }
            Branch::Resistor { .. } => n_r += 1,
            Branch::VoltageSource { port } | Branch::CurrentSource { port } => {
                if port >= MAX_PORTS { return Err(Invalid("circuit source port exceeds capacity")); }
                if matches!(branch, Branch::VoltageSource { .. }) { n_v += 1; }
            }
        }
    }
    if n_l + 2*n_c + graph.node_count.saturating_sub(1) + n_v + graph.transformers.len() > MAX_STATES {
        return Err(Invalid("circuit exceeds bounded descriptor state capacity"));
    }
    Ok(n_r)
}

fn bind(graph: &CircuitGraph, mut sources: Vec<MaterialResistor>, contacts: &[SeriesElectricalContact])
    -> Result<(CircuitGraph, Vec<MaterialResistor>, CircuitDae), MaterialCircuitError>
{
    use MaterialCircuitError::Invalid;
    let n_r = admit_graph(graph)?;
    if sources.len() != n_r { return Err(Invalid("every resistor requires exactly one material binding")); }
    sources.sort_by_key(|s| s.branch);
    let mut lowered = graph.clone();
    for (i, source) in sources.iter().enumerate() {
        if i > 0 && sources[i-1].branch == source.branch { return Err(Invalid("duplicate resistor binding")); }
        let branch = lowered.branches.get_mut(source.branch).ok_or(Invalid("material branch index out of range"))?;
        if !matches!(branch.2, Branch::Resistor { .. }) { return Err(Invalid("material binding does not name a resistor")); }
        let r = source.conductor.resistance_ohm();
        if !positive(r) || !positive(1.0/r) { return Err(Invalid("resistor conductance is unrepresentable")); }
        branch.2 = Branch::Resistor { ohms: r };
    }
    contacts::lower(&mut lowered, contacts)?;
    let dae = assemble_circuit(&lowered)?;
    Ok((lowered, sources, dae))
}
pub(crate) fn extended(dae: &CircuitDae, fluxes: &[f64], charges: &[f64], inputs: &[f64])
    -> Result<Vec<f64>, MaterialCircuitError>
{
    if fluxes.len() != dae.flux_index.len() || charges.len() != dae.charge_index.len()
        || inputs.len() != dae.ports || fluxes.iter().chain(charges).chain(inputs).any(|v| !v.is_finite()) {
        return Err(MaterialCircuitError::Invalid("finite complete branch-ordered fluxes, charges and ports required"));
    }
    let mut x = vec![0.0; dae.system.state_dim()];
    for (&index, &value) in dae.flux_index.iter().zip(fluxes).chain(dae.charge_index.iter().zip(charges)) {
        x[index] = value;
    }
    Ok(x)
}

// Shared native-effort accounting for electrical and coupled mechanical storage.
// Bulk and interface losses retain their original owners and separate receivers.
#[allow(clippy::too_many_arguments)]
fn account_step(
    graph: &CircuitGraph, resistors: &[MaterialResistor], contacts: &[SeriesElectricalContact],
    node_indices: &[usize], x0: &[f64], record: &StepRecord, h0: f64, h1: f64, dt: f64,
    controls: MaterialCircuitSolve, cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<ResistorHeat>, Vec<ContactHeat>, f64, f64), MaterialCircuitError> {
    if record.x.iter().chain(&record.y).any(|v| !v.is_finite())
        || [h0, h1, record.delta_h, record.dissipated, record.supplied, record.solver_residual]
            .iter().any(|v| !v.is_finite()) || h0 < 0.0 || h1 < 0.0 || record.dissipated < 0.0 {
        return Err(MaterialCircuitError::Balance("nonfinite or negative storage/loss"));
    }
    let voltage = |node: usize| {
        if node == 0 { 0.0 } else {
            let i = node_indices[node-1];
            f64::midpoint(x0[i], record.x[i])
        }
    };
    let mut heat_j = 0.0;
    let mut heat = Vec::with_capacity(resistors.len());
    let mut contact_heat = Vec::with_capacity(contacts.len());
    for source in resistors {
        poll(cancelled)?;
        let (a,b,branch) = graph.branches[source.branch];
        let terminal_v = voltage(a)-voltage(b);
        let has_contact = contacts.iter().any(|c| c.branch() == source.branch);
        // Preserve the exact legacy conductor arithmetic on contact-free branches.
        let (v,i) = if has_contact {
            let Branch::Resistor { ohms } = branch else { unreachable!("admitted resistor") };
            let current = terminal_v/ohms;
            if !current.is_finite() || (terminal_v != 0.0 && current == 0.0) {
                return Err(MaterialCircuitError::Balance("series current is unrepresentable"));
            }
            (current*source.conductor.resistance_ohm(),current)
        } else { (terminal_v,source.conductor.current_a_for_drive(OhmicDrive::Voltage(terminal_v))?) };
        let p = source.conductor.joule_power_w(i)?;
        let q = dt*p;
        let defect = v*i-p;
        if !q.is_finite() || (p != 0.0 && q == 0.0) || !defect.is_finite()
            || defect.abs() > controls.power_tolerance_w {
            return Err(MaterialCircuitError::Balance("resistor electrical work and material heat disagree"));
        }
        heat_j += q;
        let mut branch_power = p;
        for contact in contacts.iter().filter(|c| c.branch() == source.branch) {
            poll(cancelled)?;
            let loss = contact.heat(i, dt)?;
            branch_power += loss.power_w;
            heat_j += loss.heat_j;
            contact_heat.push(loss);
        }
        let terminal_defect = terminal_v*i-branch_power;
        if !terminal_defect.is_finite() || terminal_defect.abs() > controls.power_tolerance_w {
            return Err(MaterialCircuitError::Balance("bulk and contact powers do not close terminal work"));
        }
        heat.push(ResistorHeat { material:source.clone(), voltage_v:v, current_a:i, power_w:p, heat_j:q });
    }
    contact_heat.sort_by_key(|h| h.contact.id());
    let residual = (h1-h0)+heat_j-record.supplied;
    if !heat_j.is_finite() || !residual.is_finite()
        || (heat_j-record.dissipated).abs() > controls.energy_tolerance_j
        || ((h1-h0)-record.delta_h).abs() > controls.energy_tolerance_j
        || residual.abs() > controls.energy_tolerance_j {
        return Err(MaterialCircuitError::Balance("branch heat, storage and external supply do not close"));
    }
    Ok((heat, contact_heat, heat_j, residual))
}

impl MaterialCircuit {
    /// Bind ALL resistor branches and initialize flux/charge in their respective
    /// original branch order. Numeric resistor fields on `graph` are replaced,
    /// never an unsourced fallback. Initial algebraic consistency is checked by
    /// the existing descriptor owner; only the supplied electrical storage is kept.
    pub fn new(graph: CircuitGraph, sources: Vec<MaterialResistor>, fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>, initial_inputs: &[f64]) -> Result<Self, MaterialCircuitError>
    {
        Self::new_with_contacts(graph, sources, Vec::new(), fluxes_wb, charges_c, initial_inputs)
    }

    /// Bind bulk conductors and separately sourced contact drops on those same
    /// branches. Every contact is R''/area in SERIES, never replacement bulk
    /// resistivity. Contact heat remains distinct from bulk heat in the report.
    /// The existing descriptor owns initialization and all electrical storage.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_contacts(graph: CircuitGraph, sources: Vec<MaterialResistor>,
        mut contacts: Vec<SeriesElectricalContact>, fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>, initial_inputs: &[f64]) -> Result<Self, MaterialCircuitError>
    {
        contacts::admit(&mut contacts)?;
        let (graph, resistors, dae) = bind(&graph, sources, &contacts)?;
        let x = extended(&dae, &fluxes_wb, &charges_c, initial_inputs)?;
        let initial = dae.consistent_initial_state(&x, initial_inputs)?;
        let h = dae.system.hamiltonian(&x);
        if !h.is_finite() || h < 0.0 || initial.iter().any(|v| !v.is_finite()) {
            return Err(MaterialCircuitError::Invalid("unrepresentable initial circuit storage or constraints"));
        }
        for i in 0..dae.system.differential_dim() {
            let allowance = 64.0*f64::EPSILON*x[i].abs();
            if (initial[i]-x[i]).abs() > allowance {
                return Err(MaterialCircuitError::Invalid("algebraic initialization changed a conserved electrical coordinate"));
            }
        }
        Ok(Self { graph, resistors, contacts, fluxes_wb, charges_c, storage: Vec::new() })
    }
    /// Original graph topology with the last accepted material resistor values.
    pub const fn graph(&self) -> &CircuitGraph { &self.graph }
    /// Original source bindings at construction or the last accepted step midpoint.
    pub fn resistors(&self) -> &[MaterialResistor] { &self.resistors }
    /// Explicit electrical contact sources, in stable ID order.
    pub fn contacts(&self) -> &[SeriesElectricalContact] { &self.contacts }
    /// Conserved inductor flux linkage history [Wb], in inductor branch order.
    pub fn fluxes_wb(&self) -> &[f64] { &self.fluxes_wb }
    /// Conserved capacitor charge history [C], in capacitor branch order.
    pub fn charges_c(&self) -> &[f64] { &self.charges_c }

    /// Re-evaluate a candidate with caller-resolved midpoint resistors and held
    /// external source ports. Source-free graphs use `inputs = &[]`: all heat
    /// then comes from finite electrical storage. Each trial starts from the
    /// SAME charges/fluxes, so nonlinear thermal iteration cannot drain them twice.
    ///
    /// The original descriptor owner solves KCL/KVL and electrical storage.
    /// Its algebraic effort is `(x0_alg+x1_alg)/2`, not the static Hamiltonian
    /// gradient (which is zero there). Resistor heat is reconstructed from those
    /// exact terminal efforts through the original material conductor law.
    /// Per-branch work, summed dissipation and total electrical energy must agree.
    ///
    /// Coefficients C/L and topology are frozen. Source waveforms and resistance
    /// variation are midpoint approximations, not switching/arc, dielectric loss,
    /// induction-geometry feedback, or a time-accuracy certificate. The dense
    /// solver has its own bounded Newton work; cancellation is checked between
    /// assembly, stepping and at most 64 branch evaluations and before return.
    pub fn step_candidate(&self, sources: Vec<MaterialResistor>, inputs: &[f64], dt: f64,
        controls: MaterialCircuitSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<MaterialCircuitCandidate, MaterialCircuitError>
    {
        poll(&mut cancelled)?;
        if !self.contacts.is_empty() {
            return Err(MaterialCircuitError::Invalid("contact-bound circuits require explicit refreshed contact sources"));
        }
        self.step_with_contacts_candidate(sources, Vec::new(), inputs, dt, controls, cancelled)
    }

    /// Advance through the SAME circuit solver, retaining individual interface
    /// drops and losses as well as bulk Joule heat. Contacts are re-evaluated by
    /// the caller at the same state as the bulk material. IDs/branch assignments
    /// cannot change during a step; explicit reconfiguration is a separate call.
    /// Source and geometry swaps do not reset charge, flux, or storage provenance.
    /// Nonzero contact area and zero interfacial energy storage are assumptions;
    /// opening a circuit, arcing and dynamic interface memory are not modeled.
    pub fn step_with_contacts_candidate(&self, sources: Vec<MaterialResistor>,
        mut contacts: Vec<SeriesElectricalContact>, inputs: &[f64], dt: f64,
        controls: MaterialCircuitSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<MaterialCircuitCandidate, MaterialCircuitError>
    {
        poll(&mut cancelled)?;
        contacts::admit(&mut contacts)?;
        if contacts.len() != self.contacts.len() || contacts.iter().zip(&self.contacts)
            .any(|(a,b)| (a.id(),a.branch()) != (b.id(),b.branch())) {
            return Err(MaterialCircuitError::Invalid("electrical contact identities changed inside a physical step"));
        }
        if !positive(dt) || !positive(controls.energy_tolerance_j) || !positive(controls.power_tolerance_w) {
            return Err(MaterialCircuitError::Invalid("positive finite clock and absolute work/power tolerances required"));
        }
        let (graph, resistors, dae) = bind(&self.graph, sources, &contacts)?;
        let x0 = extended(&dae, &self.fluxes_wb, &self.charges_c, inputs)?;
        poll(&mut cancelled)?;
        let (record, _) = dae.step_audited(&x0, inputs, dt)?;
        let h0 = dae.system.hamiltonian(&x0);
        let h1 = dae.system.hamiltonian(&record.x);
        let (heat, contact_heat, heat_j, residual) = account_step(&graph, &resistors,
            &contacts, &dae.node_potential_index, &x0, &record, h0, h1, dt, controls, &mut cancelled)?;
        let next = Self { graph, resistors, contacts,
            fluxes_wb:dae.flux_index.iter().map(|&i| record.x[i]).collect(),
            charges_c:dae.charge_index.iter().map(|&i| record.x[i]).collect(),
            storage:self.storage.clone() };
        poll(&mut cancelled)?;
        Ok(MaterialCircuitCandidate { circuit:next, step:MaterialCircuitStep {
            native:record, stored_before_j:h0, stored_after_j:h1, resistors:heat,
            contacts:contact_heat, heat_j, energy_residual_j:residual } })
    }
}
