//! Reciprocal electric/magnetic field and rigid-body motion on the circuit DAE.
//!
//! A uniform, constant-area dielectric gap stores `q^2 d / (2 epsilon A)`.
//! A uniform magnetic gap stores `lambda^2 d / (2 mu N^2 A)`, where lambda is
//! flux LINKAGE, not magnetic flux. Terminal effort and attractive force are
//! derivatives of that ONE energy with respect to charge/flux linkage and gap.
//! The same descriptor step evolves the circuit and guided moving body; motion
//! is not an externally imposed C/L update or an independently integrated force.
//! Body mass comes from the existing regional material/inertia compiler.
//! Magnetic geometry is a single uniform gap with a negligible-reluctance return
//! path; it does not describe a solenoid changing the length of its entire coil.

/// Source-bound layered axial support, on the same gap and momentum coordinates.
pub mod support;
/// Material conduction at the operating fields solved with actual gap motion.
pub mod conductivity;
use support::AxialSupport;

use fs_matdb::EvaluationDecision;
use fs_material::field_storage::{ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_PROPERTY,
    ResolvedCapacitor, ResolvedInductor};
use fs_mbd_material::ResolvedBody;
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, StepRecord, Storage, step_descriptor};
use fs_phs::circuit::{Branch, CircuitDae, CircuitError};

use super::storage::MaterialStorageElement;
use super::contacts::{self, ContactHeat, SeriesElectricalContact};
use super::{MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
    ResistorHeat, account_step, bind_with_mutual, extended, poll, positive};

/// A source-bound capacitor OR inductor and one translationally guided body.
/// Positive gap velocity separates the electrodes or magnetic poles; the other
/// boundary is fixed. Regional density/geometry, field material, remaining
/// circuit storage and graph topology stay attached to this accepted state.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialFieldMotion {
    circuit: MaterialCircuit,
    body: ResolvedBody,
    support: Option<AxialSupport>,
    storage_slot: usize,
    coordinate_slot: usize,
    gap_limits_m: [f64; 2],
    momentum_kg_m_s: f64,
}

/// Held external inputs over one common electrical/mechanical interval.
#[derive(Clone, Copy, Debug)]
pub struct FieldMotionInput<'a> {
    /// Original circuit port order, unchanged by the added mechanical coordinate.
    pub electrical_ports: &'a [f64],
    /// Force applied by an external owner, positive toward increasing gap [N].
    pub external_force_n: f64,
    /// Common physical duration [s].
    pub duration_s: f64,
}

/// Joint field, mechanical, source and heat accounting from one descriptor step.
pub struct FieldMotionStep {
    /// Native state order: original fluxes/charges, gap, momentum, then original
    /// algebraic coordinates. Outputs preserve circuit ports and append the
    /// discrete-gradient velocity, which is conjugate to external force.
    pub native: StepRecord,
    /// Total capacitor/inductor energy before and after the coupled step [J].
    pub electrical_energy_j: [f64; 2],
    /// Body plus optional axial member's consistent kinetic energy [J].
    pub kinetic_energy_j: [f64; 2],
    /// Original layered support's elastic energy before and after the step [J].
    /// Zero without a support; this reversible store is never deposited as heat.
    pub elastic_energy_j: [f64; 2],
    /// Work supplied through original electrical ports [J].
    pub electrical_supply_j: f64,
    /// External mechanical work, force times actual gap change [J].
    pub mechanical_supply_j: f64,
    /// Work transferred from electrical storage to body/support motion [J]. This is internal
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
pub struct FieldMotionCandidate {
    model: MaterialFieldMotion,
    step: FieldMotionStep,
}
impl FieldMotionCandidate {
    /// Joint endpoint; electrical charge cannot be accepted without the same gap.
    pub const fn model(&self) -> &MaterialFieldMotion { &self.model }
    /// Original solver result and independently partitioned physical work.
    pub const fn step(&self) -> &FieldMotionStep { &self.step }
    /// Publish the joint model and its work/heat record together.
    pub fn into_parts(self) -> (MaterialFieldMotion, FieldMotionStep) {
        (self.model, self.step)
    }
}

fn invalid(what: &'static str) -> MaterialCircuitError { MaterialCircuitError::Invalid(what) }
fn phs(error: fs_phs::PhsError) -> MaterialCircuitError { CircuitError::from(error).into() }
fn within(gap: f64, limits: [f64; 2]) -> bool {
    gap.is_finite() && gap > limits[0] && gap < limits[1]
}
fn coupling_scale(element: &MaterialStorageElement) -> Result<f64, MaterialCircuitError> {
    let (material, property_name, area, turns) = match element {
        MaterialStorageElement::Capacitor(capacitor) =>
            (capacitor.material(), ELECTRIC_PERMITTIVITY_PROPERTY, capacitor.area_m2(), 1.0),
        MaterialStorageElement::Inductor(inductor) =>
            (inductor.material(), MAGNETIC_PERMEABILITY_PROPERTY, inductor.area_m2(), f64::from(inductor.turns())),
    };
    let property = material.property(property_name)
        .ok_or_else(|| invalid("moving field element lost its material source"))?;
    // A temperature-dependent free-energy coefficient does not by itself give
    // internal energy: the missing field entropy would need a thermal owner.
    // Do not label that missing reversible heat as conserved mechanical energy.
    if property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity {
        return Err(invalid("conservative gap motion requires a validity-wide constant field coefficient"));
    }
    let value = property.value_si() * area * turns * turns;
    if !positive(value) || !positive(1.0 / value) {
        return Err(invalid("field coupling scale is unrepresentable"));
    }
    Ok(value)
}

impl MaterialFieldMotion {
    /// Bind an existing material RLC circuit to a source-resolved rigid body.
    /// `branch` is the ORIGINAL graph index of its moving capacitor or inductor.
    /// The capacitor gap or inductor path length is the initial separation;
    /// mass is ONLY `body.inertia().mass()`.
    ///
    /// The caller supplies two positive, finite, ordered OPEN gap limits from
    /// the admitted geometry model and an initial velocity [m/s]. Area stays
    /// fixed. Body regions translate together; rotation/deformation are guided
    /// out. Source and geometry swaps before construction change the actual
    /// dynamics without changing this code. Live body/material replacement needs
    /// its own mass, momentum and thermodynamic transaction.
    ///
    /// Electric geometry assumes ideal electrodes and a uniform dielectric gap.
    /// Magnetic geometry assumes ONE uniform, constant-area gap, with all flux
    /// linking the declared winding count and negligible return-path reluctance.
    /// The inductor's permeability is the GAP medium's permeability, not an
    /// effective ferromagnetic-core coefficient. The moving body's material
    /// determines its inertia; pole rigidity and ideal magnetic return geometry
    /// remain assembly assumptions. This is not a moving-length solenoid model.
    ///
    /// Both reductions are linear and nondispersive with a fixed opposite face.
    /// No fringing, leakage, saturation, hysteresis, pole/plate contact, breakdown,
    /// magnetostriction, electrostriction or thermal field law is supplied.
    /// Temperature-dependent coefficients refuse because their reversible heat
    /// requires field entropy data. Current-dependent coil loss remains with
    /// explicitly bound resistor/contact branches, not the lossless gap energy.
    pub fn new(circuit: MaterialCircuit, branch: usize, body: ResolvedBody,
        velocity_m_s: f64, gap_limits_m: [f64; 2]) -> Result<Self, MaterialCircuitError>
    {
        Self::build(circuit, branch, body, None, velocity_m_s, gap_limits_m)
    }

    /// Couple a layered return member to the SAME electrical/mechanical step.
    /// Its original section owns elastic force and stored energy; the member's
    /// consistent axial inertia is added to the canonical guided-body mass.
    /// Source/geometry substitutions therefore change both dynamics and the
    /// electrical back-reaction without a separate force or spring integrator.
    ///
    /// All original field restrictions apply. The member has a fixed far end and
    /// a prescribed linear axial displacement shape. Ordinary supports retain
    /// frozen material state. Thermally live supports require the coupled
    /// enthalpy owner; no support damping, buckling or elastic waves are added.
    /// Construct the complete assembly before advancing history: attaching a
    /// member to a live state would require an explicit energy/momentum transfer.
    pub fn new_with_support(circuit: MaterialCircuit, branch: usize, body: ResolvedBody,
        support: AxialSupport, velocity_m_s: f64, gap_limits_m: [f64; 2])
        -> Result<Self, MaterialCircuitError>
    {
        Self::build(circuit, branch, body, Some(support), velocity_m_s, gap_limits_m)
    }

    fn build(circuit: MaterialCircuit, branch: usize, body: ResolvedBody,
        support: Option<AxialSupport>, velocity_m_s: f64, gap_limits_m: [f64; 2])
        -> Result<Self, MaterialCircuitError>
    {
        if circuit.mutual.iter().any(|link| link.henries != 0.0
            && (link.branches.0 == branch || link.branches.1 == branch)) {
            return Err(invalid("moving a mutually coupled winding requires a shared geometry-dependent inductance law"));
        }
        let slot = circuit.storage.iter().position(|source| source.branch == branch)
            .ok_or_else(|| invalid("moving field requires an original material storage binding"))?;
        let element = &circuit.storage[slot].element;
        let gap = match element {
            MaterialStorageElement::Capacitor(capacitor) => capacitor.gap_m(),
            MaterialStorageElement::Inductor(inductor) => inductor.path_length_m(),
        };
        if gap_limits_m.iter().any(|&v| !positive(v))
            || !within(gap, gap_limits_m) || !velocity_m_s.is_finite() {
            return Err(invalid("finite velocity and positive open gap limits containing the field gap required"));
        }
        let scale = coupling_scale(element)?;
        if gap_limits_m.iter().any(|gap| !positive(scale / gap) || !positive(gap / scale)) {
            return Err(invalid("field storage is unrepresentable within the declared gap domain"));
        }
        if let Some(support) = &support { support.force_and_energy(gap)?; }
        let mass = body.inertia().mass() + support.as_ref().map_or(0.0, AxialSupport::modal_mass_kg);
        let momentum = mass * velocity_m_s;
        if !positive(mass) || !positive(1.0 / mass) || !momentum.is_finite()
            || (velocity_m_s != 0.0 && momentum == 0.0) {
            return Err(invalid("source-derived body mass or momentum is unrepresentable"));
        }
        let coordinate_slot = circuit.graph.branches[..branch].iter().filter(|(_, _, b)| {
            matches!((element, b),
                (MaterialStorageElement::Capacitor(_), Branch::Capacitor { .. })
                | (MaterialStorageElement::Inductor(_), Branch::Inductor { .. }))
        }).count();
        let result = Self { circuit, body, support, storage_slot: slot, coordinate_slot,
            gap_limits_m, momentum_kg_m_s: momentum };
        result.field_energy()?;
        result.kinetic_energy()?;
        result.field_force_n()?;
        Ok(result)
    }

    /// Existing circuit at the SAME accepted field gap and electrical coordinates.
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    /// Canonical regional geometry, material receipts and spatial inertia.
    pub const fn body(&self) -> &ResolvedBody { &self.body }
    /// Optional member, retaining every original section layer and source path.
    pub const fn support(&self) -> Option<&AxialSupport> { self.support.as_ref() }
    /// Body mass plus the member's consistent axial inertia [kg].
    pub fn generalized_mass_kg(&self) -> f64 {
        self.body.inertia().mass() + self.support.as_ref().map_or(0.0, AxialSupport::modal_mass_kg)
    }
    /// Stored elastic energy of the original support at the accepted gap [J].
    pub fn elastic_energy_j(&self) -> Result<f64, MaterialCircuitError> {
        match &self.support {
            Some(support) => Ok(support.force_and_energy(self.gap_m())?.1),
            None => Ok(0.0),
        }
    }
    /// Original material snapshot and CURRENT geometry of the moving element.
    pub fn element(&self) -> &MaterialStorageElement {
        &self.circuit.storage[self.storage_slot].element
    }
    /// The moving dielectric gap, if the original branch is capacitive.
    pub fn capacitor(&self) -> Option<&ResolvedCapacitor> {
        match self.element() {
            MaterialStorageElement::Capacitor(capacitor) => Some(capacitor),
            MaterialStorageElement::Inductor(_) => None,
        }
    }
    /// The moving magnetic gap, if the original branch is inductive.
    pub fn inductor(&self) -> Option<&ResolvedInductor> {
        match self.element() {
            MaterialStorageElement::Inductor(inductor) => Some(inductor),
            MaterialStorageElement::Capacitor(_) => None,
        }
    }
    /// Positive electrode/pole separation [m], held in the field geometry owner.
    pub fn gap_m(&self) -> f64 {
        match self.element() {
            MaterialStorageElement::Capacitor(capacitor) => capacitor.gap_m(),
            MaterialStorageElement::Inductor(inductor) => inductor.path_length_m(),
        }
    }
    /// Generalized momentum of body and optional axial member [kg m/s].
    pub const fn momentum_kg_m_s(&self) -> f64 { self.momentum_kg_m_s }
    /// Physical endpoint velocity [m/s], positive toward increasing separation.
    pub fn velocity_m_s(&self) -> f64 { self.momentum_kg_m_s / self.generalized_mass_kg() }
    /// Instantaneous fixed-charge/flux-linkage force `-dH/dd` [N].
    /// This is `-q^2/(2 epsilon A)` or `-lambda^2/(2 mu N^2 A)`.
    /// Its discrete-gradient counterpart, not this endpoint sample, drives a step.
    pub fn field_force_n(&self) -> Result<f64, MaterialCircuitError> {
        let coordinate = self.field_coordinate();
        let force = -0.5 * coordinate * (coordinate / coupling_scale(self.element())?);
        if !force.is_finite() || (coordinate != 0.0 && force == 0.0) {
            return Err(invalid("field gap force is unrepresentable"));
        }
        Ok(force)
    }
    fn field_coordinate(&self) -> f64 {
        match self.element() {
            MaterialStorageElement::Capacitor(_) => self.circuit.charges_c[self.coordinate_slot],
            MaterialStorageElement::Inductor(_) => self.circuit.fluxes_wb[self.coordinate_slot],
        }
    }
    fn field_energy(&self) -> Result<f64, MaterialCircuitError> {
        Ok(match self.element() {
            MaterialStorageElement::Capacitor(capacitor) => capacitor.stored_energy_j(self.field_coordinate())?,
            MaterialStorageElement::Inductor(inductor) => inductor.stored_energy_j(self.field_coordinate())?,
        })
    }
    fn kinetic_energy(&self) -> Result<f64, MaterialCircuitError> {
        let energy = 0.5 * self.momentum_kg_m_s * self.velocity_m_s();
        if !energy.is_finite() || (self.momentum_kg_m_s != 0.0 && energy == 0.0) {
            return Err(invalid("body kinetic energy is unrepresentable"));
        }
        Ok(energy)
    }

    /// Solve electrical storage, field gap/momentum and KCL/KVL together using
    /// the original descriptor integrator. Resistor material bindings are held
    /// for this trial, as on `MaterialCircuit::step_candidate`; their heat uses
    /// the same native terminal efforts and conductor law.
    ///
    /// The moving element's quadratic term is REPLACED by its gap-dependent
    /// field energy: `q^2 d/(2 eps A)` or `lambda^2 d/(2 mu N^2 A)`.
    /// The original circuit's J/R/G stamps, transformer constraints and remaining
    /// storage are retained. An optional layered support adds its original
    /// material energy and consistent inertia. No force gain, separate circuit
    /// integrator, authored spring constant or overridden mass is introduced.
    ///
    /// The endpoint must also admit its unchanged proposed charge/flux under the
    /// original algebraic constraints at the NEW geometry. Incompatible ideal
    /// clamps refuse; no projection, switching impulse or lost charge is hidden.
    /// Out-of-domain Newton trials, gap limits, imbalance and cancellation leave
    /// the whole model unchanged. Work tolerances do not certify time accuracy.
    pub fn step_candidate(&self, sources: Vec<MaterialResistor>, input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<FieldMotionCandidate, MaterialCircuitError>
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
        contacts: Vec<SeriesElectricalContact>, input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<FieldMotionCandidate, MaterialCircuitError>
    {
        poll(&mut cancelled)?;
        if self.support.as_ref().is_some_and(AxialSupport::is_thermally_live) {
            return Err(invalid("live thermoelastic support requires the joint enthalpy and reversible-heat owner"));
        }
        self.step_bound_support_candidate(sources, contacts, input, controls, cancelled)
    }

    // Only the coupled enthalpy owner may change this member's thermal state.
    // The actual native descriptor, physical history and all branch accounting
    // are identical to the conservative path; no mechanics are solved twice.
    pub(crate) fn step_thermal_support_candidate(&self,
        section: fs_material::section::thermal::ConcentricThermoelasticSection,
        sources: Vec<MaterialResistor>, input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<FieldMotionCandidate, MaterialCircuitError>
    {
        poll(&mut cancelled)?;
        if !self.circuit.contacts.is_empty() {
            return Err(invalid("thermal support currently requires explicit bulk-only electrical heat routes"));
        }
        let support = self.support.as_ref()
            .ok_or_else(|| invalid("thermal support requires an existing material member"))?
            .with_thermal_section(section)?;
        let mut trial = self.clone();
        trial.support = Some(support);
        trial.elastic_energy_j()?;
        trial.step_bound_support_candidate(sources, Vec::new(), input, controls, cancelled)
    }

    fn step_bound_support_candidate(&self, sources: Vec<MaterialResistor>,
        mut contacts: Vec<SeriesElectricalContact>, input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<FieldMotionCandidate, MaterialCircuitError>
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
        let (graph, resistors, dae) = bind_with_mutual(&self.circuit.graph, sources,
            &contacts, &self.circuit.mutual)?;
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
        let elastic = [law.elastic_energy(&x0), law.elastic_energy(&record.x)];
        let nodes: Vec<_> = dae.node_potential_index.iter().map(|&i| map(i)).collect();
        let (heat, contact_heat, heat_j, residual) = account_step(&graph, &resistors, &contacts, &nodes, &x0, &record,
            electrical[0]+kinetic[0]+elastic[0], electrical[1]+kinetic[1]+elastic[1], dt, controls, &mut cancelled)?;
        let electrical_supply = dt * input.electrical_ports.iter().zip(&record.y)
            .map(|(&u, &y)| u*y).sum::<f64>();
        let mechanical_supply = input.external_force_n * (gap-self.gap_m());
        let native_mechanical_supply = dt * input.external_force_n * record.y[dae.ports];
        let field_work = kinetic[1]-kinetic[0]+elastic[1]-elastic[0]-mechanical_supply;
        let electrical_work = electrical_supply-heat_j-(electrical[1]-electrical[0]);
        let total_residual = (electrical[1]-electrical[0])+(kinetic[1]-kinetic[0])
            +(elastic[1]-elastic[0])+heat_j-electrical_supply-mechanical_supply;
        if [electrical_supply, mechanical_supply, native_mechanical_supply, field_work,
            electrical_work, total_residual].iter().any(|v| !v.is_finite())
            || (mechanical_supply-native_mechanical_supply).abs() > controls.energy_tolerance_j
            || (field_work-electrical_work).abs() > controls.energy_tolerance_j
            || (total_residual-residual).abs() > controls.energy_tolerance_j
            || total_residual.abs() > controls.energy_tolerance_j {
            return Err(MaterialCircuitError::Balance("reciprocal field work, displacement work and branch heat disagree"));
        }
        let mut storage = self.circuit.storage.clone();
        storage[self.storage_slot].element = match self.element() {
            MaterialStorageElement::Capacitor(capacitor) => MaterialStorageElement::Capacitor(
                ResolvedCapacitor::try_from_material(capacitor.material().clone(), capacitor.area_m2(), gap)?),
            MaterialStorageElement::Inductor(inductor) => MaterialStorageElement::Inductor(
                ResolvedInductor::try_from_material(inductor.material().clone(), inductor.area_m2(), gap, inductor.turns())?),
        };
        poll(&mut cancelled)?;
        let graph = super::storage::lower(&graph, &mut storage)?;
        let mut circuit = MaterialCircuit::new_with_contacts_and_mutual_inductance(graph, resistors,
            contacts, self.circuit.mutual.clone(),
            dae.flux_index.iter().map(|&i| record.x[i]).collect(),
            dae.charge_index.iter().map(|&i| record.x[i]).collect(), input.electrical_ports)?;
        circuit.storage = storage;
        let model = Self { circuit, body: self.body.clone(), support: self.support.clone(), storage_slot: self.storage_slot,
            coordinate_slot: self.coordinate_slot, gap_limits_m: self.gap_limits_m,
            momentum_kg_m_s: record.x[nd+1] };
        model.field_energy()?;
        model.kinetic_energy()?;
        model.field_force_n()?;
        model.elastic_energy_j()?;
        poll(&mut cancelled)?;
        Ok(FieldMotionCandidate { model, step: FieldMotionStep { native: record,
            electrical_energy_j: electrical, kinetic_energy_j: kinetic, elastic_energy_j: elastic,
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
            return Err(invalid("moving field element exceeds bounded descriptor capacity"));
        }
        let (coordinate_index, coefficient) = match self.element() {
            MaterialStorageElement::Capacitor(capacitor) =>
                (dae.charge_index[self.coordinate_slot], capacitor.capacitance_f()),
            MaterialStorageElement::Inductor(inductor) =>
                (dae.flux_index[self.coordinate_slot], inductor.inductance_h()),
        };
        let old_q = dae.system.storage().exact_quadratic_matrix()
            .ok_or_else(|| invalid("gap coupling requires the original quadratic circuit storage"))?;
        if old_q.len() != old_n*old_n
            // Scaled preparation of unrelated mutual windings may round the
            // independent diagonal by an ULP. Off-diagonal coupling to the
            // moving coordinate still must be exactly absent.
            || (old_q[coordinate_index*old_n+coordinate_index]*coefficient-1.0).abs()
                > 64.0*f64::EPSILON
            || (0..old_n).any(|i| i != coordinate_index &&
                (old_q[coordinate_index*old_n+i] != 0.0 || old_q[i*old_n+coordinate_index] != 0.0)) {
            return Err(invalid("active field element does not match the original circuit storage stamp"));
        }
        let mut q = vec![0.0; nd*nd];
        for i in 0..old_d {
            for k in 0..old_d { q[i*nd+k] = old_q[i*old_n+k]; }
        }
        // Remove only the active field element's original energy; all other C/L
        // terms remain exactly as stamped by the existing Kirchhoff owner.
        q[coordinate_index*nd+coordinate_index] = 0.0;
        let law = GapStorage { other: QuadraticStorage::new(q, nd).map_err(phs)?,
            coordinate_index, gap_index: old_d, momentum_index: old_d+1,
            mass: self.generalized_mass_kg(), scale: coupling_scale(self.element())?,
            limits: self.gap_limits_m, support: self.support.clone() };
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
    coordinate_index: usize,
    gap_index: usize,
    momentum_index: usize,
    mass: f64,
    scale: f64,
    limits: [f64; 2],
    support: Option<AxialSupport>,
}
impl GapStorage {
    fn electrical_energy(&self, x: &[f64]) -> f64 {
        let coordinate = x[self.coordinate_index];
        self.other.hamiltonian(x) + 0.5*coordinate*(coordinate/(self.scale/x[self.gap_index]))
    }
    fn kinetic_energy(&self, x: &[f64]) -> f64 {
        0.5*x[self.momentum_index]*(x[self.momentum_index]/self.mass)
    }
    fn elastic_energy(&self, x: &[f64]) -> f64 {
        self.support.as_ref().map_or(0.0, |support| support.force_and_energy(x[self.gap_index])
            .map_or(f64::NAN, |(_, energy)| energy))
    }
}
impl Storage for GapStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if !within(x[self.gap_index], self.limits) { return f64::NAN; }
        self.electrical_energy(x)+self.kinetic_energy(x)+self.elastic_energy(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if !within(x[self.gap_index], self.limits) { out.fill(f64::NAN); return; }
        self.other.gradient(x, out);
        let coordinate = x[self.coordinate_index];
        out[self.coordinate_index] += coordinate/(self.scale/x[self.gap_index]);
        out[self.gap_index] += 0.5*coordinate*(coordinate/self.scale);
        out[self.momentum_index] += x[self.momentum_index]/self.mass;
        if let Some(support) = &self.support {
            match support.force_and_energy(x[self.gap_index]) {
                Ok((force, _)) => out[self.gap_index] -= force,
                Err(_) => out.fill(f64::NAN),
            }
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], _: &[f64]) -> Option<f64> {
        if !within(a[self.gap_index], self.limits) || !within(b[self.gap_index], self.limits) {
            return Some(f64::NAN);
        }
        if !self.elastic_energy(a).is_finite() || !self.elastic_energy(b).is_finite() {
            return Some(f64::NAN);
        }
        let dq = b[self.coordinate_index]-a[self.coordinate_index];
        let dd = b[self.gap_index]-a[self.gap_index];
        // Exact cubic midpoint remainder. The retained circuit and kinetic
        // and frozen-prestrain support terms are quadratic and contribute zero;
        // avoid subtracting large Hs.
        Some(0.125*dq*(dq/self.scale)*dd)
    }
}
