//! Growing, resolved-deposition thermal networks (issue #1, bounded R1 rung).
//!
//! Cells are isothermal material parcels, not printer commands or a mandatory
//! mesh. Their mass and specific enthalpy are authoritative; equilibrium chart
//! density never changes reference mass. Caller-owned geometry and stable IDs
//! survive cooling and reheating. Future cells do not exist in the operator.
//!
//! This composes the existing nonlinear lumped enthalpy and finite-contact
//! solvers. Each interval applies environment/source backward Euler followed by
//! contacts in stable-ID order (first-order Lie splitting). Pairwise transfer
//! is debited and credited once. A small algebraic residual is NOT a time,
//! splitting, spatial, or material-error bound; refine time and cells separately.
//! All active transfers share a conservative aggregate Biot <= 0.1 gate.
//! That gate is not a proof of cell uniformity: the caller must justify the
//! parcel partition and internal diffusion regime for its actual shapes.
//!
//! Geometry grows by birth, but existing parcels are not remeshed or deformed.
//! No glass-transition, crystallization kinetics, stress, healing, adhesion,
//! fracture, sag, contact detection, slicer, or controller is inferred. Contact
//! resistance is thermal, never a mechanical bond-strength claim. Birth and
//! contact switches have no claimed derivative. Clone is an in-memory checkpoint;
//! there is no persistent checkpoint codec or cross-ISA reproducibility claim.

use std::collections::BTreeMap;
use std::fmt;

mod schedule;

use fs_exec::Cx;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};

use crate::lumped::{
    BiotGate, LUMPED_BIOT_CEILING, LumpedContactStepConfig, LumpedEnthalpyBody,
    LumpedEnthalpyMarchConfig, solve_lumped_contact_step, solve_lumped_enthalpy,
};
use crate::{ConductionError, SeriesThermalResistance};

macro_rules! identity {
    ($name:ident, $doc:literal) => {
        #[doc = $doc]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(pub u64);
    };
}
identity!(CellId, "Stable material-parcel identity, independent of array position.");
identity!(BodyId, "Stable part, bed, or support identity.");
identity!(SegmentId, "Stable source deposition-segment identity; refinement may share it.");
identity!(MaterialId, "Stable material-law identity, not a table index.");
identity!(InterfaceId, "Stable thermal-interface identity; not a mechanical bond.");
identity!(EventId, "Strictly increasing resolved-schedule identity for replay protection.");

/// Explicit constant transport over the entire supplied enthalpy chart.
/// Constants are declarations, not calibrated printer/material parameters.
#[derive(Debug, Clone, Copy)]
pub struct DepositionMaterial<'a> {
    /// The same constitutive chart accepted by spatial enthalpy conduction.
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    /// Positive conductivity, W/(m K).
    pub conductivity_w_m_k: f64,
    /// Diffuse-gray hemispherical emissivity in [0,1]. Zero disables radiation.
    pub emissivity: f64,
}

/// Reference geometry in one caller-declared, right-handed world frame.
#[derive(Debug, Clone, PartialEq)]
pub struct CellGeometry<G> {
    /// Located sample/centroid in metres.
    pub position_m: [f64; 3],
    /// Orthonormal right-handed material axes in world coordinates (row vectors).
    /// Axis zero is the deposition direction; it does not initialize stress.
    pub material_axes: [[f64; 3]; 3],
    /// Positive reference volume, m3. Reference density is incoming mass / volume.
    pub volume_m3: f64,
    /// Whole reference surface area, m2, BEFORE subtracting thermal contacts.
    pub surface_area_m2: f64,
    /// Application-owned shape, mesh handle, feature correspondence, etc.
    pub detail: G,
}

/// Piecewise-constant environment between resolved schedule events.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Environment {
    /// Convection-fluid absolute temperature, K.
    pub ambient_k: f64,
    /// Large-enclosure effective radiation temperature, K.
    pub radiation_k: f64,
    /// Nonnegative convection coefficient, W/(m2 K), on EXPOSED area only.
    pub convection_w_m2_k: f64,
    /// Signed internally deposited power, W; not incoming material enthalpy.
    pub source_w: f64,
}

/// Material delivered instantaneously at a resolved event time.
/// Refine the schedule to represent continuous extrusion; there is no inferred
/// motion timing, flow rate, bead-width law, or pre-birth thermal capacity.
#[derive(Debug, Clone)]
pub struct MaterialBirth<G> {
    pub id: CellId,
    pub body: BodyId,
    pub segment: SegmentId,
    pub material: MaterialId,
    /// Delivered mass, kg, strictly positive.
    pub mass_kg: f64,
    /// Incoming specific enthalpy, J/kg, using THIS material chart's datum.
    pub specific_enthalpy_j_kg: f64,
    pub geometry: CellGeometry<G>,
    pub environment: Environment,
}

/// Caller-resolved contact. Areas must be disjoint on each parcel; the owner
/// of geometry must establish that fact. This layer checks total area coverage,
/// not geometric intersection. A finite bed/support can be an ordinary cell.
#[derive(Debug, Clone)]
pub struct ThermalContact {
    pub id: InterfaceId,
    pub a: CellId,
    pub b: CellId,
    /// Actual shared area, m2; removed from BOTH parcels' exposed surfaces.
    pub area_m2: f64,
    /// Total resistance, K/W, including any declared internal spreading terms.
    /// An area change requires an explicitly recomputed resistance.
    pub resistance: SeriesThermalResistance,
}

#[derive(Debug, Clone)]
pub enum ContactChange {
    /// Create, update, or reopen a contact, retaining its identity and heat history.
    Set(ThermalContact),
    /// Disable heat transfer and re-expose both surfaces, retaining history.
    Open(InterfaceId),
}

/// All births, contact changes and environmental changes at one physical time.
/// Births are resolved together before contacts, so same-time neighbours work.
#[derive(Debug, Clone)]
pub struct DepositionEvent<G> {
    pub id: EventId,
    /// Absolute simulation time, seconds, finite and nonnegative.
    pub time_s: f64,
    pub births: Vec<MaterialBirth<G>>,
    pub contacts: Vec<ContactChange>,
    pub environments: Vec<(CellId, Environment)>,
}

/// Explicit retained-state and per-call work limits; no hidden full-field history.
#[derive(Debug, Clone, Copy)]
pub struct DepositionBudget {
    pub max_cells: usize,
    /// Includes opened contacts whose history is retained.
    pub max_contacts: usize,
    pub max_materials: usize,
    pub max_events: usize,
    pub max_steps: usize,
    /// Environment and pair-contact solves, each with the existing bounded root solve.
    pub max_operator_solves: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct DepositionStepConfig {
    pub maximum_step_s: f64,
    /// Per-operator algebraic energy allowance, J.
    pub solver_energy_tolerance_j: f64,
    /// Total stored minus supplied energy allowance for this call, J.
    pub energy_tolerance_j: f64,
}

/// Immutable observation suitable for subsequent temperature-dependent mechanics.
#[derive(Debug, Clone)]
pub struct DepositedCell<G> {
    /// Original arrival declaration, retained unchanged through boundary updates.
    pub birth: MaterialBirth<G>,
    pub birth_time_s: f64,
    /// Incoming reference temperature is retained even after reheating.
    pub reference_temperature_k: f64,
    pub phase: EquilibriumPhaseState,
    /// Current environment; the original one remains in `birth.environment`.
    pub environment: Environment,
}

/// Interface thermal history is deliberately distinct from healing or damage.
#[derive(Debug, Clone)]
pub struct ContactState {
    pub contact: ThermalContact,
    pub active: bool,
    pub first_contact_time_s: f64,
    /// Signed cumulative heat from A to B, J, across committed steps only.
    pub heat_a_to_b_j: f64,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct DepositionReport {
    pub events: usize,
    pub steps: usize,
    pub operator_solves: usize,
    pub incoming_mass_kg: f64,
    pub incoming_enthalpy_j: f64,
    pub boundary_heat_j: f64,
    pub source_heat_j: f64,
    pub stored_energy_change_j: f64,
    pub energy_residual_j: f64,
    /// Sum of absolute pairwise contact-law residuals, not splitting error.
    pub contact_residual_l1_j: f64,
    pub maximum_biot: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum DepositionError {
    Invalid(&'static str),
    Budget(&'static str),
    Conduction(ConductionError),
    Phase { cell: CellId, reason: String },
    Energy { residual_j: f64, tolerance_j: f64 },
}
impl fmt::Display for DepositionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "deposition refused: {self:?}")
    }
}
impl std::error::Error for DepositionError {}
impl From<ConductionError> for DepositionError {
    fn from(error: ConductionError) -> Self { Self::Conduction(error) }
}

/// Growing-domain owner. Mutating public operations stage a private clone and
/// publish only on success; failure/cancellation never ages or double-adds cells.
/// Stable BTreeMap iteration fixes operator order on the same toolchain/ISA.
#[derive(Debug, Clone)]
pub struct DepositionSimulation<'a, G> {
    time_s: f64,
    last_event: Option<EventId>,
    materials: BTreeMap<MaterialId, DepositionMaterial<'a>>,
    cells: BTreeMap<CellId, DepositedCell<G>>,
    contacts: BTreeMap<InterfaceId, ContactState>,
    budget: DepositionBudget,
}

impl<'a, G: Clone> DepositionSimulation<'a, G> {
    pub fn new(
        materials: Vec<(MaterialId, DepositionMaterial<'a>)>,
        budget: DepositionBudget,
    ) -> Result<Self, DepositionError> {
        if materials.len() > budget.max_materials { return Err(DepositionError::Budget("materials")); }
        let mut table = BTreeMap::new();
        for (id, material) in materials {
            positive(material.conductivity_w_m_k)?;
            if !material.emissivity.is_finite() || !(0.0..=1.0).contains(&material.emissivity) {
                return Err(DepositionError::Invalid("emissivity must be in [0,1]"));
            }
            if table.insert(id, material).is_some() { return Err(DepositionError::Invalid("duplicate material ID")); }
        }
        Ok(Self { time_s: 0.0, last_event: None, materials: table,
            cells: BTreeMap::new(), contacts: BTreeMap::new(), budget })
    }

    pub fn time_s(&self) -> f64 { self.time_s }
    pub fn cells(&self) -> &BTreeMap<CellId, DepositedCell<G>> { &self.cells }
    pub fn contacts(&self) -> &BTreeMap<InterfaceId, ContactState> { &self.contacts }
    pub fn total_mass_kg(&self) -> f64 { sum(self.cells.values().map(|c| c.birth.mass_kg)) }
    pub fn stored_enthalpy_j(&self) -> f64 {
        sum(self.cells.values().map(|c| c.birth.mass_kg * c.phase.specific_enthalpy_j_kg()))
    }

    /// Exposed areas in stable cell-ID order. No inactive contact masks a surface.
    pub fn exposed_areas_m2(&self) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        let mut covered: BTreeMap<_, f64> = self.cells.keys().map(|id| (*id, 0.0)).collect();
        for state in self.contacts.values().filter(|s| s.active) {
            for id in [state.contact.a, state.contact.b] {
                let area = covered.get_mut(&id).ok_or(DepositionError::Invalid("contact endpoint is not born"))?;
                *area = finite(*area + state.contact.area_m2)?;
            }
        }
        for (id, area) in &mut covered {
            let total = self.cells[id].birth.geometry.surface_area_m2;
            let exposed = total - *area;
            // Only a few ulps of area arithmetic may be corrected, never overlap.
            if exposed < -16.0 * f64::EPSILON * total {
                return Err(DepositionError::Invalid("contact areas exceed a parcel surface"));
            }
            *area = exposed.max(0.0);
        }
        Ok(covered)
    }

    /// Activate an event at the CURRENT time. Use advance_to before a future
    /// event. Replaying an accepted event is refused, not silently double-counted.
    pub fn apply_event(
        &mut self, cx: &Cx<'_>, event: &DepositionEvent<G>, energy_tolerance_j: f64,
    ) -> Result<DepositionReport, DepositionError> {
        positive(energy_tolerance_j)?;
        poll(cx, 0)?;
        let mut trial = self.clone();
        let before = trial.stored_enthalpy_j();
        let mut report = DepositionReport::default();
        trial.activate(cx, event, &mut report)?;
        finish(&trial, before, &mut report, energy_tolerance_j)?;
        poll(cx, report.events)?;
        *self = trial;
        Ok(report)
    }

    /// Advance a fixed active domain. Environments remain constant; no future
    /// births are guessed. Calls return only a report and retain only final state.
    pub fn advance_to(
        &mut self, cx: &Cx<'_>, end_time_s: f64, config: DepositionStepConfig,
    ) -> Result<DepositionReport, DepositionError> {
        validate_config(config)?;
        poll(cx, 0)?;
        let mut trial = self.clone();
        let before = trial.stored_enthalpy_j();
        let mut report = DepositionReport::default();
        trial.march(cx, end_time_s, config, &mut report)?;
        finish(&trial, before, &mut report, config.energy_tolerance_j)?;
        poll(cx, report.steps)?;
        *self = trial;
        Ok(report)
    }

    fn activate(
        &mut self, cx: &Cx<'_>, event: &DepositionEvent<G>, report: &mut DepositionReport,
    ) -> Result<(), DepositionError> {
        poll(cx, report.events)?;
        if !event.time_s.is_finite() || event.time_s != self.time_s
            || self.last_event.is_some_and(|last| event.id <= last) {
            return Err(DepositionError::Invalid("events require current time and increasing IDs"));
        }
        charge(&mut report.events, self.budget.max_events, "events")?;
        if event.births.len() > self.budget.max_cells.saturating_sub(self.cells.len()) {
            return Err(DepositionError::Budget("cells"));
        }
        if event.contacts.len() > self.budget.max_contacts || event.environments.len() > self.budget.max_cells {
            return Err(DepositionError::Budget("event changes"));
        }
        for birth in &event.births {
            poll(cx, self.cells.len())?;
            if self.cells.contains_key(&birth.id) { return Err(DepositionError::Invalid("cell identity already exists")); }
            positive(birth.mass_kg)?;
            validate_geometry(&birth.geometry)?;
            positive(finite(birth.mass_kg / birth.geometry.volume_m3)?)?;
            validate_environment(birth.environment)?;
            let material = self.materials.get(&birth.material).ok_or(DepositionError::Invalid("unknown material ID"))?;
            let phase = material.curve.state_at_specific_enthalpy(birth.specific_enthalpy_j_kg)
                .map_err(|e| DepositionError::Phase { cell: birth.id, reason: e.to_string() })?;
            report.incoming_mass_kg = finite(report.incoming_mass_kg + birth.mass_kg)?;
            report.incoming_enthalpy_j = finite(report.incoming_enthalpy_j
                + finite(birth.mass_kg * birth.specific_enthalpy_j_kg)?)?;
            self.cells.insert(birth.id, DepositedCell { birth: birth.clone(), birth_time_s: self.time_s,
                reference_temperature_k: phase.temperature_k(), phase, environment: birth.environment });
        }
        let mut seen = std::collections::BTreeSet::new();
        for change in &event.contacts {
            poll(cx, seen.len())?;
            let id = match change { ContactChange::Set(c) => c.id, ContactChange::Open(id) => *id };
            if !seen.insert(id) { return Err(DepositionError::Invalid("duplicate interface change in event")); }
            match change {
                ContactChange::Set(contact) => {
                    positive(contact.area_m2)?;
                    if contact.a == contact.b || !self.cells.contains_key(&contact.a) || !self.cells.contains_key(&contact.b) {
                        return Err(DepositionError::Invalid("contact needs two distinct born cells"));
                    }
                    if let Some(state) = self.contacts.get_mut(&contact.id) {
                        if state.contact.a != contact.a || state.contact.b != contact.b {
                            return Err(DepositionError::Invalid("interface ID cannot change its oriented endpoints"));
                        }
                        state.contact = contact.clone();
                        state.active = true;
                    } else {
                        if self.contacts.len() >= self.budget.max_contacts { return Err(DepositionError::Budget("contacts")); }
                        self.contacts.insert(contact.id, ContactState { contact: contact.clone(), active: true,
                            first_contact_time_s: self.time_s, heat_a_to_b_j: 0.0 });
                    }
                }
                ContactChange::Open(id) => {
                    self.contacts.get_mut(id).ok_or(DepositionError::Invalid("unknown interface ID"))?.active = false;
                }
            }
        }
        let mut seen = std::collections::BTreeSet::new();
        for &(id, environment) in &event.environments {
            poll(cx, seen.len())?;
            if !seen.insert(id) { return Err(DepositionError::Invalid("duplicate environment change")); }
            validate_environment(environment)?;
            self.cells.get_mut(&id).ok_or(DepositionError::Invalid("environment targets unborn cell"))?.environment = environment;
        }
        self.exposed_areas_m2()?;
        finite(self.total_mass_kg())?;
        finite(self.stored_enthalpy_j())?;
        self.last_event = Some(event.id);
        Ok(())
    }

    fn body(&self, id: CellId, area_m2: f64, convection: f64, emissivity: f64)
        -> Result<LumpedEnthalpyBody<'a>, DepositionError> {
        let cell = &self.cells[&id];
        let material = self.materials[&cell.birth.material];
        Ok(LumpedEnthalpyBody::try_new(format!("deposition-cell-{}", id.0), cell.birth.mass_kg,
            area_m2, convection, emissivity,
            cell.birth.geometry.volume_m3 / cell.birth.geometry.surface_area_m2,
            material.conductivity_w_m_k, material.curve)?)
    }

    fn march(&mut self, cx: &Cx<'_>, end: f64, config: DepositionStepConfig,
        report: &mut DepositionReport) -> Result<(), DepositionError> {
        if !end.is_finite() || end < self.time_s { return Err(DepositionError::Invalid("time must be finite and nondecreasing")); }
        if end == self.time_s { return Ok(()); }
        if self.cells.is_empty() { self.time_s = end; return Ok(()); }
        let areas = self.exposed_areas_m2()?;
        let mut environments = Vec::with_capacity(self.cells.len());
        let mut biot = BTreeMap::new();
        for (&id, cell) in &self.cells {
            poll(cx, environments.len())?;
            let exposed = areas[&id];
            let environment = cell.environment;
            let material = self.materials[&cell.birth.material];
            // A completely buried parcel can still have a source. In that case
            // the dummy whole area has exactly zero boundary-transfer coefficients.
            let body = self.body(id, if exposed > 0.0 { exposed } else { cell.birth.geometry.surface_area_m2 },
                if exposed > 0.0 { environment.convection_w_m2_k } else { 0.0 },
                if exposed > 0.0 { material.emissivity } else { 0.0 })?;
            biot.insert(id, body.maximum_biot(environment.radiation_k));
            environments.push((id, body, environment));
        }
        let mut edges = Vec::new();
        for (&id, state) in self.contacts.iter().filter(|(_, s)| s.active) {
            poll(cx, edges.len())?;
            let c = &state.contact;
            let a = self.body(c.a, c.area_m2, 0.0, 0.0)?;
            let b = self.body(c.b, c.area_m2, 0.0, 0.0)?;
            for (cell_id, body) in [(c.a, &a), (c.b, &b)] {
                let material = self.materials[&self.cells[&cell_id].birth.material];
                let contribution = body.characteristic_length_m()
                    / c.resistance.budget().value_k_per_w / c.area_m2 / material.conductivity_w_m_k;
                let bound = biot.get_mut(&cell_id).expect("active endpoint has a bound");
                *bound = finite(*bound + contribution)?;
            }
            edges.push((id, c.a, c.b, a, b, c.resistance.clone()));
        }
        for bound in biot.values() {
            if !bound.is_finite() || *bound > LUMPED_BIOT_CEILING {
                return Err(DepositionError::Invalid("aggregate Biot exceeds 0.1; use spatial enthalpy transport"));
            }
            report.maximum_biot = report.maximum_biot.max(*bound);
        }
        while self.time_s < end {
            poll(cx, report.steps)?;
            charge(&mut report.steps, self.budget.max_steps, "steps")?;
            let remaining = end - self.time_s;
            let dt = config.maximum_step_s.min(remaining);
            if self.time_s + dt <= self.time_s { return Err(DepositionError::Invalid("timestep cannot advance representable time")); }
            for (id, body, environment) in &environments {
                charge(&mut report.operator_solves, self.budget.max_operator_solves, "operator solves")?;
                let old_h = self.cells[id].phase.specific_enthalpy_j_kg();
                let march = solve_lumped_enthalpy(cx, body, BiotGate::corpus_default(), LumpedEnthalpyMarchConfig {
                    initial_specific_enthalpy_j_kg: old_h,
                    ambient_temperature_k: environment.ambient_k,
                    radiation_temperature_k: environment.radiation_k,
                    internal_power_w: environment.source_w,
                    duration_s: dt, maximum_step_s: dt, maximum_steps: 1,
                    enthalpy_tolerance_j_kg: config.solver_energy_tolerance_j / body.mass_kg(),
                })?;
                let sample = march.samples().last().ok_or(DepositionError::Invalid("empty thermal endpoint"))?;
                if !sample.step_energy_residual_j.is_finite()
                    || sample.step_energy_residual_j.abs() > config.solver_energy_tolerance_j {
                    return Err(DepositionError::Energy { residual_j: sample.step_energy_residual_j,
                        tolerance_j: config.solver_energy_tolerance_j });
                }
                report.boundary_heat_j = finite(report.boundary_heat_j + dt
                    * (sample.convection_into_body_w + sample.radiation_into_body_w))?;
                report.source_heat_j = finite(report.source_heat_j + dt * environment.source_w)?;
                self.cells.get_mut(id).expect("prepared cell exists").phase = sample.phase_state;
            }
            for (id, a_id, b_id, a, b, resistance) in &edges {
                charge(&mut report.operator_solves, self.budget.max_operator_solves, "operator solves")?;
                let step = solve_lumped_contact_step(cx, a, b, resistance, LumpedContactStepConfig {
                    specific_enthalpy_a_j_kg: self.cells[a_id].phase.specific_enthalpy_j_kg(),
                    specific_enthalpy_b_j_kg: self.cells[b_id].phase.specific_enthalpy_j_kg(),
                    supplied_heat_a_j: 0.0, supplied_heat_b_j: 0.0,
                    duration_s: dt, energy_tolerance_j: config.solver_energy_tolerance_j,
                })?;
                self.cells.get_mut(a_id).expect("prepared cell exists").phase = step.state_a;
                self.cells.get_mut(b_id).expect("prepared cell exists").phase = step.state_b;
                let state = self.contacts.get_mut(id).expect("prepared contact exists");
                state.heat_a_to_b_j = finite(state.heat_a_to_b_j + step.heat_a_to_b_j)?;
                report.contact_residual_l1_j = finite(report.contact_residual_l1_j + step.contact_residual_j.abs())?;
            }
            self.time_s = if dt == remaining { end } else { self.time_s + dt };
        }
        Ok(())
    }
}

fn finish<G: Clone>(simulation: &DepositionSimulation<'_, G>, before: f64,
    report: &mut DepositionReport, tolerance: f64) -> Result<(), DepositionError> {
    report.stored_energy_change_j = finite(simulation.stored_enthalpy_j() - before)?;
    report.energy_residual_j = finite(report.stored_energy_change_j
        - report.incoming_enthalpy_j - report.boundary_heat_j - report.source_heat_j)?;
    if report.energy_residual_j.abs() > tolerance {
        return Err(DepositionError::Energy { residual_j: report.energy_residual_j, tolerance_j: tolerance });
    }
    Ok(())
}
fn validate_config(config: DepositionStepConfig) -> Result<(), DepositionError> {
    positive(config.maximum_step_s)?;
    positive(config.solver_energy_tolerance_j)?;
    positive(config.energy_tolerance_j)
}
fn validate_environment(e: Environment) -> Result<(), DepositionError> {
    positive(e.ambient_k)?; positive(e.radiation_k)?; finite(e.source_w)?;
    if !e.convection_w_m2_k.is_finite() || e.convection_w_m2_k < 0.0 {
        return Err(DepositionError::Invalid("convection must be finite and nonnegative"));
    }
    Ok(())
}
fn validate_geometry<G>(g: &CellGeometry<G>) -> Result<(), DepositionError> {
    positive(g.volume_m3)?; positive(g.surface_area_m2)?;
    positive(finite(g.volume_m3 / g.surface_area_m2)?)?;
    for value in g.position_m { finite(value)?; }
    for i in 0..3 {
        for j in 0..3 {
            let dot = sum((0..3).map(|k| g.material_axes[i][k] * g.material_axes[j][k]));
            if !dot.is_finite() || (dot - if i == j { 1.0 } else { 0.0 }).abs() > 1e-10 {
                return Err(DepositionError::Invalid("material axes must be orthonormal"));
            }
        }
    }
    let [a, b, c] = g.material_axes;
    let det = a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
    if (det - 1.0).abs() > 1e-10 { return Err(DepositionError::Invalid("material axes must be right-handed")); }
    Ok(())
}
fn poll(cx: &Cx<'_>, at: usize) -> Result<(), DepositionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage: "deposition", at })?;
    Ok(())
}
fn charge(counter: &mut usize, limit: usize, name: &'static str) -> Result<(), DepositionError> {
    if *counter >= limit { return Err(DepositionError::Budget(name)); }
    *counter += 1; Ok(())
}
fn positive(value: f64) -> Result<(), DepositionError> {
    if value.is_finite() && value > 0.0 { Ok(()) }
    else { Err(DepositionError::Invalid("positive finite SI quantity required")) }
}
fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() { Ok(value) } else { Err(DepositionError::Invalid("nonfinite arithmetic")) }
}
fn sum(values: impl Iterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values { let y = value - correction; let next = total + y;
        correction = (next - total) - y; total = next; }
    total
}
