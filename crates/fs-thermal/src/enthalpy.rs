//! Growing-domain thermal transport for resolved material-extrusion schedules.
//!
//! Cells are caller-resolved control volumes, not a slicer or bead-forming law.
//! Each has one uniform enthalpy and temperature. Internal links exchange the
//! SAME signed joules on their two sides. Birth adds exactly `mass * h_in`;
//! previously deposited cells retain their history. Contact and bed areas are
//! removed from the free surface, including for zero-conductance contacts.
//! Geometry is fixed between explicit events; no deformation or remeshing is
//! inferred. Cell, segment, body and material identities survive thermal steps.
//!
//! Forward Euler is subcycled with a declared maximum step and a monotonicity
//! bound `dt <= 0.8 m / (max(dT/dh) * sum(outgoing conductance))`. All event
//! times are hit exactly. No future cell conducts heat. A failed/cancelled
//! advance publishes no candidate. Cost is sparse, O(cells + links + events) storage,
//! O((cells + links) log(cells)) per substep; no dense global solve or full
//! history retention. Time and spatial resolution are independent caller
//! choices. This reduced model needs a uniform-cell/Biot justification and
//! independently supplied conductances; it is not a validated print predictor.

use std::collections::BTreeMap;

// Keep the std-only kernel testable as a standalone crate as well as a module.
#[path = "deposition/observe.rs"]
mod observe;
pub use observe::{ContactUpdate, DepositionObserver, PowerBalance};
#[path = "deposition/coverage.rs"]
mod coverage;
pub use coverage::{ContactPatch, ContactPatchUpdate};

/// Endpoint residuals sharing the existing storage, flux and boundary owner.
#[path = "enthalpy/backward_euler.rs"]
pub mod backward_euler;

/// Conservative enthalpy transport driven by a balanced, prescribed mass flow.
#[path = "enthalpy/advection.rs"]
pub mod advection;

/// Opt-in entropy accounting on the original enthalpy and heat-transfer laws.
#[path = "enthalpy/entropy.rs"]
pub mod entropy;

/// Conservative sublimation and deposition on this same thermal timeline.
#[path = "enthalpy/sublimation.rs"]
pub mod sublimation;

/// Persistent physical cell identity; distinct from its source segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CellId(pub u64);

/// Persistent material identity, never an index into a reordered table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct MaterialId(pub u64);

/// Replaceable enthalpy/temperature law. All materials must use a consistent
/// specific-energy datum for an energy audit to have physical meaning. `T(h)`
/// must be continuous and nondecreasing over the admitted domain, and the law
/// must remain immutable throughout an advance and checkpoint/restart.
pub trait ThermalLaw {
    /// Temperature in kelvin at specific enthalpy in J/kg. Refuse extrapolation.
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError>;
    /// Global upper bound on dT/dh [kg K/J] over the admitted enthalpy domain.
    /// Must be finite and nonnegative; an underestimated bound invalidates CFL.
    /// Admitted once at timeline construction and reused across advances and
    /// checkpoints; the bound must remain valid for the immutable law.
    fn maximum_temperature_slope(&self) -> f64;
    /// Exact relative thermal entropy integral `integral dh / T(h)` [J/(kg K)]
    /// on this same admitted caloric law, at its fixed pressure/composition.
    /// `None` explicitly means unavailable; the transport owner never invents
    /// a heat capacity, entropy datum, mixing law or extrapolated chart.
    /// Implementations must retain latent heat in this integral exactly once.
    fn specific_entropy_change_j_kg_k(
        &self, _initial_h_j_kg: f64, _final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> { Ok(None) }
    /// Optional frozen mixture/reference binding for stationary species advection.
    /// The original single-material path does not require or infer this data.
    fn advected_composition(&self) -> Option<advection::species::CaloricComposition<'_>> { None }
}

/// Declared reservoir contact, independent of mechanical adhesion or strength.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Reservoir {
    /// Covered surface area, m².
    pub area_m2: f64,
    /// Total contact conductance, W/K (zero means thermally insulating).
    pub conductance_w_k: f64,
    /// Reservoir temperature, K.
    pub temperature_k: f64,
}

/// Endpoint environment held constant until an explicit boundary event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Boundary {
    /// Free-surface ambient temperature, K.
    pub ambient_k: f64,
    /// Free-surface convection coefficient, W/(m² K), including zero.
    pub htc_w_m2_k: f64,
    /// Optional bed/support reservoir. Its area is excluded from ambient loss.
    pub reservoir: Option<Reservoir>,
    /// Applied total cell heat source, W, positive inward.
    pub source_w: f64,
}

/// Located deposition geometry in one caller-declared world frame, SI metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellGeometry {
    /// Centreline start, m.
    pub start_m: [f64; 3],
    /// Centreline end, m. Direction supplies a later constitutive material frame.
    pub end_m: [f64; 3],
    /// Deposited control-volume volume, m³. Mass is independently supplied.
    pub volume_m3: f64,
    /// Total surface area before subtracting contacts, m².
    pub surface_area_m2: f64,
}

/// Incoming cell, activated once at its event time. A segment can contain many
/// cells. Volume and delivered mass imply a reference density; equilibrium
/// density never silently changes stored mass.
#[derive(Debug, Clone, PartialEq)]
pub struct DepositedCell {
    /// Stable cell identifier.
    pub id: CellId,
    /// Application's source segment identifier.
    pub segment_id: u64,
    /// Application's body/part/support identifier.
    pub body_id: u64,
    /// Explicit material-law identity.
    pub material_id: MaterialId,
    /// Delivered mass, kg, strictly positive.
    pub mass_kg: f64,
    /// Incoming specific enthalpy, J/kg. Uses the selected law's energy datum.
    pub specific_h_j_kg: f64,
    /// Located geometry.
    pub geometry: CellGeometry,
    /// Birth-time environment.
    pub boundary: Boundary,
}

/// Thermal contact between two active cells. Conductance, geometric contact
/// and future healing/damage are deliberately separate physical quantities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalLink {
    /// Stable interface identifier for later bond/contact histories.
    pub id: u64,
    /// First incident cell.
    pub a: CellId,
    /// Second incident cell.
    pub b: CellId,
    /// Default covered area on EACH side, m². A `ContactPatch` can instead
    /// supply separate covered areas; then this is its nominal interface area.
    pub area_m2: f64,
    /// Total thermal conductance, W/K.
    pub conductance_w_k: f64,
}

/// Bulk conductance [W/K] of two homogeneous isotropic layers in series.
///
/// Distances are positive layer thicknesses [m], conductivity is [W/(m K)],
/// and area is the shared transfer area [m2]. Geometry and material/source
/// admission belong to the caller; this adds no contact resistance. Explicit
/// zero area or conductivity insulates, but does not excuse invalid inputs.
/// Unrepresentable positive resistances or conductance refuse.
pub fn series_conductance_w_k(
    area_m2: f64,
    distances_m: [f64; 2],
    conductivity_w_m_k: [f64; 2],
) -> Result<f64, DepositionError> {
    series_contact_conductance_w_k(area_m2, distances_m, conductivity_w_m_k, 0.0)
}

/// Conductance of the same two bulk layers plus one interface in series.
///
/// `area_resistance_m2_k_w` is the source-resolved area-specific contact
/// resistance [m² K/W]. Zero explicitly omits that interface resistance and
/// retains the original two-layer arithmetic. The caller owns the ordered
/// materials, source support and physical contact geometry. Invalid resistance
/// refuses even for an explicitly insulated face.
pub fn series_contact_conductance_w_k(
    area_m2: f64,
    distances_m: [f64; 2],
    conductivity_w_m_k: [f64; 2],
    area_resistance_m2_k_w: f64,
) -> Result<f64, DepositionError> {
    if !nonnegative(area_m2)
        || !distances_m.into_iter().all(positive)
        || !conductivity_w_m_k.into_iter().all(nonnegative)
        || !nonnegative(area_resistance_m2_k_w)
    {
        return Err(DepositionError::InvalidInput(
            "series conduction needs finite area/conductivity and positive distances",
        ));
    }
    if area_m2 == 0.0 || conductivity_w_m_k.contains(&0.0) {
        return Ok(0.0);
    }
    let ra = finite(distances_m[0] / conductivity_w_m_k[0])?;
    let rb = finite(distances_m[1] / conductivity_w_m_k[1])?;
    let bulk = finite(ra + rb)?;
    let total = if area_resistance_m2_k_w == 0.0 {
        bulk
    } else {
        finite(bulk + area_resistance_m2_k_w)?
    };
    let conductance = finite(area_m2 / total)?;
    if ra <= 0.0 || rb <= 0.0 || total <= 0.0 || conductance <= 0.0 {
        return Err(DepositionError::InvalidInput(
            "unrepresentable material conduction",
        ));
    }
    Ok(conductance)
}

/// A signed INTERNAL energy transfer, independent of contact area/coverage.
/// The identical watts are subtracted from `a` and added to `b`; they never
/// enter the external-heat account. Both endpoints must be distinct active cells.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatTransfer {
    /// First receiving control volume.
    pub a: CellId,
    /// Second receiving control volume.
    pub b: CellId,
    /// Signed power [W], positive from a to b.
    pub power_w: f64,
}

/// Resolved physical event. Equal-time events execute in supplied order; put
/// births before their contacts. No printer command interpretation is done.
#[derive(Debug, Clone, PartialEq)]
pub enum Event {
    /// Add a new control volume, including incoming mass and energy.
    Birth(DepositedCell),
    /// Add a thermal interface and cover its ambient-facing surface areas.
    Contact(ThermalLink),
    /// Remove a thermal interface and expose its two surfaces again.
    Separation(u64),
    /// Replace endpoint environment, including reservoir removal or fan change.
    Boundary(CellId, Boundary),
}

/// One absolute physical time in seconds, nonnegative and finite.
#[derive(Debug, Clone, PartialEq)]
pub struct TimedEvent {
    /// Activation/change time, s.
    pub time_s: f64,
    /// Resolved change.
    pub event: Event,
}

/// Explicit work/memory ceilings for the schedule and each advance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DepositionBudget {
    /// Maximum active cells.
    pub max_cells: usize,
    /// Maximum active interfaces.
    pub max_links: usize,
    /// Maximum retained resolved schedule events.
    pub max_events: usize,
    /// Maximum thermal substeps in one advance (events do not consume steps).
    pub max_substeps: usize,
}

/// Timestep and absolute conservation admission, separate from model accuracy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepConfig {
    /// Maximum thermal step, s. CFL may reduce it further.
    pub max_dt_s: f64,
    /// Maximum stored-minus-incoming-minus-external defect, J.
    /// Also bounds the cumulative absolute local work too small to change a
    /// cell's enthalpy coordinate during an explicit advance. Opposing losses
    /// count separately, so a zero global residual cannot conceal them.
    pub energy_tolerance_j: f64,
}

/// A refusal leaves the original timeline and all material history unchanged.
#[derive(Debug, Clone, PartialEq)]
pub enum DepositionError {
    /// Bad units, topology, duplicate IDs, range, or timeline.
    InvalidInput(&'static str),
    /// Declared allocation/work ceiling exhausted.
    Budget(&'static str),
    /// Caller cancellation gate refused work.
    Cancelled,
    /// Constitutive chart refused an enthalpy query.
    MaterialDomain,
    /// Finite inputs produced nonfinite arithmetic or stalled time.
    Arithmetic,
    /// Absolute conservation gate failed, with signed residual in joules.
    EnergyBalance(f64),
}

impl std::fmt::Display for DepositionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "deposition thermal refusal: {self:?}")
    }
}
impl std::error::Error for DepositionError {}

/// Per-advance conservation and cost observations (not physical error bounds).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyAudit {
    /// Mass added by activated births, kg.
    pub incoming_mass_kg: f64,
    /// Energy carried by activated births, J.
    pub incoming_energy_j: f64,
    /// Integrated ambient/reservoir/source input, J, positive inward.
    pub external_heat_j: f64,
    /// Actual change in stored `sum(m*h)`, J.
    pub stored_energy_change_j: f64,
    /// Stored change minus incoming enthalpy minus external input, J.
    /// Evaluated from compensated per-cell initial/birth-relative changes.
    /// Recomputing it from the rounded aggregate report fields can lose small
    /// heat increments against large formation-enthalpy reference energies.
    pub residual_j: f64,
    /// Number of accepted thermal steps (explicit substeps or implicit intervals).
    pub substeps: usize,
    /// Number of accepted physical events.
    pub events: usize,
}

/// Immutable-candidate timeline. Clone is an in-memory checkpoint of the
/// schedule cursor, active geometry, interfaces and enthalpy, not a history log.
#[derive(Clone)]
pub struct DepositionSimulation<'a> {
    time_s: f64,
    cursor: usize,
    schedule: Vec<TimedEvent>,
    materials: BTreeMap<MaterialId, &'a dyn ThermalLaw>,
    // Immutable global CFL bounds, admitted once rather than rescanning a
    // material chart for every active cell at every thermal substep.
    material_temperature_slopes: BTreeMap<MaterialId, f64>,
    cells: BTreeMap<CellId, DepositedCell>,
    links: BTreeMap<u64, ThermalLink>,
    // Overrides only; symmetric links retain their existing area semantics.
    covered_areas: BTreeMap<u64, [f64; 2]>,
    budget: DepositionBudget,
}

fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

impl<'a> DepositionSimulation<'a> {
    /// Start an empty physical timeline at t=0. Initial bed/support/existing
    /// material can be declared as births at t=0 with their known enthalpy.
    /// `advance_to(0)` activates those events without thermal aging.
    pub fn new(
        schedule: Vec<TimedEvent>,
        materials: &[(MaterialId, &'a dyn ThermalLaw)],
        budget: DepositionBudget,
    ) -> Result<Self, DepositionError> {
        if schedule.len() > budget.max_events || materials.len() > budget.max_cells {
            return Err(DepositionError::Budget("schedule/materials"));
        }
        let mut previous = 0.0;
        for event in &schedule {
            if !nonnegative(event.time_s) || event.time_s < previous {
                return Err(DepositionError::InvalidInput(
                    "events must be ordered in absolute seconds",
                ));
            }
            previous = event.time_s;
        }
        let mut laws = BTreeMap::new();
        let mut material_temperature_slopes = BTreeMap::new();
        for &(id, law) in materials {
            let slope = law.maximum_temperature_slope();
            if !nonnegative(slope) || laws.insert(id, law).is_some() {
                return Err(DepositionError::InvalidInput(
                    "duplicate material or invalid temperature slope",
                ));
            }
            material_temperature_slopes.insert(id, slope);
        }
        Ok(Self {
            time_s: 0.0,
            cursor: 0,
            schedule,
            materials: laws,
            material_temperature_slopes,
            cells: BTreeMap::new(),
            links: BTreeMap::new(),
            covered_areas: BTreeMap::new(),
            budget,
        })
    }

    /// Current accepted physical time, s.
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Current material state, located geometry and source correspondence.
    #[must_use]
    pub fn cells(&self) -> &BTreeMap<CellId, DepositedCell> {
        &self.cells
    }
    /// Current exchange links. Use `contact_patch` for authoritative per-side
    /// coverage: a link alone does not encode an asymmetric material footprint.
    #[must_use]
    pub fn links(&self) -> &BTreeMap<u64, ThermalLink> {
        &self.links
    }
    /// Stored reference mass, kg.
    #[must_use]
    pub fn mass_kg(&self) -> f64 {
        self.cells.values().map(|c| c.mass_kg).sum()
    }
    /// Stored thermal energy relative to declared material data, J.
    #[must_use]
    pub fn energy_j(&self) -> f64 {
        self.cells
            .values()
            .map(|c| c.mass_kg * c.specific_h_j_kg)
            .sum()
    }
    /// Query temperature without retaining every timestep's full field.
    pub fn temperature_k(&self, id: CellId) -> Result<f64, DepositionError> {
        let cell = self
            .cells
            .get(&id)
            .ok_or(DepositionError::InvalidInput("inactive cell"))?;
        let t = self.materials[&cell.material_id].temperature_k(cell.specific_h_j_kg)?;
        if !positive(t) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(t)
    }
    /// Exposed ambient surface, m², after subtracting current contacts.
    pub fn exposed_area_m2(&self, id: CellId) -> Result<f64, DepositionError> {
        let cell = self
            .cells
            .get(&id)
            .ok_or(DepositionError::InvalidInput("inactive cell"))?;
        let mut area =
            cell.geometry.surface_area_m2 - cell.boundary.reservoir.map_or(0.0, |r| r.area_m2);
        for link in self.links.values() {
            if link.a == id || link.b == id {
                let covered = self.contact_areas(link);
                area -= if link.a == id { covered[0] } else { covered[1] };
            }
        }
        if !nonnegative(area) {
            return Err(DepositionError::InvalidInput(
                "contacts exceed surface area",
            ));
        }
        Ok(area)
    }

    /// Advance to an absolute time, hitting every birth/boundary event exactly.
    /// `checkpoint` is polled every 256 cells/links in thermal loops and between events/steps;
    /// use the Cx adapter in `enthalpy::deposition` for executor integration.
    /// Only the returned candidate commits aging, births and schedule progress.
    pub fn advance_to(
        &self,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, EnergyAudit), DepositionError> {
        self.advance_observed(&(), target_s, config, checkpoint)
            .map(|(state, (), audit)| (state, audit))
    }

    fn apply(&mut self, event: Event, audit: &mut EnergyAudit) -> Result<(), DepositionError> {
        match event {
            Event::Birth(cell) => {
                if self.cells.len() >= self.budget.max_cells {
                    return Err(DepositionError::Budget("active cells"));
                }
                if self.cells.contains_key(&cell.id)
                    || !self.materials.contains_key(&cell.material_id)
                    || !positive(cell.mass_kg)
                    || !positive(cell.geometry.volume_m3)
                    || !positive(cell.geometry.surface_area_m2)
                    || !cell
                        .geometry
                        .start_m
                        .iter()
                        .chain(&cell.geometry.end_m)
                        .all(|x| x.is_finite())
                    || cell.geometry.start_m == cell.geometry.end_m
                    || !cell.specific_h_j_kg.is_finite()
                {
                    return Err(DepositionError::InvalidInput(
                        "invalid birth geometry/mass/material or reused cell ID",
                    ));
                }
                validate_boundary(cell.boundary)?;
                audit.incoming_mass_kg = finite(audit.incoming_mass_kg + cell.mass_kg)?;
                audit.incoming_energy_j =
                    finite(audit.incoming_energy_j + finite(cell.mass_kg * cell.specific_h_j_kg)?)?;
                let id = cell.id;
                self.cells.insert(id, cell);
                self.temperature_k(id)?;
                self.exposed_area_m2(id)?;
            }
            Event::Contact(link) => {
                self.insert_contact_patch(ContactPatch::symmetric(link))?;
            }
            Event::Separation(id) => {
                if self.links.remove(&id).is_none() {
                    return Err(DepositionError::InvalidInput("inactive interface"));
                }
                self.covered_areas.remove(&id);
            }
            Event::Boundary(id, boundary) => {
                validate_boundary(boundary)?;
                self.cells
                    .get_mut(&id)
                    .ok_or(DepositionError::InvalidInput("inactive boundary cell"))?
                    .boundary = boundary;
                self.exposed_area_m2(id)?;
            }
        }
        Ok(())
    }

    // A single owner for physical watts and the explicit stability bound.
    // Endpoint residuals consume these watts directly, not differences between
    // two advanced states or a second implementation of boundary/contact heat.
    fn heat_rates(
        &self,
        maximum_explicit_dt: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<HeatRates, DepositionError> {
        // Per-cell rates and outgoing conductances; no quadratic cell searches.
        let mut temperatures = BTreeMap::new();
        let mut rates = BTreeMap::new();
        let mut outgoing = BTreeMap::new();
        let mut covered = BTreeMap::new();
        for (i, (&id, cell)) in self.cells.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            temperatures.insert(id, self.temperature_k(id)?);
            covered.insert(id, cell.boundary.reservoir.map_or(0.0, |r| r.area_m2));
            rates.insert(id, 0.0);
            outgoing.insert(id, 0.0);
        }
        for (i, link) in self.links.values().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let power =
                finite(link.conductance_w_k * (temperatures[&link.b] - temperatures[&link.a]))?;
            let coverage = self.contact_areas(link);
            for (id, sign, side_area) in [(link.a, 1.0, coverage[0]), (link.b, -1.0, coverage[1])] {
                let rate = rates.get_mut(&id).ok_or(DepositionError::Arithmetic)?;
                *rate = finite(*rate + sign * power)?;
                let g = outgoing.get_mut(&id).ok_or(DepositionError::Arithmetic)?;
                *g = finite(*g + link.conductance_w_k)?;
                let area = covered.get_mut(&id).ok_or(DepositionError::Arithmetic)?;
                *area = finite(*area + side_area)?;
            }
        }
        let mut external_w = 0.0;
        let mut dt = maximum_explicit_dt;
        for (i, (&id, cell)) in self.cells.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let area = cell.geometry.surface_area_m2 - covered[&id];
            if !nonnegative(area) {
                return Err(DepositionError::InvalidInput(
                    "contacts exceed surface area",
                ));
            }
            let ambient_g = finite(area * cell.boundary.htc_w_m2_k)?;
            let mut g = finite(outgoing[&id] + ambient_g)?;
            let mut external = finite(
                cell.boundary.source_w + ambient_g * (cell.boundary.ambient_k - temperatures[&id]),
            )?;
            if let Some(r) = cell.boundary.reservoir {
                g = finite(g + r.conductance_w_k)?;
                external =
                    finite(external + r.conductance_w_k * (r.temperature_k - temperatures[&id]))?;
            }
            external_w = finite(external_w + external)?;
            let rate = rates.get_mut(&id).ok_or(DepositionError::Arithmetic)?;
            *rate = finite(*rate + external)?;
            let slope = self.material_temperature_slopes[&cell.material_id];
            let decay = finite(g * slope / cell.mass_kg)?;
            if decay > 0.0 {
                dt = dt.min(0.8 / decay);
            }
        }
        if !positive(dt) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(HeatRates { temperatures, watts: rates, external_w, explicit_dt_s: dt })
    }

    fn thermal_step(
        &mut self,
        rates: HeatRates,
        audit: &mut EnergyAudit,
        unrepresented_heat_j: &mut f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(f64, BTreeMap<CellId, f64>), DepositionError> {
        let HeatRates { temperatures, watts: rates, external_w, explicit_dt_s: dt } =
            rates;
        for (i, (&id, cell)) in self.cells.iter_mut().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let old_h = cell.specific_h_j_kg;
            let next_h = finite(old_h + finite(dt * rates[&id] / cell.mass_kg)?)?;
            // Bound local work that disappears on the enthalpy coordinate.
            // Both sides of an internal exchange count; cancellation cannot
            // conceal loss. The caller applies its absolute joule tolerance.
            if rates[&id] != 0.0 && next_h == old_h {
                let lost = finite(dt * rates[&id])?.abs();
                *unrepresented_heat_j = finite((*unrepresented_heat_j + lost).next_up())?;
            }
            cell.specific_h_j_kg = next_h;
            let t = self.materials[&cell.material_id].temperature_k(cell.specific_h_j_kg)?;
            if !positive(t) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        audit.external_heat_j = finite(audit.external_heat_j + finite(dt * external_w)?)?;
        Ok((dt, temperatures))
    }
}

struct HeatRates {
    temperatures: BTreeMap<CellId, f64>,
    watts: BTreeMap<CellId, f64>,
    external_w: f64,
    explicit_dt_s: f64,
}

fn validate_boundary(b: Boundary) -> Result<(), DepositionError> {
    if !positive(b.ambient_k)
        || !nonnegative(b.htc_w_m2_k)
        || !b.source_w.is_finite()
        || b.reservoir.is_some_and(|r| {
            !nonnegative(r.area_m2) || !nonnegative(r.conductance_w_k) || !positive(r.temperature_k)
        })
    {
        return Err(DepositionError::InvalidInput(
            "invalid thermal boundary in SI units",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn series_conduction_retains_both_layers_and_explicit_insulation() {
        let g = series_conductance_w_k(6.0, [0.25, 0.5], [12.0, 10.0]).unwrap();
        close(g, 6.0 / (0.25 / 12.0 + 0.5 / 10.0), 1e-12);
        assert_eq!(series_conductance_w_k(0.0, [1.0; 2], [1.0; 2]), Ok(0.0));
        assert_eq!(series_conductance_w_k(1.0, [1.0; 2], [0.0, 1.0]), Ok(0.0));
        assert!(series_conductance_w_k(0.0, [f64::NAN, 1.0], [0.0, 1.0]).is_err());
        assert!(series_conductance_w_k(1.0, [1.0; 2], [-1.0, 1.0]).is_err());
        assert!(series_conductance_w_k(1.0, [f64::MIN_POSITIVE; 2], [f64::MAX; 2],).is_err());
    }

    #[test]
    fn contact_resistance_composes_once_with_bulk_and_checks_insulated_inputs() {
        let bulk = series_conductance_w_k(6.0, [0.25, 0.5], [12.0, 10.0]).unwrap();
        let zero = series_contact_conductance_w_k(6.0, [0.25, 0.5], [12.0, 10.0], 0.0).unwrap();
        assert_eq!(bulk.to_bits(), zero.to_bits());
        let contact = series_contact_conductance_w_k(6.0, [0.25, 0.5], [12.0, 10.0], 0.1).unwrap();
        close(
            1.0 / contact,
            0.25 / (6.0 * 12.0) + 0.5 / (6.0 * 10.0) + 0.1 / 6.0,
            1e-14,
        );
        assert!(contact < bulk);
        assert_eq!(
            series_contact_conductance_w_k(0.0, [1.0; 2], [1.0; 2], 0.1),
            Ok(0.0)
        );
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(series_contact_conductance_w_k(0.0, [1.0; 2], [1.0; 2], invalid).is_err());
        }
    }

    struct Sensible;
    impl ThermalLaw for Sensible {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if !(0.0..=1_000_000.0).contains(&h) {
                return Err(DepositionError::MaterialDomain);
            }
            Ok(300.0 + h / 2000.0)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0 / 2000.0
        }
    }
    static LAW: Sensible = Sensible;
    fn budget() -> DepositionBudget {
        DepositionBudget {
            max_cells: 10_000,
            max_links: 10_000,
            max_events: 30_000,
            max_substeps: 100_000,
        }
    }
    fn config(dt: f64) -> StepConfig {
        StepConfig {
            max_dt_s: dt,
            energy_tolerance_j: 1e-8,
        }
    }
    fn boundary(htc: f64) -> Boundary {
        Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: htc,
            reservoir: None,
            source_w: 0.0,
        }
    }
    fn cell(id: u64, h: f64) -> DepositedCell {
        DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: 1,
            material_id: MaterialId(7),
            mass_kg: 0.001,
            specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [0.0, 0.0, id as f64 * 0.0002],
                end_m: [0.01, 0.0, id as f64 * 0.0002],
                volume_m3: 1e-6,
                surface_area_m2: 0.001,
            },
            boundary: boundary(0.0),
        }
    }
    fn event(t: f64, event: Event) -> TimedEvent {
        TimedEvent { time_s: t, event }
    }
    fn birth(t: f64, id: u64, h: f64) -> TimedEvent {
        event(t, Event::Birth(cell(id, h)))
    }
    fn link(a: u64, b: u64) -> ThermalLink {
        ThermalLink {
            id: a,
            a: CellId(a),
            b: CellId(b),
            area_m2: 0.0002,
            conductance_w_k: 0.2,
        }
    }
    fn simulation(events: Vec<TimedEvent>) -> DepositionSimulation<'static> {
        DepositionSimulation::new(events, &[(MaterialId(7), &LAW)], budget()).unwrap()
    }
    fn run(
        sim: &DepositionSimulation<'static>,
        t: f64,
        dt: f64,
    ) -> (DepositionSimulation<'static>, EnergyAudit) {
        sim.advance_to(t, config(dt), || Ok(())).unwrap()
    }
    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a} != {b}");
    }

    #[test]
    fn immutable_material_slopes_are_reused_across_cells_steps_and_checkpoints() {
        struct CountedLaw(std::cell::Cell<usize>);
        impl ThermalLaw for CountedLaw {
            fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
                LAW.temperature_k(h)
            }
            fn maximum_temperature_slope(&self) -> f64 {
                self.0.set(self.0.get() + 1);
                LAW.maximum_temperature_slope()
            }
        }
        let laws = [
            CountedLaw(std::cell::Cell::new(0)),
            CountedLaw(std::cell::Cell::new(0)),
        ];
        let mut other = cell(2, 400_000.0);
        other.material_id = MaterialId(8);
        let schedule = vec![
            birth(0.0, 1, 0.0),
            event(0.0, Event::Birth(other)),
            event(0.0, Event::Contact(link(1, 2))),
            event(0.5, Event::Boundary(CellId(2), boundary(100.0))),
        ];
        let mut counted = DepositionSimulation::new(
            schedule.clone(),
            &[(MaterialId(7), &laws[0]), (MaterialId(8), &laws[1])],
            budget(),
        )
        .unwrap();
        let mut reference = DepositionSimulation::new(
            schedule,
            &[(MaterialId(7), &LAW), (MaterialId(8), &LAW)],
            budget(),
        )
        .unwrap();
        assert!(laws.iter().all(|law| law.0.get() == 1));
        for target in [0.5, 1.0] {
            let (next, audit) = counted.advance_to(target, config(0.01), || Ok(())).unwrap();
            let (expected, expected_audit) =
                reference.advance_to(target, config(0.01), || Ok(())).unwrap();
            assert!(audit.substeps > 1);
            assert_eq!(next.cells(), expected.cells());
            assert_eq!(next.links(), expected.links());
            assert_eq!(audit, expected_audit);
            counted = next.clone();
            reference = expected;
        }
        assert!(laws.iter().all(|law| law.0.get() == 1));
    }

    // G0: birth energy is exactly once, incoming reference mass never inferred
    // from a material's changing equilibrium density, and no cell ages at birth.
    #[test]
    fn single_bead_mass_energy_and_no_preactivation() {
        let sim = simulation(vec![birth(0.75, 1, 400_000.0)]);
        let (before, a) = run(&sim, 0.5, 9.0);
        assert!(before.cells().is_empty());
        assert_eq!(a.incoming_energy_j, 0.0);
        let (born, a) = run(&before, 0.75, 9.0);
        close(born.mass_kg(), 0.001, 1e-15);
        close(born.energy_j(), 400.0, 1e-12);
        close(born.temperature_k(CellId(1)).unwrap(), 500.0, 1e-12);
        assert_eq!(a.events, 1);
        let (same, a) = run(&born, 0.75, 9.0);
        assert_eq!(a.events, 0);
        assert_eq!(same.cells(), born.cells());
    }

    // G1: independent lumped cooling solution; timestep refinement is distinct
    // from spatial resolution. tau=m cp/(h A)=20 seconds here.
    #[test]
    fn single_bead_cooling_matches_exponential_with_first_order_refinement() {
        let mut c = cell(1, 400_000.0);
        c.boundary = boundary(100.0);
        let sim = simulation(vec![event(0.0, Event::Birth(c))]);
        let exact = 300.0 + 200.0 * (-1.0_f64).exp();
        let mut errors = Vec::new();
        for dt in [1.0, 0.5, 0.25] {
            let (s, a) = run(&sim, 20.0, dt);
            errors.push((s.temperature_k(CellId(1)).unwrap() - exact).abs());
            close(
                a.stored_energy_change_j,
                a.incoming_energy_j + a.external_heat_j,
                1e-9,
            );
        }
        assert!((1.9..2.1).contains(&(errors[0] / errors[1])));
        assert!((1.9..2.1).contains(&(errors[1] / errors[2])));
    }

    #[test]
    fn successive_layer_reheats_old_material_and_changes_exposure() {
        let sim = simulation(vec![
            birth(0.0, 1, 40_000.0),
            birth(1.0, 2, 400_000.0),
            event(1.0, Event::Contact(link(1, 2))),
            event(3.0, Event::Separation(1)),
        ]);
        let (first, _) = run(&sim, 0.9, 0.01);
        assert_eq!(first.cells().len(), 1);
        let (second, a) = run(&first, 2.0, 0.01);
        assert!(second.temperature_k(CellId(1)).unwrap() > first.temperature_k(CellId(1)).unwrap());
        close(second.exposed_area_m2(CellId(1)).unwrap(), 0.0008, 1e-15);
        close(second.energy_j(), 440.0, 1e-9);
        close(a.incoming_energy_j, 400_000.0 * 0.001, 1e-12);
        close(a.residual_j, 0.0, 1e-9);
        assert_eq!(second.cells()[&CellId(1)].segment_id, 1);
        let (separated, a) = run(&second, 3.0, 0.01);
        close(separated.exposed_area_m2(CellId(1)).unwrap(), 0.001, 1e-15);
        close(a.incoming_energy_j, 0.0, 0.0);
    }

    // G1: independent two-body fixed-domain exchange. Equal capacities C=2 J/K
    // give a temperature difference decaying as exp(-2 G t/C).
    #[test]
    fn fixed_domain_contact_limit_matches_two_body_reference() {
        let sim = simulation(vec![
            birth(0.0, 1, 0.0),
            birth(0.0, 2, 400_000.0),
            event(0.0, Event::Contact(link(1, 2))),
        ]);
        let (s, a) = run(&sim, 10.0, 0.001);
        let exact = 400.0 - 100.0 * (-2.0_f64).exp();
        close(s.temperature_k(CellId(1)).unwrap(), exact, 0.004);
        close(s.energy_j(), 400.0, 1e-8);
        close(a.external_heat_j, 0.0, 0.0);
        close(a.residual_j, 0.0, 1e-8);
    }

    #[test]
    fn cfl_subcycling_prevents_contact_temperature_overshoot() {
        let sim = simulation(vec![
            birth(0.0, 1, 0.0),
            birth(0.0, 2, 400_000.0),
            event(0.0, Event::Contact(link(1, 2))),
        ]);
        let (s, a) = run(&sim, 100.0, 1000.0);
        assert!(a.substeps > 1);
        for id in [CellId(1), CellId(2)] {
            assert!((300.0..=500.0).contains(&s.temperature_k(id).unwrap()));
        }
        close(s.energy_j(), 400.0, 1e-9);
    }

    #[test]
    fn bed_release_and_fan_change_hit_exact_timeline() {
        let mut c = cell(1, 0.0);
        c.boundary.reservoir = Some(Reservoir {
            area_m2: 0.0004,
            conductance_w_k: 0.1,
            temperature_k: 400.0,
        });
        let sim = simulation(vec![
            event(0.0, Event::Birth(c)),
            event(0.35, Event::Boundary(CellId(1), boundary(100.0))),
        ]);
        let (released, _) = run(&sim, 0.35, 0.2);
        assert!(released.temperature_k(CellId(1)).unwrap() > 300.0);
        assert_eq!(released.time_s(), 0.35);
        close(released.exposed_area_m2(CellId(1)).unwrap(), 0.001, 1e-15);
        let (cooled, a) = run(&released, 1.0, 0.2);
        assert!(
            cooled.temperature_k(CellId(1)).unwrap() < released.temperature_k(CellId(1)).unwrap()
        );
        assert!(a.external_heat_j < 0.0);
    }

    // G3: partitioning a uniform bead changes spatial cell count but not its
    // mass, energy or temperature. This proves this uniform limit only.
    #[test]
    fn uniform_bead_spatial_partition_preserves_solution() {
        let mut final_t = Vec::new();
        for count in [1, 2, 8] {
            let events = (0..count)
                .map(|id| {
                    let mut c = cell(id, 400_000.0);
                    c.mass_kg /= count as f64;
                    c.geometry.volume_m3 /= count as f64;
                    c.geometry.surface_area_m2 /= count as f64;
                    c.boundary = boundary(100.0);
                    event(0.0, Event::Birth(c))
                })
                .collect();
            let (s, a) = run(&simulation(events), 10.0, 0.05);
            close(s.mass_kg(), 0.001, 1e-15);
            close(a.incoming_energy_j, 400.0, 1e-10);
            final_t.push(s.temperature_k(CellId(0)).unwrap());
        }
        close(final_t[0], final_t[1], 1e-10);
        close(final_t[0], final_t[2], 1e-10);
    }

    // G1: insulated 1-D strip, cell-centred cosine mode. Unlike the uniform
    // partition test this exercises a spatial temperature gradient. Conductance
    // is k*A/dx, mass is rho*A*dx; continuum alpha=k/(rho*cp).
    #[test]
    fn nonuniform_strip_spatial_refinement_matches_diffusion_mode() {
        let length_m = 0.01;
        let cross_section_m2 = 1e-6;
        let density_kg_m3 = 1000.0;
        let conductivity_w_m_k = 0.2;
        let cp_j_kg_k = 2000.0;
        let end_s = 10.0;
        let alpha = conductivity_w_m_k / (density_kg_m3 * cp_j_kg_k);
        let mut errors = Vec::new();
        for count in [4, 8, 16] {
            let dx = length_m / f64::from(count);
            let mut events = Vec::new();
            for i in 0..count {
                let x = (f64::from(i) + 0.5) * dx;
                let initial_t = 400.0 + 50.0 * (std::f64::consts::PI * x / length_m).cos();
                let mut c = cell(i as u64, cp_j_kg_k * (initial_t - 300.0));
                c.mass_kg = density_kg_m3 * cross_section_m2 * dx;
                c.geometry.start_m = [f64::from(i) * dx, 0.0, 0.0];
                c.geometry.end_m = [(f64::from(i) + 1.0) * dx, 0.0, 0.0];
                c.geometry.volume_m3 = cross_section_m2 * dx;
                c.geometry.surface_area_m2 = 2.0 * cross_section_m2 + 0.004 * dx;
                events.push(event(0.0, Event::Birth(c)));
            }
            for i in 1..count {
                let mut l = link((i - 1) as u64, i as u64);
                l.area_m2 = cross_section_m2;
                l.conductance_w_k = conductivity_w_m_k * cross_section_m2 / dx;
                events.push(event(0.0, Event::Contact(l)));
            }
            let (s, a) = run(&simulation(events.clone()), end_s, 0.001);
            let (finer_time, _) = run(&simulation(events), end_s, 0.0005);
            let amplitude =
                50.0 * (-alpha * (std::f64::consts::PI / length_m).powi(2) * end_s).exp();
            let mut error: f64 = 0.0;
            let mut time_change: f64 = 0.0;
            for i in 0..count {
                let x = (f64::from(i) + 0.5) * dx;
                let exact = 400.0 + amplitude * (std::f64::consts::PI * x / length_m).cos();
                let t = s.temperature_k(CellId(i as u64)).unwrap();
                error = error.max((t - exact).abs());
                time_change = time_change
                    .max((t - finer_time.temperature_k(CellId(i as u64)).unwrap()).abs());
            }
            assert!(
                time_change < error * 0.01,
                "time error must not mask spatial error"
            );
            close(a.residual_j, 0.0, 1e-10);
            errors.push(error);
        }
        for pair in errors.windows(2) {
            assert!(
                (3.5..4.5).contains(&(pair[0] / pair[1])),
                "errors={errors:?}"
            );
        }
    }

    // An independently specified piecewise enthalpy plateau checks that the
    // network carries h, rather than replacing it by temperature/cp storage.
    #[test]
    fn latent_plateau_stores_energy_before_temperature_rises() {
        struct Plateau;
        impl ThermalLaw for Plateau {
            fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
                if !(0.0..=400_000.0).contains(&h) {
                    return Err(DepositionError::MaterialDomain);
                }
                Ok(if h < 100_000.0 {
                    300.0 + h / 2000.0
                } else if h <= 200_000.0 {
                    350.0
                } else {
                    350.0 + (h - 200_000.0) / 2000.0
                })
            }
            fn maximum_temperature_slope(&self) -> f64 {
                0.0005
            }
        }
        let law = Plateau;
        let mut c = cell(1, 100_000.0);
        c.boundary.source_w = 10.0;
        let s = DepositionSimulation::new(
            vec![event(0.0, Event::Birth(c))],
            &[(MaterialId(7), &law)],
            budget(),
        )
        .unwrap();
        let (middle, a) = s.advance_to(5.0, config(0.125), || Ok(())).unwrap();
        close(middle.temperature_k(CellId(1)).unwrap(), 350.0, 1e-12);
        close(middle.cells()[&CellId(1)].specific_h_j_kg, 150_000.0, 1e-8);
        close(a.external_heat_j, 50.0, 1e-10);
        let (after, a) = middle.advance_to(15.0, config(0.125), || Ok(())).unwrap();
        close(after.temperature_k(CellId(1)).unwrap(), 375.0, 1e-10);
        close(a.residual_j, 0.0, 1e-10);
    }

    #[test]
    fn cancelled_trial_and_exhausted_work_do_not_commit_births() {
        let sim = simulation(vec![birth(0.0, 1, 0.0), birth(0.5, 2, 400_000.0)]);
        let mut calls = 0;
        let err = sim.advance_to(1.0, config(0.1), || {
            calls += 1;
            if calls > 15 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(err, Err(DepositionError::Cancelled)));
        assert_eq!(sim.time_s(), 0.0);
        assert!(sim.cells().is_empty());
        let (restarted, a) = run(&sim, 1.0, 0.1);
        assert_eq!(a.events, 2);
        close(restarted.mass_kg(), 0.002, 1e-15);
        let mut limited = sim.clone();
        limited.budget.max_substeps = 1;
        assert!(matches!(
            limited.advance_to(1.0, config(0.1), || Ok(())),
            Err(DepositionError::Budget("thermal substeps"))
        ));
        assert!(limited.cells().is_empty());
    }

    #[test]
    fn checkpoint_restart_replays_identical_step_grid() {
        let sim = simulation(vec![
            birth(0.0, 1, 40_000.0),
            birth(1.0, 2, 400_000.0),
            event(1.0, Event::Contact(link(1, 2))),
        ]);
        let (whole, _) = run(&sim, 2.0, 0.125);
        let (mid, _) = run(&sim, 1.0, 0.125);
        let checkpoint = mid.clone();
        let (resumed, _) = run(&checkpoint, 2.0, 0.125);
        assert_eq!(whole.cells(), resumed.cells());
        assert_eq!(whole.links(), resumed.links());
    }

    #[test]
    fn invalid_topology_mass_units_material_and_history_are_refused() {
        for events in [
            vec![birth(1.0, 1, 0.0), birth(0.0, 2, 0.0)],
            vec![birth(0.0, 1, 0.0), birth(0.0, 1, 0.0)],
            vec![birth(0.0, 1, 0.0), event(0.0, Event::Contact(link(1, 2)))],
            vec![birth(0.0, 1, 2_000_000.0)],
        ] {
            match DepositionSimulation::new(events, &[(MaterialId(7), &LAW)], budget()) {
                Err(_) => {}
                Ok(s) => {
                    assert!(s.advance_to(2.0, config(0.1), || Ok(())).is_err());
                    assert!(s.cells().is_empty());
                }
            }
        }
        for mut c in [cell(1, 0.0), cell(2, 0.0), cell(3, 0.0)] {
            match c.id.0 {
                1 => c.mass_kg = 0.0,
                2 => c.material_id = MaterialId(999),
                _ => c.geometry.volume_m3 = f64::NAN,
            }
            assert!(
                simulation(vec![event(0.0, Event::Birth(c))])
                    .advance_to(0.0, config(1.0), || Ok(()))
                    .is_err()
            );
        }
        let s = simulation(vec![birth(0.0, 1, 0.0)]);
        assert!(s.advance_to(-1.0, config(1.0), || Ok(())).is_err());
        assert!(s.advance_to(1.0, config(0.0), || Ok(())).is_err());
    }

    #[test]
    fn overlapping_contact_area_is_refused_even_when_insulating() {
        let mut l = link(1, 2);
        l.area_m2 = 0.002;
        l.conductance_w_k = 0.0;
        let sim = simulation(vec![
            birth(0.0, 1, 0.0),
            birth(0.0, 2, 0.0),
            event(0.0, Event::Contact(l)),
        ]);
        assert!(matches!(
            sim.advance_to(0.0, config(1.0), || Ok(())),
            Err(DepositionError::InvalidInput(
                "contacts exceed surface area"
            ))
        ));
    }

    #[test]
    fn heterogeneous_material_identity_and_source_energy_are_retained() {
        struct Different;
        impl ThermalLaw for Different {
            fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
                Ok(300.0 + h / 1000.0)
            }
            fn maximum_temperature_slope(&self) -> f64 {
                0.001
            }
        }
        let law = Different;
        let mut c = cell(2, 200_000.0);
        c.material_id = MaterialId(9);
        c.boundary.source_w = 2.0;
        let s = DepositionSimulation::new(
            vec![birth(0.0, 1, 400_000.0), event(0.0, Event::Birth(c))],
            &[(MaterialId(7), &LAW), (MaterialId(9), &law)],
            budget(),
        )
        .unwrap();
        let (s, a) = s.advance_to(2.0, config(0.1), || Ok(())).unwrap();
        close(a.external_heat_j, 4.0, 1e-12);
        close(s.energy_j(), 604.0, 1e-9);
        assert_eq!(s.cells()[&CellId(2)].material_id, MaterialId(9));
    }

    /// Cost observation, not a hardware-independent performance acceptance.
    #[test]
    #[ignore = "run optimized with --ignored --nocapture for schedule cost"]
    fn representative_layered_strip_cost() {
        let mut events = Vec::new();
        for id in 0..1000 {
            events.push(birth(id as f64 * 0.001, id, 400_000.0));
            if id > 0 {
                events.push(event(id as f64 * 0.001, Event::Contact(link(id - 1, id))));
            }
        }
        let start = std::time::Instant::now();
        let (s, a) = run(&simulation(events), 2.0, 0.001);
        close(s.mass_kg(), 1.0, 1e-12);
        close(s.energy_j(), 400_000.0, 1e-8);
        println!(
            "1000 cells, {} links, {} steps, {} events: {:?}; residual={} J",
            s.links().len(),
            a.substeps,
            a.events,
            start.elapsed(),
            a.residual_j
        );
    }
}
