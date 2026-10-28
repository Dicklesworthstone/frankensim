//! Conservative upwind enthalpy transport on a prescribed mass-flow graph.
//!
//! This is an observer of the original enthalpy timeline, not another thermal
//! integrator. Cell mass, volume and material identity stay fixed. Each flowing
//! component must use one immutable caloric law, or explicitly admitted stationary
//! species/caloric bindings, and balance mass at every cell,
//! including explicit inlet/outlet reservoirs. Boundary enthalpy is external
//! energy; interior transfers are internal and never counted as a heat source.
//! Flux arithmetic subtracts one common enthalpy reference per connected flow
//! component. The fixed-mass balance removes that reference analytically, so
//! formation enthalpy cannot swamp sensible heating through cancellation of
//! gross inlet/outlet powers. Admitted mass-balance roundoff and its reference-
//! energy equivalent remain explicit diagnostics, never physical heat.
//! Species mixing, density/volume changes, pressure work and phase-dependent
//! mobility need their own models; no velocity or melt motion is inferred here.

/// Stationary species balances on the same prescribed flow graph, including
/// bounded inlet-fed recirculation.
// enthalpy.rs loads this file by #[path], so a bare `mod` would resolve beside it.
#[path = "advection/species.rs"]
pub mod species;

/// Opt-in entropy transport on the same accepted upwind trajectory.
#[path = "advection/entropy.rs"]
pub mod entropy;

use std::collections::{BTreeMap, BTreeSet};
use super::{CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, HeatTransfer, MaterialId};

/// One physical connection, independent of conductive-contact identity/area.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassFlow {
    /// Stable connection identity. Parallel connections need distinct IDs.
    pub id: u64,
    /// First active thermal control volume.
    pub a: CellId,
    /// Second active thermal control volume.
    pub b: CellId,
    /// Signed mass rate [kg/s], positive from a to b.
    pub kg_s: f64,
}

impl MassFlow {
    // One direction convention shared by enthalpy and stationary species.
    fn directed(self) -> Option<(CellId, CellId, f64)> {
        if self.kg_s == 0.0 { None }
        else if self.kg_s > 0.0 { Some((self.a, self.b, self.kg_s)) }
        else { Some((self.b, self.a, -self.kg_s)) }
    }
}

/// External mass connection. Rates are positive magnitudes [kg/s]. Inlet
/// material identity and enthalpy datum must match the receiving cell's law.
/// Outlets always carry the CURRENT local enthalpy, never a prescribed value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MassBoundary {
    /// Supplied material entering an active control volume.
    Inlet {
        /// Stable connection identity, distinct from all interior/boundary IDs.
        id: u64,
        /// Receiving control volume.
        cell: CellId,
        /// Positive inward mass rate [kg/s].
        kg_s: f64,
        /// Caloric law and energy-datum identity.
        material_id: MaterialId,
        /// Incoming specific enthalpy [J/kg], within that law's domain.
        specific_h_j_kg: f64,
    },
    /// Material leaving an active control volume at its current enthalpy.
    Outlet {
        /// Stable connection identity.
        id: u64,
        /// Outgoing control volume.
        cell: CellId,
        /// Positive outward mass rate [kg/s].
        kg_s: f64,
    },
}
impl MassBoundary {
    fn identity(self) -> (u64, CellId, f64) {
        match self {
            Self::Inlet { id, cell, kg_s, .. } | Self::Outlet { id, cell, kg_s } => (id, cell, kg_s),
        }
    }
}

/// Numerical observations for the last successful advance. Boundary energy
/// duplicates an observation of the thermal audit, never an additional source.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AdvectionAudit {
    /// Sum of absolute mass transported across all connections [kg]. A mass
    /// parcel may cross several connections; this is NOT added/stored mass.
    pub transported_mass_kg: f64,
    /// Sum of absolute connection enthalpy transfers [J]. This depends on the
    /// common enthalpy datum; it is NOT irreversible dissipation or external heat.
    pub absolute_enthalpy_transfer_j: f64,
    /// Gross externally supplied mass [kg]; not added to fixed stored mass.
    pub inlet_mass_kg: f64,
    /// Gross externally discharged mass [kg].
    pub outlet_mass_kg: f64,
    /// Net inlet minus outlet enthalpy [J], already included exactly once in
    /// the thermal external-input account. Do not add it to that audit again.
    pub external_enthalpy_j: f64,
    /// Integral of the sum of absolute cell incoming-minus-outgoing mass rates
    /// admitted by the fixed-mass roundoff test [kg]. The evolution enforces
    /// zero continuity; this is numerical input roundoff, not stored/leaked mass.
    pub absolute_mass_balance_roundoff_kg: f64,
    /// Integral of sum |h_reference * (incoming-outgoing mass rate)| [J].
    /// Reports the coordinate-dependent term removed by that zero-continuity
    /// convention. It is not external heat or an extra energy debit.
    pub absolute_datum_energy_roundoff_j: f64,
    /// Largest outgoing mass / stored cell mass on an accepted substep.
    pub maximum_courant: f64,
    /// Actual accepted thermal substeps carrying the prescribed flow graph.
    pub substeps: usize,
}

#[derive(Clone)]
struct Pending {
    time_s: f64,
    decay_per_s: f64,
    mass_kg_s: f64,
    absolute_power_w: f64,
    inlet_kg_s: f64,
    outlet_kg_s: f64,
    external_power_w: f64,
    mass_balance_roundoff_kg_s: f64,
    datum_energy_roundoff_w: f64,
}

#[derive(Clone, Default)]
struct BoundaryPending {
    time_s: f64,
    rates: BTreeMap<CellId, [f64; 2]>,
    inlet_kg_s: f64,
    outlet_kg_s: f64,
    external_power_w: f64,
    absolute_power_w: f64,
}

/// Add prescribed flow to ANY existing thermal source/history observer.
///
/// The wrapped observer still owns its sources, conductances, events and work
/// accounts. Advection contributes pairwise INTERNAL watts, explicit reservoir
/// EXTERNAL watts, and a combined CFL restriction to that very same step. Material T(h), including sensible-heat
/// curvature and latent plateaux, is always evaluated by the original timeline.
/// Recreate this wrapper at an accepted endpoint to change the prescribed flow.
#[derive(Clone)]
pub struct Advected<O> {
    inner: O,
    flows: Vec<MassFlow>,
    boundaries: Vec<MassBoundary>,
    // Stable component representatives, compiled once for nonzero connections.
    // Boundary-only cells use themselves and need no topology entry.
    reference_cells: BTreeMap<CellId, CellId>,
    stationary_species: Option<species::StationarySpecies>,
    boundary_pending: Option<BoundaryPending>,
    pending: Option<Pending>,
    audit: AdvectionAudit,
    entropy: Option<entropy::EntropyState>,
}

fn finite(x: f64) -> Result<f64, DepositionError> {
    if x.is_finite() { Ok(x) } else { Err(DepositionError::Arithmetic) }
}
fn invalid(s: &'static str) -> DepositionError { DepositionError::InvalidInput(s) }

fn component_references(flows: &[MassFlow]) -> BTreeMap<CellId, CellId> {
    fn root(parents: &mut BTreeMap<CellId, CellId>, mut id: CellId) -> CellId {
        loop {
            let next = parents[&id];
            if next == id { return id; }
            let grandparent = parents[&next];
            parents.insert(id, grandparent);
            id = next;
        }
    }
    let mut parents = BTreeMap::new();
    for flow in flows.iter().filter(|flow| flow.kg_s != 0.0) {
        parents.entry(flow.a).or_insert(flow.a);
        parents.entry(flow.b).or_insert(flow.b);
        let a = root(&mut parents, flow.a);
        let b = root(&mut parents, flow.b);
        parents.insert(a.max(b), a.min(b));
    }
    let ids: Vec<_> = parents.keys().copied().collect();
    for id in ids {
        let representative = root(&mut parents, id);
        parents.insert(id, representative);
    }
    parents
}

impl<O> Advected<O> {
    /// Admit bounded, uniquely named connections. Sort by stable ID so input
    /// order does not choose floating-point reduction order. Active-cell, common
    /// material and mass-balance admission happens at each actual thermal step.
    pub fn new(inner: O, mut flows: Vec<MassFlow>, max_flows: usize) -> Result<Self, DepositionError> {
        if flows.len() > max_flows { return Err(DepositionError::Budget("advective connections")); }
        let mut ids = BTreeSet::new();
        for flow in &flows {
            if flow.a == flow.b || !flow.kg_s.is_finite() || !ids.insert(flow.id) {
                return Err(invalid("invalid or duplicate advective connection"));
            }
        }
        flows.sort_by_key(|f| f.id);
        let reference_cells = component_references(&flows);
        Ok(Self { inner, flows, boundaries: Vec::new(), reference_cells,
            stationary_species: None, boundary_pending: None,
            pending: None, audit: AdvectionAudit::default(), entropy: None })
    }
    /// Replace the complete inlet/outlet data at an accepted endpoint. The
    /// declared ceiling covers interior AND boundary connections. Cell-level
    /// fixed-mass balance and inlet material support are checked on each step.
    /// This does not infer reservoir pressure, velocity, mixing or momentum.
    pub fn with_boundaries(mut self, mut boundaries: Vec<MassBoundary>, max_connections: usize)
        -> Result<Self, DepositionError> {
        if self.flows.len().checked_add(boundaries.len()).is_none_or(|n| n > max_connections) {
            return Err(DepositionError::Budget("mass-flow connections"));
        }
        let mut ids: BTreeSet<_> = self.flows.iter().map(|flow| flow.id).collect();
        for boundary in &boundaries {
            let (id, _, kg_s) = boundary.identity();
            if !kg_s.is_finite() || kg_s <= 0.0 || !ids.insert(id) {
                return Err(invalid("invalid or duplicate mass-flow boundary"));
            }
            if let MassBoundary::Inlet { specific_h_j_kg, .. } = boundary {
                finite(*specific_h_j_kg)?;
            }
        }
        boundaries.sort_by_key(|boundary| boundary.identity().0);
        self.boundaries = boundaries;
        self.stationary_species = None;
        self.boundary_pending = None;
        self.pending = None;
        Ok(self)
    }
    /// Immutable inlet/outlet data in stable connection-ID order.
    pub fn boundaries(&self) -> &[MassBoundary] { &self.boundaries }
    /// Original source/history owner, including any accepted material receipts.
    pub fn inner(&self) -> &O { &self.inner }
    /// Consume the wrapper without erasing the original source/history type.
    pub fn into_inner(self) -> O { self.inner }
    /// Immutable prescribed connection data in stable-ID order.
    pub fn flows(&self) -> &[MassFlow] { &self.flows }
    /// Observations reset by begin(), on the tentative observer only.
    pub fn audit(&self) -> AdvectionAudit { self.audit }

    fn reference_h(&self, thermal: &DepositionSimulation<'_>, id: CellId)
        -> Result<f64, DepositionError> {
        let reference = self.reference_cells.get(&id).unwrap_or(&id);
        thermal.cells().get(reference).map(|cell| cell.specific_h_j_kg)
            .ok_or(invalid("inactive advective reference cell"))
    }
}

impl<O: DepositionObserver> DepositionObserver for Advected<O> {
    type Error = O::Error;

    fn requires_step_interval(&self) -> bool { self.inner.requires_step_interval() }

    fn step_sources(&mut self, thermal_start: &DepositionSimulation<'_>,
        start_s: f64, dt_s: f64, end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner.step_sources(thermal_start, start_s, dt_s, end_s, check)
    }

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner.begin()?;
        self.pending = None;
        self.boundary_pending = None;
        self.audit = AdvectionAudit::default();
        if let Some(entropy) = &mut self.entropy { entropy.reset(); }
        Ok(())
    }
    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        if self.entropy.is_some() && self.stationary_species.is_some() {
            return Err(invalid("thermal advection entropy does not cover stationary species mixing").into());
        }
        if let Some(entropy) = &mut self.entropy {
            entropy.capture(thermal, check)?;
        }
        if let Some(species) = &self.stationary_species {
            species.validate_laws(thermal, &self.boundaries)?;
        }
        let mut sources = self.inner.sources(thermal, check)?;
        let mut pending = BoundaryPending { time_s: thermal.time_s(), ..BoundaryPending::default() };
        for (index, boundary) in self.boundaries.iter().enumerate() {
            if index % 256 == 0 { check()?; }
            let (_, id, kg_s) = boundary.identity();
            let cell = thermal.cells().get(&id).ok_or(invalid("inactive mass-flow boundary cell"))?;
            let reference_h = self.reference_h(thermal, id)?;
            let (power_w, gross_power_w, side) = match *boundary {
                MassBoundary::Inlet { material_id, specific_h_j_kg, .. } => {
                    if material_id != cell.material_id && self.stationary_species.is_none() {
                        return Err(invalid("inlet material/datum does not match receiving cell").into());
                    }
                    // Stationary mixing validates a common reference, but the
                    // inlet still carries its OWN species-weighted enthalpy.
                    let law = thermal.materials.get(&material_id)
                        .ok_or(invalid("unregistered inlet caloric law"))?;
                    let t = law.temperature_k(specific_h_j_kg)?;
                    if !t.is_finite() || t <= 0.0 { return Err(DepositionError::MaterialDomain.into()); }
                    pending.inlet_kg_s = finite(pending.inlet_kg_s + kg_s)?;
                    (finite(kg_s * finite(specific_h_j_kg - reference_h)?)?,
                        finite(kg_s * specific_h_j_kg)?, 0)
                }
                MassBoundary::Outlet { .. } => {
                    pending.outlet_kg_s = finite(pending.outlet_kg_s + kg_s)?;
                    (-finite(kg_s * finite(cell.specific_h_j_kg - reference_h)?)?,
                        -finite(kg_s * cell.specific_h_j_kg)?, 1)
                }
            };
            let rates = pending.rates.entry(id).or_default();
            rates[side] = finite(rates[side] + kg_s)?;
            pending.external_power_w = finite(pending.external_power_w + power_w)?;
            pending.absolute_power_w = finite(pending.absolute_power_w + gross_power_w.abs())?;
            if power_w != 0.0 {
                let watts = sources.entry(id).or_default();
                *watts = finite(*watts + power_w)?;
            }
        }
        check()?;
        self.boundary_pending = Some(pending);
        Ok(sources)
    }
    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner.conductances(thermal, check)
    }
    fn heat_transfers(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<Vec<HeatTransfer>, Self::Error> {
        check()?;
        let mut transfers = self.inner.heat_transfers(thermal, check)?;
        // Separate incoming/outgoing sums avoid cancellation in the admission
        // scale. Small roundoff is admitted, not an adjustable mass-leak budget.
        let boundary = self.boundary_pending.take().ok_or(invalid("mass-flow boundaries have no source epoch"))?;
        if boundary.time_s != thermal.time_s() {
            return Err(invalid("mass-flow boundary epoch does not match thermal state").into());
        }
        let mut rates = boundary.rates;
        let mut mass_kg_s = finite(boundary.inlet_kg_s + boundary.outlet_kg_s)?;
        let mut absolute_power_w = boundary.absolute_power_w;
        for (index, flow) in self.flows.iter().enumerate() {
            if index % 256 == 0 { check()?; }
            let a = thermal.cells().get(&flow.a).ok_or(invalid("inactive advective cell"))?;
            let b = thermal.cells().get(&flow.b).ok_or(invalid("inactive advective cell"))?;
            let Some((up_id, _, rate)) = flow.directed() else { continue; };
            if a.material_id != b.material_id && self.stationary_species.is_none() {
                return Err(invalid("advection across different material IDs requires a composition model").into());
            }
            let (up, down) = if up_id == a.id { (a, b) } else { (b, a) };
            let reference_h = self.reference_h(thermal, up.id)?;
            let power_w = finite(rate * finite(up.specific_h_j_kg - reference_h)?)?;
            let gross_power_w = finite(rate * up.specific_h_j_kg)?;
            let out = rates.entry(up.id).or_default();
            out[1] = finite(out[1] + rate)?;
            let incoming = rates.entry(down.id).or_default();
            incoming[0] = finite(incoming[0] + rate)?;
            mass_kg_s = finite(mass_kg_s + rate)?;
            absolute_power_w = finite(absolute_power_w + gross_power_w.abs())?;
            transfers.push(HeatTransfer { a: up.id, b: down.id, power_w });
        }
        let mut decay_per_s: f64 = 0.0;
        let mut mass_balance_roundoff_kg_s = 0.0;
        let mut datum_energy_roundoff_w = 0.0;
        for (index, (&id, &[incoming, outgoing])) in rates.iter().enumerate() {
            if index % 256 == 0 { check()?; }
            let scale = incoming.max(outgoing);
            if (incoming-outgoing).abs() > 64.0*f64::EPSILON*scale {
                return Err(invalid("fixed-mass advection requires balanced incoming and outgoing mass at every cell").into());
            }
            let mass_roundoff = (incoming-outgoing).abs();
            mass_balance_roundoff_kg_s = finite(mass_balance_roundoff_kg_s + mass_roundoff)?;
            datum_energy_roundoff_w = finite(datum_energy_roundoff_w
                + finite(self.reference_h(thermal, id)? * mass_roundoff)?.abs())?;
            let decay = finite(outgoing / thermal.cells()[&id].mass_kg)?;
            if outgoing != 0.0 && decay == 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            decay_per_s = decay_per_s.max(decay);
        }
        check()?;
        self.pending = Some(Pending { time_s: thermal.time_s(), decay_per_s,
            mass_kg_s, absolute_power_w, inlet_kg_s: boundary.inlet_kg_s,
            outlet_kg_s: boundary.outlet_kg_s, external_power_w: boundary.external_power_w,
            mass_balance_roundoff_kg_s, datum_energy_roundoff_w });
        Ok(transfers)
    }
    fn limit_timestep_s(&mut self, thermal: &DepositionSimulation<'_>, proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<f64, Self::Error> {
        let dt = self.inner.limit_timestep_s(thermal, proposed_s, check)?;
        if !dt.is_finite() || dt <= 0.0 || dt > proposed_s {
            return Err(invalid("wrapped observer returned an invalid timestep bound").into());
        }
        let pending = self.pending.as_ref().ok_or(invalid("advection has no prepared source epoch"))?;
        if pending.time_s != thermal.time_s() {
            return Err(invalid("advection source epoch does not match thermal state").into());
        }
        if pending.decay_per_s == 0.0 { return Ok(dt); }
        // dt <= 0.8/a from the existing conductive kernel; b is the largest
        // advective outgoing mass/mass. The harmonic limit obeys dt <= 0.8/(a+b),
        // unlike independently imposing two CFL limits (which can overshoot).
        // Stable reciprocal evaluation avoids overflowing dt*b or 1/tiny_dt.
        let advective_dt = 0.8 / pending.decay_per_s;
        let result = if dt <= advective_dt {
            dt / (1.0 + dt/advective_dt)
        } else {
            advective_dt / (1.0 + advective_dt/dt)
        };
        if !result.is_finite() || result <= 0.0 {
            return Err(DepositionError::Arithmetic.into());
        }
        Ok(result)
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.inner.event(event, thermal)
    }
    fn step(&mut self, start_s: f64, dt: f64, old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>, check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error> {
        check()?;
        let pending = self.pending.take().ok_or(invalid("advection step has no prepared flow"))?;
        if pending.time_s != start_s || !dt.is_finite() || dt <= 0.0 {
            return Err(invalid("advection step/source epoch mismatch").into());
        }
        self.audit.transported_mass_kg = finite(self.audit.transported_mass_kg + finite(dt*pending.mass_kg_s)?)?;
        self.audit.absolute_enthalpy_transfer_j = finite(self.audit.absolute_enthalpy_transfer_j
            + finite(dt*pending.absolute_power_w)?)?;
        self.audit.inlet_mass_kg = finite(self.audit.inlet_mass_kg + finite(dt*pending.inlet_kg_s)?)?;
        self.audit.outlet_mass_kg = finite(self.audit.outlet_mass_kg + finite(dt*pending.outlet_kg_s)?)?;
        self.audit.external_enthalpy_j = finite(self.audit.external_enthalpy_j + finite(dt*pending.external_power_w)?)?;
        self.audit.absolute_mass_balance_roundoff_kg = finite(self.audit.absolute_mass_balance_roundoff_kg
            + finite(dt*pending.mass_balance_roundoff_kg_s)?)?;
        self.audit.absolute_datum_energy_roundoff_j = finite(self.audit.absolute_datum_energy_roundoff_j
            + finite(dt*pending.datum_energy_roundoff_w)?)?;
        self.audit.maximum_courant = self.audit.maximum_courant.max(finite(dt*pending.decay_per_s)?);
        self.audit.substeps = self.audit.substeps.checked_add(1).ok_or(DepositionError::Arithmetic)?;
        if let Some(entropy) = &mut self.entropy {
            entropy.record(thermal, &self.flows, &self.boundaries, &self.reference_cells, old, dt, check)?;
        }
        self.inner.step(start_s, dt, old, thermal, check)
    }
}

#[cfg(test)]
#[path = "advection/tests.rs"]
mod tests;
