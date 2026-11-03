//! Staggered mechanical separation/reclosure coupled to real thermal contact.
//!
//! A bonded original patch retains its supplied thermal conductance until the
//! cohesive law fails. A failed patch opens at a declared normal-gap threshold
//! and recontacts at a separate closure threshold. Closed failed patches use a
//! separately supplied conductance; recontact NEVER resets cohesive damage.
//! Partial damage is not interpreted as lost geometric area or scaled heat flow.
//!
//! Each accepted mechanical endpoint changes the NEXT thermal interval. This
//! does not localize a fracture time inside the preceding step, solve mechanical
//! equilibrium or move geometry. `synchronize_with_heat` optionally deposits
//! an explicit fraction of newly dissipated fracture work into the same cells.
//! The timed synchronization methods also consume physical bond viscosity from
//! the original constant-rate path, with separate fracture/viscous heat fractions
//! and a retained mechanical clock independent of thermal integration stops.
//! Fixed local frames,
//! small tangential travel and uniform original contact area are explicit limits.
//! The opt-in [`prescribed`] path instead resolves switching times inside an
//! explicitly supplied straight separation segment before advancing its heat.

/// Event-resolved heat transfer for explicitly prescribed separation histories.
pub mod prescribed;

use std::collections::BTreeMap;
use fs_exec::Cx;
use fs_material::cohesive::{CohesiveError, CohesivePoint};
use crate::transient::deposition::{CellId, ContactUpdate, DepositionError,
    DepositionObserver, DepositionSimulation, EnergyAudit, Event, StepConfig, ThermalLink};

/// Explicit reduced closure model, not a measured universal adhesion relation.
#[derive(Debug, Clone, Copy)]
pub struct ContactPolicy {
    /// A failed patch opens at or above this positive normal gap, m.
    pub open_gap_m: f64,
    /// An open failed patch recloses at or below this gap in [0, open_gap), m.
    pub reclose_gap_m: f64,
    /// Independent total W/K after debonding while still/re-again closed.
    pub reclosed_conductance_w_k: f64,
    /// Fixed-area/frame applicability limit for the local tangential jump, m.
    pub maximum_slip_m: f64,
}

/// One interface at the initial accepted mechanical/thermal stage. The point
/// may already carry damage; it is retained, never replaced by a virgin point.
#[derive(Debug, Clone, Copy)]
pub struct Interface<'a> {
    /// Original area/endpoints and the separately supplied bonded conductance.
    pub reference: ThermalLink,
    pub point: CohesivePoint<'a>,
    pub closed: bool,
    pub policy: ContactPolicy,
}
impl Interface<'_> {
    fn next_closed(&self, point: CohesivePoint<'_>) -> Result<bool, Error> {
        let p = self.policy;
        let j = point.jump_m();
        // Normalize before squaring to avoid under/overflow of metre squares.
        let a = j[1].abs() / p.maximum_slip_m;
        let b = j[2].abs() / p.maximum_slip_m;
        if !a.is_finite() || !b.is_finite() || a*a+b*b > 1.0 {
            return Err(Error::Invalid("tangential travel exceeds the fixed-contact-area model"));
        }
        Ok(if !point.failed() { true }
            else if j[0] >= p.open_gap_m { false }
            else if j[0] <= p.reclose_gap_m { true }
            else { self.closed })
    }

    fn thermal_link(&self) -> Option<ThermalLink> {
        self.closed.then(|| ThermalLink {
            conductance_w_k: if self.point.failed() { self.policy.reclosed_conductance_w_k }
                else { self.reference.conductance_w_k },
            ..self.reference
        })
    }
}

#[derive(Debug, Clone)]
pub enum Error {
    Thermal(DepositionError),
    Material { interface: u64, error: CohesiveError },
    Invalid(&'static str),
}
impl From<DepositionError> for Error {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fracture/thermal contact refused: {self:?}")
    }
}
impl std::error::Error for Error {}

/// Observations for this endpoint, not a closed thermodynamic balance.
#[derive(Debug, Clone, Copy, Default)]
pub struct ContactReport {
    pub opened: usize,
    pub reclosed: usize,
    pub changed_conductances: usize,
    /// Incremental fracture dissipation over the bound original areas, J.
    /// Reported only: not automatically added to thermal enthalpy.
    pub fracture_dissipation_j: f64,
    /// Physical surviving-bond dashpot loss on the timed straight jump path, J.
    /// Separate from intrinsic fracture; zero for untimed synchronization.
    pub viscous_dissipation_j: f64,
}

/// Explicit partition of NEW physical cohesive dissipation on one interface.
/// Unconverted work remains a separately reported nonthermal fracture sink.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FractureHeatPartition {
    pub interface: u64,
    /// Fraction of new fracture work converted to immediate heat, in [0, 1].
    pub fraction_to_heat: f64,
    /// Fraction of that heat entering `reference.a`, in [0, 1]. The remainder
    /// enters `reference.b`, even when the mechanical patch has just opened.
    pub fraction_to_first: f64,
}

/// Independent conversion and receiving-side splits for the two physical losses.
/// Fracture energy can remain in newly created surfaces while bond viscosity
/// follows a different explicitly supplied heat conversion. No default is used.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CohesiveHeatPartition {
    pub interface: u64,
    pub fracture_fraction_to_heat: f64,
    pub fracture_fraction_to_first: f64,
    pub viscous_fraction_to_heat: f64,
    pub viscous_fraction_to_first: f64,
}

/// One jointly accepted contact/history endpoint and its physical heat input.
#[derive(Debug, Clone)]
pub struct FractureHeatReport {
    pub contact: ContactReport,
    /// Actual supplied heat [J] at the original cells, in stable ID order.
    pub heat_by_cell_j: BTreeMap<CellId, f64>,
    /// Fracture work not converted into heat [J]; never numerical dissipation.
    pub nonthermal_fracture_energy_j: f64,
    /// Physical bond-viscosity loss not converted by the explicit partition, J.
    pub nonthermal_viscous_energy_j: f64,
    /// The original thermal owner's storage/input account. `external_heat_j`
    /// is internal to the combined mechanical/thermal system.
    pub thermal: EnergyAudit,
    /// Reported thermal storage + reported nonthermal energy - reported
    /// fracture and viscous work [J], resolved through compensated differences.
    pub transfer_residual_j: f64,
}

#[derive(Clone)]
struct Contacts<'a>(BTreeMap<u64, Interface<'a>>);
impl DepositionObserver for Contacts<'_> {
    type Error = Error;
    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Error> {
        let id = match event { Event::Contact(link) => Some(link.id),
            Event::Separation(id) => Some(*id), _ => None };
        if id.is_some_and(|id| self.0.contains_key(&id)) {
            return Err(Error::Invalid("scheduled contact change conflicts with mechanical ownership"));
        }
        Ok(())
    }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>, _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Error> {
        checkpoint()?;
        Ok(())
    }
}

/// Accepted thermal state and mechanical contact endpoints. Thermal transport
/// uses the previous accepted closure, then `synchronize` changes the following
/// interval. A caller stages a structural solve and commits that proposal only
/// after synchronization succeeds. Arbitrary material-law changes are refused.
#[derive(Clone)]
pub struct FractureCooling<'t, 'm> {
    thermal: DepositionSimulation<'t>,
    contacts: Contacts<'m>,
    // Thermal transport can advance before the caller supplies its mechanical
    // endpoint. Keep the actual interval separately; neither thermal substeps
    // nor a failed synchronization may reset its starting time.
    last_synchronized_time_s: f64,
}
impl<'t, 'm> FractureCooling<'t, 'm> {
    /// Bind already activated thermal cells and exactly matching current links.
    /// Already-open/damaged interfaces may be imported with their real history.
    /// Construction and all sweeps are bounded by `max_interfaces`; cloning is
    /// bounded by the thermal owner's budget, but is not internally cancellable.
    pub fn new(thermal: DepositionSimulation<'t>, interfaces: Vec<Interface<'m>>,
        max_interfaces: usize, mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        if interfaces.is_empty() || interfaces.len() > max_interfaces {
            return Err(Error::Invalid("interface count is empty or exceeds the declared budget"));
        }
        let mut contacts = BTreeMap::new();
        let mut updates = Vec::with_capacity(interfaces.len());
        for interface in interfaces {
            checkpoint()?;
            let link = interface.reference;
            let p = interface.policy;
            if !positive(link.area_m2) || !nonnegative(link.conductance_w_k)
                || link.a == link.b || !thermal.cells().contains_key(&link.a)
                || !thermal.cells().contains_key(&link.b)
                || !positive(p.open_gap_m) || !nonnegative(p.reclose_gap_m)
                || p.reclose_gap_m >= p.open_gap_m
                || p.open_gap_m > interface.point.law().maximum_separation_m()
                || !nonnegative(p.reclosed_conductance_w_k) || !positive(p.maximum_slip_m)
            { return Err(Error::Invalid("invalid contact binding or explicit closure policy")); }
            if [link.a,link.b].iter().any(|id| link.area_m2 > thermal.cells()[id].geometry.surface_area_m2) {
                return Err(Error::Invalid("original contact area exceeds an incident cell surface"));
            }
            if interface.next_closed(interface.point)? != interface.closed {
                return Err(Error::Invalid("initial closure contradicts the supplied mechanical endpoint"));
            }
            let expected = interface.thermal_link();
            updates.push(ContactUpdate { id: link.id, expected, replacement: expected });
            if contacts.insert(link.id, interface).is_some() {
                return Err(Error::Invalid("duplicate mechanical interface identity"));
            }
        }
        // Also checks that due scheduled events have been processed first.
        let thermal = thermal.replace_contacts(&updates, &mut checkpoint)?;
        let last_synchronized_time_s = thermal.time_s();
        Ok(Self { thermal, contacts: Contacts(contacts), last_synchronized_time_s })
    }

    pub fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub fn interfaces(&self) -> &BTreeMap<u64, Interface<'m>> { &self.contacts.0 }
    /// Time of the retained mechanical endpoints, unchanged by heat transport.
    pub fn last_synchronized_time_s(&self) -> f64 { self.last_synchronized_time_s }

    /// Advance heat with accepted closure held fixed. A scheduled event touching
    /// any mechanically owned link refuses, even if another equal-time event
    /// would restore its old value. Unrelated births/boundary events still work.
    pub fn advance_to(&self, target_s: f64, config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, EnergyAudit), Error> {
        let (thermal, contacts, audit) = self.thermal.advance_observed(
            &self.contacts, target_s, config, checkpoint)?;
        Ok((Self { thermal, contacts, last_synchronized_time_s: self.last_synchronized_time_s }, audit))
    }

    /// Accept exactly one proposed point for every bound ID at this thermal time.
    /// Points from `fs-solid::cohesive` can be passed directly. Re-evaluating their
    /// straight jump increments against the retained history verifies the SAME
    /// damage path and immutable law; it is not a second equilibrium solve.
    /// Split paths at intermediate reversals: inconsistent retained damage
    /// refuses, but undamaged closed loops cannot be recovered from endpoints.
    ///
    /// Missing/extra/duplicate IDs, stale time, foreign laws, lost damage,
    /// thermal-contact conflicts, excessive slip, or cancellation refuse the
    /// whole proposal. Reclosing a failed interface restores only heat transfer,
    /// never a tensile/shear bond. No thermal time passes in this operation.
    pub fn synchronize(&self, time_s: f64, proposed: &[(u64, CohesivePoint<'m>)],
        local_energy_tolerance_j_m2: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, ContactReport), Error> {
        self.synchronize_recording(time_s, None, proposed, local_energy_tolerance_j_m2,
            |_, _, _, _| Ok(()), checkpoint)
    }

    /// Accept one constant-rate straight mechanical segment at the current
    /// thermal time. `duration_s` must equal the elapsed time since the last
    /// accepted synchronization and be positive; transport alone does not age
    /// the mechanical clock. The original timed cohesive owner supplies both
    /// fracture and physical viscous work, never endpoint-force quadrature.
    /// Compute `duration_s = time_s - self.last_synchronized_time_s()` from the
    /// accepted absolute times, using that same represented interval in the
    /// mechanical update; a separately rounded nominal step can differ.
    ///
    /// The caller must supply the actual segment, splitting reversals. Complete
    /// original-law history is replayed and compared before publication.
    /// Undamaged closed loops leave no point-history signature; their viscous
    /// work cannot be reconstructed or detected from the endpoints. This API
    /// neither assigns loss to heat nor solves mechanical equilibrium.
    /// Repeating the accepted time refuses.
    pub fn synchronize_over_time(&self, time_s: f64, duration_s: f64,
        proposed: &[(u64, CohesivePoint<'m>)], local_energy_tolerance_j_m2: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, ContactReport), Error> {
        self.synchronize_recording(time_s, Some(duration_s), proposed,
            local_energy_tolerance_j_m2, |_, _, _, _| Ok(()), checkpoint)
    }

    /// Accept contact/history changes and heat from their new fracture work.
    /// Exactly one explicit partition is required for every bound interface.
    /// The SAME authoritative cohesive trial supplies the work: stored elastic
    /// energy, initial damage, unloading/reclosure and numerical defects never
    /// become heat. Synchronizing an already accepted point produces no heat.
    ///
    /// The joule tolerance bounds partition/aggregation roundoff, actual cell
    /// storage defects and final report rounding together. The original T(h)
    /// law owns temperature and latent storage. All changes remain a candidate
    /// until the final checkpoint; any failure returns no partial update.
    ///
    /// This is staggered endpoint heating at fixed original areas, with no
    /// thermal time advancement or within-step feedback into immutable toughness.
    /// The caller still stages mechanical equilibrium and refines coupling times.
    pub fn synchronize_with_heat(&self, time_s: f64,
        proposed: &[(u64, CohesivePoint<'m>)], partitions: &[FractureHeatPartition],
        local_energy_tolerance_j_m2: f64, heat_energy_tolerance_j: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, FractureHeatReport), Error> {
        self.synchronize_heating(time_s, None, proposed, partitions.iter().map(|p|
            CohesiveHeatPartition { interface: p.interface,
                fracture_fraction_to_heat: p.fraction_to_heat,
                fracture_fraction_to_first: p.fraction_to_first,
                viscous_fraction_to_heat: 0.0, viscous_fraction_to_first: 0.0 }),
            local_energy_tolerance_j_m2, heat_energy_tolerance_j, checkpoint)
    }

    /// Timed counterpart with independently selected fracture and bond-viscous
    /// heat conversion. The interval obeys `synchronize_over_time`; the single
    /// authoritative replay also supplies the exact integrated dashpot loss.
    /// Each bound interface needs one complete partition, including explicit
    /// zero fractions when a channel is retained as a nonthermal sink.
    ///
    /// Contact, mechanical time/history and both caloric inputs publish together
    /// only after the original thermal law and shared joule budget accept them.
    /// No heat is invented from stored energy, numerical defects or prior damage.
    /// Temperature-dependent toughness feedback and fracture-event localization
    /// remain separate operations; this is staggered endpoint heating.
    pub fn synchronize_over_time_with_heat(&self, time_s: f64, duration_s: f64,
        proposed: &[(u64, CohesivePoint<'m>)], partitions: &[CohesiveHeatPartition],
        local_energy_tolerance_j_m2: f64, heat_energy_tolerance_j: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, FractureHeatReport), Error> {
        self.synchronize_heating(time_s, Some(duration_s), proposed, partitions.iter().copied(),
            local_energy_tolerance_j_m2, heat_energy_tolerance_j, checkpoint)
    }

    // Both public paths use the same compensated physical partition and native
    // caloric update; the optional interval selects the original material owner.
    fn synchronize_heating(&self, time_s: f64, duration_s: Option<f64>,
        proposed: &[(u64, CohesivePoint<'m>)],
        partitions: impl ExactSizeIterator<Item = CohesiveHeatPartition>,
        local_energy_tolerance_j_m2: f64, heat_energy_tolerance_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, FractureHeatReport), Error> {
        checkpoint()?;
        if partitions.len() != self.contacts.0.len() || !nonnegative(heat_energy_tolerance_j) {
            return Err(Error::Invalid("complete heat partitions and finite nonnegative joule tolerance required"));
        }
        let mut ordered = BTreeMap::new();
        for partition in partitions {
            checkpoint()?;
            if [partition.fracture_fraction_to_heat, partition.fracture_fraction_to_first,
                partition.viscous_fraction_to_heat, partition.viscous_fraction_to_first]
                .into_iter().any(|fraction| !(0.0..=1.0).contains(&fraction))
                || !self.contacts.0.contains_key(&partition.interface)
                || ordered.insert(partition.interface, partition).is_some()
            { return Err(Error::Invalid("unknown, duplicate or invalid fracture heat partition")); }
        }
        let mut heat: BTreeMap<CellId, HeatSum> = BTreeMap::new();
        let mut fracture_work = HeatSum::default();
        let mut viscous_work = HeatSum::default();
        let mut nonthermal_fracture = HeatSum::default();
        let mut nonthermal_viscous = HeatSum::default();
        let mut residual = HeatSum::default();
        let mut spent = HeatSum::default();
        let (mut next, mut contact) = self.synchronize_recording(time_s, duration_s, proposed,
            local_energy_tolerance_j_m2, |id, link, fracture, viscous| {
                let partition = ordered[&id];
                for (dissipation, fraction, split, work, nonthermal) in [
                    (fracture, partition.fracture_fraction_to_heat, partition.fracture_fraction_to_first,
                        &mut fracture_work, &mut nonthermal_fracture),
                    (viscous, partition.viscous_fraction_to_heat, partition.viscous_fraction_to_first,
                        &mut viscous_work, &mut nonthermal_viscous),
                ] {
                    let converted = fraction * dissipation;
                    let first = split * converted;
                    let second = converted - first;
                    let remainder = dissipation - converted;
                    heat.entry(link.a).or_default().add(first)?;
                    heat.entry(link.b).or_default().add(second)?;
                    work.add(dissipation)?;
                    nonthermal.add(remainder)?;
                    let mut local = HeatSum::default();
                    for term in [first, second, remainder, -dissipation] { local.add(term)?; }
                    let defect = local.value()?;
                    residual.add(defect)?;
                    spent.add(defect.abs())?;
                }
                Ok(())
            }, &mut checkpoint)?;
        let mut heat_by_cell_j = BTreeMap::new();
        for (id, sum) in heat {
            checkpoint()?;
            let (rounded, tail) = sum.rounded()?;
            heat_by_cell_j.insert(id, rounded);
            residual.add(-tail)?;
            spent.add(tail.abs())?;
        }
        // Unlike the legacy report, the opt-in heat report carries compensated
        // aggregate work. Every lost final bit spends the SAME joule budget.
        let (fracture_work, fracture_tail) = fracture_work.rounded()?;
        let (viscous_work, viscous_tail) = viscous_work.rounded()?;
        let (nonthermal_fracture_energy_j, nonthermal_fracture_tail) = nonthermal_fracture.rounded()?;
        let (nonthermal_viscous_energy_j, nonthermal_viscous_tail) = nonthermal_viscous.rounded()?;
        contact.fracture_dissipation_j = fracture_work;
        contact.viscous_dissipation_j = viscous_work;
        for term in [fracture_tail, viscous_tail, -nonthermal_fracture_tail, -nonthermal_viscous_tail] {
            residual.add(term)?;
            spent.add(term.abs())?;
        }
        let used = spent.value()?;
        if used > heat_energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(used).into());
        }
        let remaining = heat_energy_tolerance_j - used;
        let (thermal, mut audit) = next.thermal.apply_heat_j(
            &heat_by_cell_j, remaining, remaining, &mut checkpoint)?;

        // Reduce only touched cells, using exactly the original thermal owner's
        // represented m*delta_h. This is no second caloric evaluation or update.
        // Its report totals can themselves lose small terms beside large ones.
        let mut stored = HeatSum::default();
        let mut supplied = HeatSum::default();
        let mut thermal_defect = HeatSum::default();
        for (&id, &input) in &heat_by_cell_j {
            checkpoint()?;
            let old = &self.thermal.cells()[&id];
            let cell = &thermal.cells()[&id];
            let change = cell.mass_kg * (cell.specific_h_j_kg - old.specific_h_j_kg);
            let defect = change - input;
            stored.add(change)?;
            supplied.add(input)?;
            thermal_defect.add(defect)?;
            spent.add(defect.abs())?;
        }
        let (stored_change, stored_tail) = stored.rounded()?;
        let (supplied_heat, supplied_tail) = supplied.rounded()?;
        spent.add(stored_tail.abs())?;
        spent.add(supplied_tail.abs())?;
        let used = spent.value()?;
        if used > heat_energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(used).into());
        }
        // Reported thermal balance = sum(local defects) - storage tail + input
        // tail. The input-report tail cancels again in the combined transfer.
        thermal_defect.add(-stored_tail)?;
        thermal_defect.add(supplied_tail)?;
        audit.stored_energy_change_j = stored_change;
        audit.external_heat_j = supplied_heat;
        audit.residual_j = thermal_defect.value()?;
        residual.add(audit.residual_j)?;
        residual.add(-supplied_tail)?;
        let transfer_residual_j = residual.value()?;
        if transfer_residual_j.abs() > heat_energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(transfer_residual_j).into());
        }
        let report = FractureHeatReport {
            contact, heat_by_cell_j, nonthermal_fracture_energy_j, nonthermal_viscous_energy_j,
            thermal: audit, transfer_residual_j,
        };
        next.thermal = thermal;
        checkpoint()?;
        Ok((next, report))
    }

    // One authoritative trial serves both the legacy contact update and heat.
    fn synchronize_recording(&self, time_s: f64, duration_s: Option<f64>,
        proposed: &[(u64, CohesivePoint<'m>)],
        local_energy_tolerance_j_m2: f64,
        mut record: impl FnMut(u64, ThermalLink, f64, f64) -> Result<(), Error>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, ContactReport), Error> {
        checkpoint()?;
        if time_s != self.thermal.time_s() || proposed.len() != self.contacts.0.len()
            || !positive(local_energy_tolerance_j_m2)
        { return Err(Error::Invalid("require current thermal time, complete endpoints and positive energy gate")); }
        if duration_s.is_some_and(|dt| !positive(dt) || dt != time_s-self.last_synchronized_time_s) {
            return Err(Error::Invalid("timed synchronization requires the actual positive interval since the accepted mechanical endpoint"));
        }
        let mut contacts = BTreeMap::new();
        let mut updates = Vec::with_capacity(proposed.len());
        let mut report = ContactReport::default();
        // Sorted accumulation makes output independent of caller arrival order.
        let mut ordered = BTreeMap::new();
        for &(id, point) in proposed {
            checkpoint()?;
            if !self.contacts.0.contains_key(&id) || ordered.insert(id, point).is_some() {
                return Err(Error::Invalid("unknown or duplicate mechanical endpoint"));
            }
        }
        for (id, point) in ordered {
            checkpoint()?;
            let old = self.contacts.0[&id];
            if !std::ptr::eq(old.point.law(), point.law()) {
                return Err(Error::Invalid("mechanical endpoint belongs to a different immutable law"));
            }
            let (trial, viscous_density) = match duration_s {
                Some(dt) => {
                    let timed = old.point.trial_over_time(point.jump_m(), dt, local_energy_tolerance_j_m2)
                        .map_err(|error| Error::Material { interface: id, error })?;
                    (timed.equilibrium, timed.viscous_dissipated_j_m2)
                }
                None => (old.point.trial(point.jump_m(), local_energy_tolerance_j_m2)
                    .map_err(|error| Error::Material { interface: id, error })?, 0.0),
            };
            if trial.point != point {
                return Err(Error::Invalid("mechanical endpoint lost history or contains an unresolved excursion"));
            }
            let closed = old.next_closed(point)?;
            let next = Interface { point, closed, ..old };
            let expected = old.thermal_link();
            let replacement = next.thermal_link();
            if old.closed && !closed { report.opened += 1; }
            if !old.closed && closed { report.reclosed += 1; }
            if old.closed && closed && expected != replacement { report.changed_conductances += 1; }
            let dissipated = old.reference.area_m2 * trial.energy.dissipated;
            let viscous = old.reference.area_m2 * viscous_density;
            if !nonnegative(dissipated) || !nonnegative(viscous) {
                return Err(Error::Invalid("unrepresentable local cohesive dissipation"));
            }
            record(id, old.reference, dissipated, viscous)?;
            report.fracture_dissipation_j += dissipated;
            report.viscous_dissipation_j += viscous;
            if !nonnegative(report.fracture_dissipation_j) || !nonnegative(report.viscous_dissipation_j) {
                return Err(Error::Invalid("unrepresentable aggregate cohesive dissipation"));
            }
            updates.push(ContactUpdate { id, expected, replacement });
            contacts.insert(id, next);
        }
        let thermal = self.thermal.replace_contacts(&updates, &mut checkpoint)?;
        checkpoint()?;
        Ok((Self { thermal, contacts: Contacts(contacts), last_synchronized_time_s: time_s }, report))
    }

    /// Executor-scoped endpoint synchronization with identical atomic semantics.
    pub fn synchronize_in(&self, cx: &Cx<'_>, time_s: f64,
        proposed: &[(u64, CohesivePoint<'m>)], local_energy_tolerance_j_m2: f64,
    ) -> Result<(Self, ContactReport), Error> {
        self.synchronize(time_s, proposed, local_energy_tolerance_j_m2,
            || cx.checkpoint().map_err(|_| DepositionError::Cancelled))
    }

    /// Executor-scoped fracture heating with identical joint commit semantics.
    pub fn synchronize_with_heat_in(&self, cx: &Cx<'_>, time_s: f64,
        proposed: &[(u64, CohesivePoint<'m>)], partitions: &[FractureHeatPartition],
        local_energy_tolerance_j_m2: f64, heat_energy_tolerance_j: f64,
    ) -> Result<(Self, FractureHeatReport), Error> {
        self.synchronize_with_heat(time_s, proposed, partitions,
            local_energy_tolerance_j_m2, heat_energy_tolerance_j,
            || cx.checkpoint().map_err(|_| DepositionError::Cancelled))
    }
}

// Deterministic physical-energy accumulation through the existing EFT owner.
// Final rounding tails are retained for the caller's heat-conversion budget.
#[derive(Default)]
struct HeatSum { high: f64, low: f64 }
impl HeatSum {
    fn add(&mut self, value: f64) -> Result<(), Error> {
        let (high, error) = fs_math::eft::two_sum(self.high, value);
        let low = self.low + error;
        if !high.is_finite() || !low.is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        self.high = high;
        self.low = low;
        Ok(())
    }
    fn rounded(&self) -> Result<(f64, f64), Error> {
        let result = fs_math::eft::two_sum(self.high, self.low);
        if !result.0.is_finite() || !result.1.is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        Ok(result)
    }
    fn value(&self) -> Result<f64, Error> { Ok(self.rounded()?.0) }
}
fn positive(v: f64) -> bool { v.is_finite() && v > 0.0 }
fn nonnegative(v: f64) -> bool { v.is_finite() && v >= 0.0 }

#[cfg(test)]
mod tests;
