//! Electrical contact response evaluated from the SAME material/phase epoch.
//!
//! Bulk resistivity, heat transport, and each ordered contact now consume one
//! set of material samples. This module supplies no contact constitutive law:
//! the caller owns surface condition, interface temperature, pressure, phase
//! dependence and heat partition. All current/heat numerics remain unchanged.

/// Closed endpoint feedback through the shared thermal and nonlinear owners.
#[cfg(feature = "implicit-enthalpy-coupling")]
pub mod implicit;

use super::*;
use fs_conduction::observers::{CoupledObservers, ObserverCouplingError};

/// Source-owned contact evaluation, without requiring Clone on its evidence.
#[derive(Debug, PartialEq)]
pub struct ElectricalInterfaceResponse<J> {
    material: J,
    resistance_ohm_m2: f64,
    heat_fraction_to_first: f64,
}
impl<J> ElectricalInterfaceResponse<J> {
    /// Explicit area-specific resistance and heat partition in surface order.
    /// Zero resistance means perfect electrical contact, not missing data or an
    /// open circuit. This response cannot change either bulk's conductivity.
    pub fn try_new(material: J, resistance_ohm_m2: f64, heat_fraction_to_first: f64)
        -> Result<Self, JouleEnthalpyError> {
        if !resistance_ohm_m2.is_finite() || resistance_ohm_m2 < 0.0
            || !heat_fraction_to_first.is_finite()
            || !(0.0..=1.0).contains(&heat_fraction_to_first) {
            return Err(invalid("interface needs finite nonnegative resistance and heat fraction in [0,1]"));
        }
        Ok(Self { material, resistance_ohm_m2, heat_fraction_to_first })
    }
    /// Original constitutive payload and any retained usage receipts.
    pub const fn material(&self) -> &J { &self.material }
    /// Area-specific resistance, distinct from bulk resistivity [ohm m2].
    pub const fn resistance_ohm_m2(&self) -> f64 { self.resistance_ohm_m2 }
    /// Fraction of interface heat into the first authored surface.
    pub const fn heat_fraction_to_first(&self) -> f64 { self.heat_fraction_to_first }
}

impl ElectricalInterfaceResponse<crate::conjugate::material_coupling::interfaces::MaterialContact> {
    /// Bind an existing resolved interface to the ACTUAL adjacent surfaces.
    /// `MaterialContact` owns material/phase/process, texture-frame and ordered
    /// surface admission. Resistance is extracted from the same retained source;
    /// no separately copied coefficient can disagree with its usage receipt.
    ///
    /// The template fixes only the ordered cells. Heat partition is explicitly
    /// supplied, independently of resistance, and may evolve with the caller's
    /// interface model. Surface arguments must come from actual scenario/material
    /// samples, not copies of the source card's own declarations.
    pub fn try_from_source(
        template: &ElectricalContact,
        source: fs_material::electrical_interface::ResolvedElectricalInterface,
        surfaces: [&fs_matdb::SurfaceSpec; 2],
        heat_fraction_to_first: f64,
    ) -> Result<Self, JouleEnthalpyError> {
        let material = crate::conjugate::material_coupling::interfaces::MaterialContact::try_new(
            source, template.cells, surfaces, heat_fraction_to_first,
        )?;
        material.into_enthalpy_response(template.cells).map_err(Into::into)
    }
}

/// Explicit electrical connectivity supplied by the interface material owner.
/// Open is a supported constitutive result, never a fallback for missing data.
/// Both variants retain the original payload without requiring Clone.
#[derive(Debug, PartialEq)]
pub enum ElectricalInterfaceState<J> {
    /// A finite, possibly perfect, ohmic contact with explicit heat partition.
    Closed(ElectricalInterfaceResponse<J>),
    /// Electrically open; no leakage or contact heat is manufactured. The
    /// original material and thermal/mechanical interfaces remain present.
    Open(J),
}
impl<J> ElectricalInterfaceState<J> {
    /// Original constitutive/source payload, including while electrically open.
    pub const fn material(&self) -> &J {
        match self { Self::Closed(response) => response.material(), Self::Open(material) => material }
    }
    /// True only for a deliberately open electrical face.
    pub const fn is_open(&self) -> bool { matches!(self, Self::Open(_)) }
}
impl<J> From<ElectricalInterfaceResponse<J>> for ElectricalInterfaceState<J> {
    fn from(response: ElectricalInterfaceResponse<J>) -> Self { Self::Closed(response) }
}

// Only these checked responses can lower into the common electrical owner.
trait InterfaceEvaluation {
    fn applied_contact(&self, authored: ElectricalContact) -> Option<ElectricalContact>;
}
impl<J> InterfaceEvaluation for ElectricalInterfaceResponse<J> {
    fn applied_contact(&self, authored: ElectricalContact) -> Option<ElectricalContact> {
        Some(ElectricalContact { cells: authored.cells, resistance_ohm_m2: self.resistance_ohm_m2,
            heat_fraction_to_first: self.heat_fraction_to_first })
    }
}
impl<J> InterfaceEvaluation for ElectricalInterfaceState<J> {
    fn applied_contact(&self, authored: ElectricalContact) -> Option<ElectricalContact> {
        match self { Self::Closed(response) => response.applied_contact(authored), Self::Open(_) => None }
    }
}

/// One unchanged contact location and its evaluated numerical/source response.
#[derive(Debug, PartialEq)]
pub struct InterfaceSample<R> {
    /// Authored topology and reference values, never rewritten as new inputs.
    pub authored: ElectricalContact,
    /// Actual response and its original material payload.
    pub response: R,
}
/// Closed-contact observation; the existing public response semantics are unchanged.
pub type ElectricalInterfaceSample<J> = InterfaceSample<ElectricalInterfaceResponse<J>>;
/// Open/closed observation at one declared interface location.
pub type SwitchingInterfaceSample<J> = InterfaceSample<ElectricalInterfaceState<J>>;

impl<J> InterfaceSample<ElectricalInterfaceResponse<J>> {
    /// Descriptor consumed by the existing voxel current/contact solver.
    pub fn applied(&self) -> ElectricalContact {
        ElectricalContact { cells: self.authored.cells,
            resistance_ohm_m2: self.response.resistance_ohm_m2,
            heat_fraction_to_first: self.response.heat_fraction_to_first }
    }
}
impl<J> InterfaceSample<ElectricalInterfaceState<J>> {
    /// The applied finite contact, or None for an explicit opening. An opening
    /// has no fabricated numeric resistance or heat-partition coefficient.
    pub fn applied(&self) -> Option<ElectricalContact> {
        self.response.applied_contact(self.authored)
    }
}

/// Electrical interfaces in input order at one common physical evaluation time.
#[derive(Debug, PartialEq)]
pub struct InterfaceEpoch<R> {
    /// Same evaluation time as the paired bulk/thermal-transport epoch [s].
    pub time_s: f64,
    /// Complete authored contact set, each evaluated exactly once this epoch.
    pub contacts: Vec<InterfaceSample<R>>,
}
/// Closed-contact epoch used by the original material interface entry point.
pub type ElectricalInterfaceEpoch<J> = InterfaceEpoch<ElectricalInterfaceResponse<J>>;
/// Open/closed epoch from the same shared physical evaluator.
pub type SwitchingInterfaceEpoch<J> = InterfaceEpoch<ElectricalInterfaceState<J>>;

/// Original Joule/thermal accounts plus the electrical interfaces they used.
#[derive(Debug)]
pub struct InterfaceReport<F, I, R> {
    /// Existing source, bulk/thermal material epochs and energy accounting.
    pub transport: JouleTransportReport<F, I>,
    /// Last electrical interfaces that supplied heat; None at zero duration.
    pub last_applied: Option<Arc<InterfaceEpoch<R>>>,
    /// Fresh endpoint support, not a source applied during the preceding step.
    pub endpoint: Arc<InterfaceEpoch<R>>,
}
/// Original finite-contact report and accessors.
pub type JouleInterfaceReport<F, I, J> = InterfaceReport<F, I, ElectricalInterfaceResponse<J>>;
/// Material-dependent open/closed contacts with the same current/heat accounts.
pub type JouleSwitchingReport<F, I, J> = InterfaceReport<F, I, ElectricalInterfaceState<J>>;

type InterfaceLaw<'a, F, R> = dyn Fn(usize, &ElectricalContact,
    [&EnthalpyTransportSample<F>; 2]) -> Result<R, JouleEnthalpyError> + 'a;

struct InterfaceSource<'a, F, I, R> {
    transport: TransportSource<'a, F, I>,
    interface: &'a InterfaceLaw<'a, F, R>,
    last: Option<Arc<InterfaceEpoch<R>>>,
}
impl<F, I, R> Clone for InterfaceSource<'_, F, I, R> {
    fn clone(&self) -> Self {
        Self { transport: self.transport.clone(), interface: self.interface, last: self.last.clone() }
    }
}

// Establish bounded, safe ordered indexing before any constitutive callback.
// Conducting/electrode admissibility still belongs to the original solver.
pub(in crate::conjugate::enthalpy) fn admit_interfaces(input: JouleEnthalpyInput<'_>, gate: &CancelGate) -> Result<(), JouleEnthalpyError> {
    checkpoint(gate)?;
    let n = input.domain.cell_count();
    if input.contacts.len() > n.saturating_mul(3) {
        return Err(invalid("too many electrical interfaces for the voxel domain"));
    }
    let mut faces = BTreeSet::new();
    for contact in input.contacts {
        checkpoint(gate)?;
        let [a, b] = contact.cells;
        if a >= n || b >= n || !(0..6).any(|f| input.domain.neighbor(a, f) == Some(b))
            || !faces.insert((a.min(b), a.max(b))) {
            return Err(invalid("electrical interfaces need unique in-domain face-adjacent cell pairs"));
        }
    }
    Ok(())
}

impl<F, I, R: InterfaceEvaluation> InterfaceSource<'_, F, I, R> {
    fn evaluate(&mut self, thermal: &DepositionSimulation<'_>)
        -> Result<(EnthalpyCurrentState, Arc<EnthalpyTransportEpoch<F, I>>,
            Arc<InterfaceEpoch<R>>), JouleEnthalpyError> {
        let bulk = Arc::new(self.transport.prepare(thermal)?);
        let input = self.transport.account.input;
        let gate = self.transport.account.gate;
        let mut samples = Vec::with_capacity(input.contacts.len());
        let mut applied = Vec::with_capacity(input.contacts.len());
        let mut openings = Vec::new();
        for (index, &authored) in input.contacts.iter().enumerate() {
            checkpoint(gate)?;
            let sides = authored.cells.map(|c| &bulk.cells[c]);
            let response = (self.interface)(index, &authored, sides)?;
            match response.applied_contact(authored) {
                Some(contact) => applied.push(contact),
                None => openings.push(authored.cells),
            }
            samples.push(InterfaceSample { authored, response });
        }
        checkpoint(gate)?;
        let resolve = |c: usize, _: &DepositedCell, _: f64| Ok(bulk.cells[c].response.resistivity_ohm_m());
        let current = self.transport.account.evaluate_with_contact_topology(thermal, &resolve, &applied, &openings)?;
        let interfaces = Arc::new(InterfaceEpoch { time_s: bulk.time_s, contacts: samples });
        Ok((current, bulk, interfaces))
    }
}
impl<F, I, R: InterfaceEvaluation> DepositionObserver for InterfaceSource<'_, F, I, R> {
    type Error = JouleEnthalpyError;
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.transport.event(event, thermal)
    }
    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        let (current, bulk, interfaces) = self.evaluate(thermal)?;
        let power = self.transport.account.apply(current)?;
        self.transport.last = Some(bulk);
        self.last = Some(interfaces);
        Ok(power)
    }
    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.transport.conductances(thermal, check)
    }
    fn step(&mut self, start: f64, dt: f64, old: &BTreeMap<CellId, f64>, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> {
        self.transport.step(start, dt, old, thermal, check)
    }
}

/// Advance material-dependent electrical CONTACTS, bulk current and thermal
/// transport together, optionally composing other physical histories (`&()`
/// means no additional owner). This uses the existing enthalpy time integrator.
///
/// `input.contacts` supplies the fixed ordered contact topology. The interface
/// callback replaces its resistance and heat partition, once per contact per
/// epoch, and receives the SAME ordered bulk samples used by thermal links and
/// current. Authored coefficients remain available only as explicit reference
/// data; they are not added to the resolved resistance. Unlisted faces retain
/// the original perfect-contact semantics. Contacts to insulating bulk refuse.
/// No pressure, interface temperature, surface law or partition is inferred.
///
/// One actual current solve supplies all bulk/contact watts exactly once to the
/// shared thermal step. Whole-advance current budgets also bound sampling epochs.
/// The final material/link/contact/current evaluation checks endpoint support;
/// it produces no extra heat. Any late refusal or cancellation returns NO new
/// thermal or observer history. Concrete payloads need not implement Clone.
///
/// Electrical and thermal contacts are distinct physical laws and selections.
/// Other observers may supply signed heat or unowned thermal-link conductances.
/// A duplicate thermal conductance owner refuses rather than silently winning.
/// Callbacks must be pure with respect to physical history. This remains an
/// explicit, fixed-volume, stationary Ohmic model, not contact mechanics, phase
/// kinetics, induction, strong coupling or a source-stiffness stability claim.
/// Discontinuities require advance endpoints; refine the declared maximum dt.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_joule_enthalpy_with_material_interfaces_observed<'t, F, I, J, O: DepositionObserver>(
    thermal: &DepositionSimulation<'t>, target_s: f64, config: StepConfig,
    input: JouleEnthalpyInput<'_>, link_ids: &[u64], gate: &CancelGate, observer: &O,
    resolve: impl Fn(usize, &DepositedCell, f64) -> Result<EnthalpyTransport<F>, JouleEnthalpyError>,
    resolve_link: impl Fn(&ContactPatch, [&EnthalpyTransportSample<F>; 2]) -> Result<(f64, I), JouleEnthalpyError>,
    resolve_interface: impl Fn(usize, &ElectricalContact, [&EnthalpyTransportSample<F>; 2])
        -> Result<ElectricalInterfaceResponse<J>, JouleEnthalpyError>,
) -> Result<(DepositionSimulation<'t>, O, JouleInterfaceReport<F, I, J>),
    ObserverCouplingError<JouleEnthalpyError, O::Error>> {
    advance_interfaces(thermal, target_s, config, input, link_ids, gate, observer,
        resolve, resolve_link, resolve_interface)
}

/// Evolve electrically opening/reclosing interfaces on the SAME material and
/// enthalpy timeline as bulk current, contact heating and thermal transport.
///
/// The callback returns an explicit supported Open or Closed state per authored
/// location, using the same adjacent material samples as the closed-contact API.
/// Open faces participate in electrical reachability, matrix assembly and heat
/// deposition consistently; they carry zero current and zero contact heat. Bulk
/// properties are still evaluated. Closed zero resistance is never an opening.
///
/// Only electrical connectivity changes. Mass, enthalpy, geometry and thermal
/// contact coverage survive opening/reclosing. Independently owned heat sources
/// and conductance laws continue to evolve the body. The payload for every
/// interface, including an opening, is retained in applied/endpoint epochs.
/// Missing data must be an error, never Open. Neither a phase name nor a
/// temperature threshold is hard-coded in this driver.
///
/// Connectivity is sampled at actual explicit substep starts and at the final
/// endpoint. State-triggered switches are therefore resolved to that time grid,
/// NOT localized by an event root. Refine max_dt_s to resolve switching times;
/// prescribed discontinuities require advance boundaries or scheduled events.
/// The endpoint check cannot retroactively change heat delivered by a prior
/// step. This API does not supply contact mechanics, fusion/arc kinetics,
/// electromagnetic storage, charge transients or a general stability theorem.
///
/// Disconnected nonzero ideal current sources refuse. Zero-current floating
/// terminals retain undefined voltage and exactly zero work. Any material,
/// endpoint, budget or cancellation failure discards all candidate histories.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_joule_enthalpy_with_switching_interfaces_observed<'t, F, I, J, O: DepositionObserver>(
    thermal: &DepositionSimulation<'t>, target_s: f64, config: StepConfig,
    input: JouleEnthalpyInput<'_>, link_ids: &[u64], gate: &CancelGate, observer: &O,
    resolve: impl Fn(usize, &DepositedCell, f64) -> Result<EnthalpyTransport<F>, JouleEnthalpyError>,
    resolve_link: impl Fn(&ContactPatch, [&EnthalpyTransportSample<F>; 2]) -> Result<(f64, I), JouleEnthalpyError>,
    resolve_interface: impl Fn(usize, &ElectricalContact, [&EnthalpyTransportSample<F>; 2])
        -> Result<ElectricalInterfaceState<J>, JouleEnthalpyError>,
) -> Result<(DepositionSimulation<'t>, O, JouleSwitchingReport<F, I, J>),
    ObserverCouplingError<JouleEnthalpyError, O::Error>> {
    advance_interfaces(thermal, target_s, config, input, link_ids, gate, observer,
        resolve, resolve_link, resolve_interface)
}

// All finite-contact and switching entry points share this transaction and
// exactly one material evaluator, electrical solve and thermal time integrator.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance_interfaces<'t, F, I, R: InterfaceEvaluation, O: DepositionObserver>(
    thermal: &DepositionSimulation<'t>, target_s: f64, config: StepConfig,
    input: JouleEnthalpyInput<'_>, link_ids: &[u64], gate: &CancelGate, observer: &O,
    resolve: impl Fn(usize, &DepositedCell, f64) -> Result<EnthalpyTransport<F>, JouleEnthalpyError>,
    resolve_link: impl Fn(&ContactPatch, [&EnthalpyTransportSample<F>; 2]) -> Result<(f64, I), JouleEnthalpyError>,
    resolve_interface: impl Fn(usize, &ElectricalContact, [&EnthalpyTransportSample<F>; 2])
        -> Result<R, JouleEnthalpyError>,
) -> Result<(DepositionSimulation<'t>, O, InterfaceReport<F, I, R>),
    ObserverCouplingError<JouleEnthalpyError, O::Error>> {
    admit_map(thermal, input, gate).map_err(ObserverCouplingError::First)?;
    admit_interfaces(input, gate).map_err(ObserverCouplingError::First)?;
    let transport = TransportSource { account: Account::new(input, gate), selected: link_ids,
        indices: input.cell_ids.iter().copied().enumerate().map(|(c, id)| (id, c)).collect(),
        bulk: &resolve, link: &resolve_link, last: None };
    let source = InterfaceSource { transport, interface: &resolve_interface, last: None };
    let owners = CoupledObservers::new(source, observer.clone());
    let (next, owners, audit) = thermal.advance_observed(&owners, target_s, config, || checkpoint(gate))?;
    let (mut source, observer) = owners.into_parts();
    let (current, bulk, endpoint) = source.evaluate(&next).map_err(ObserverCouplingError::First)?;
    let joule = source.transport.account.finish(audit, current).map_err(ObserverCouplingError::First)?;
    let transport = JouleTransportReport { joule, last_applied: source.transport.last, endpoint: bulk };
    checkpoint(gate)?;
    Ok((next, observer, InterfaceReport { transport, last_applied: source.last, endpoint }))
}
