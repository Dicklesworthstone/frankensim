//! Electrical contact response evaluated from the SAME material/phase epoch.
//!
//! Bulk resistivity, heat transport, and each ordered contact now consume one
//! set of material samples. This module supplies no contact constitutive law:
//! the caller owns surface condition, interface temperature, pressure, phase
//! dependence and heat partition. All current/heat numerics remain unchanged.

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

/// One unchanged contact location and its evaluated numerical/source response.
#[derive(Debug, PartialEq)]
pub struct ElectricalInterfaceSample<J> {
    /// Authored topology and reference values, never rewritten as new inputs.
    pub authored: ElectricalContact,
    /// Actual resistance/partition and their original material payload.
    pub response: ElectricalInterfaceResponse<J>,
}
impl<J> ElectricalInterfaceSample<J> {
    /// Descriptor consumed by the existing voxel current/contact solver.
    pub fn applied(&self) -> ElectricalContact {
        ElectricalContact { cells: self.authored.cells,
            resistance_ohm_m2: self.response.resistance_ohm_m2,
            heat_fraction_to_first: self.response.heat_fraction_to_first }
    }
}

/// Electrical interfaces in input order at one common physical evaluation time.
#[derive(Debug, PartialEq)]
pub struct ElectricalInterfaceEpoch<J> {
    /// Same evaluation time as the paired bulk/thermal-transport epoch [s].
    pub time_s: f64,
    /// Complete authored contact set, each evaluated exactly once this epoch.
    pub contacts: Vec<ElectricalInterfaceSample<J>>,
}

/// Original Joule/thermal accounts plus the electrical interfaces they used.
#[derive(Debug)]
pub struct JouleInterfaceReport<F, I, J> {
    /// Existing source, bulk/thermal material epochs and energy accounting.
    pub transport: JouleTransportReport<F, I>,
    /// Last electrical interfaces that supplied heat; None at zero duration.
    pub last_applied: Option<Arc<ElectricalInterfaceEpoch<J>>>,
    /// Fresh endpoint support, not a source applied during the preceding step.
    pub endpoint: Arc<ElectricalInterfaceEpoch<J>>,
}

type InterfaceLaw<'a, F, J> = dyn Fn(usize, &ElectricalContact,
    [&EnthalpyTransportSample<F>; 2]) -> Result<ElectricalInterfaceResponse<J>, JouleEnthalpyError> + 'a;

struct InterfaceSource<'a, F, I, J> {
    transport: TransportSource<'a, F, I>,
    interface: &'a InterfaceLaw<'a, F, J>,
    last: Option<Arc<ElectricalInterfaceEpoch<J>>>,
}
impl<F, I, J> Clone for InterfaceSource<'_, F, I, J> {
    fn clone(&self) -> Self {
        Self { transport: self.transport.clone(), interface: self.interface, last: self.last.clone() }
    }
}

// Establish bounded, safe ordered indexing before any constitutive callback.
// Conducting/electrode admissibility still belongs to the original solver.
fn admit_interfaces(input: JouleEnthalpyInput<'_>, gate: &CancelGate) -> Result<(), JouleEnthalpyError> {
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

impl<F, I, J> InterfaceSource<'_, F, I, J> {
    fn evaluate(&mut self, thermal: &DepositionSimulation<'_>)
        -> Result<(EnthalpyCurrentState, Arc<EnthalpyTransportEpoch<F, I>>,
            Arc<ElectricalInterfaceEpoch<J>>), JouleEnthalpyError> {
        let bulk = Arc::new(self.transport.prepare(thermal)?);
        let input = self.transport.account.input;
        let gate = self.transport.account.gate;
        let mut samples = Vec::with_capacity(input.contacts.len());
        let mut applied = Vec::with_capacity(input.contacts.len());
        for (index, &authored) in input.contacts.iter().enumerate() {
            checkpoint(gate)?;
            let sides = authored.cells.map(|c| &bulk.cells[c]);
            let response = (self.interface)(index, &authored, sides)?;
            let sample = ElectricalInterfaceSample { authored, response };
            applied.push(sample.applied());
            samples.push(sample);
        }
        checkpoint(gate)?;
        let resolve = |c: usize, _: &DepositedCell, _: f64| Ok(bulk.cells[c].response.resistivity_ohm_m());
        let current = self.transport.account.evaluate_with_contacts(thermal, &resolve, &applied)?;
        let interfaces = Arc::new(ElectricalInterfaceEpoch { time_s: bulk.time_s, contacts: samples });
        Ok((current, bulk, interfaces))
    }
}
impl<F, I, J> DepositionObserver for InterfaceSource<'_, F, I, J> {
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
    Ok((next, observer, JouleInterfaceReport { transport, last_applied: source.last, endpoint }))
}
