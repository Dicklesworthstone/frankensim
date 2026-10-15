//! Several closed gas chambers and many heat contacts on ONE finite network.
//! The circuit/enthalpy owner supplies the outer residual. Each gas is advanced
//! once per trial with all its contacts, never sequentially once per contact.

/// Data-only material contacts with aggregate surface-area admission.
pub mod material;

use super::*;
use crate::circuit_enthalpy::{self, CircuitEnthalpyCandidate, CircuitEnthalpyError,
    CircuitEnthalpyInput, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients};
use crate::material_circuit::MaterialCircuit;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::cell::Cell;

/// One independent accepted gas and translating wall, with its original sources.
#[derive(Clone, Copy, Debug)]
pub struct Chamber<'a> {
    /// Stable caller-owned physical chamber ID.
    pub id: u64,
    /// Immutable accepted gas, volume and momentum through every trial.
    pub piston: &'a GasPiston,
    /// Prescribed gas heat and mechanical loads, on the network's same interval.
    pub input: PistonInput,
    /// Original momentum and gas-caloric controls.
    pub piston_solve: PistonSolve,
    /// Signed TOTAL heat bracket for this chamber's contacts, not for each one.
    pub exchange: HeatExchangeSolve,
}
/// One distinct heat contact; several may share a chamber or a thermal cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Contact {
    /// Stable contact ID, unique across the entire coupling.
    pub id: u64,
    /// Existing physical chamber ID, not a vector index.
    pub chamber: u64,
    /// Existing active thermal cell to debit by this contact's signed heat.
    pub cell: CellId,
}
/// Global bounds, in addition to each existing solver's bounded work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkSolve {
    /// Existing circuit/enthalpy Newton and physical residual controls.
    pub network: CircuitEnthalpySolve,
    /// Maximum independent gas states, checked before allocation.
    pub maximum_chambers: usize,
    /// Maximum total contacts, checked before allocation and material queries.
    pub maximum_contacts: usize,
    /// All chamber solves, including rejected outer trials and final admission.
    pub maximum_boundary_solves: usize,
    /// All-owner absolute energy residual and heat-transfer allowance [J].
    pub total_energy_tolerance_j: f64,
}
/// Original source errors retain their domain and chamber identity.
#[derive(Debug)]
pub enum NetworkError<E, F> {
    /// Circuit/thermal material callback refusal.
    Material(E),
    /// Original gas/momentum/contact refusal, with its physical chamber.
    Boundary { chamber: u64, error: ThermalPistonError<F> },
    /// Thermal owner refused a contact temperature.
    Thermal(DepositionError),
    /// Global boundary work limit exhausted.
    Budget,
    /// Contact heat overflowed, vanished or failed its receiving-cell balance.
    Transfer,
}
impl<E: core::fmt::Debug, F: core::fmt::Debug> core::fmt::Display for NetworkError<E,F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f,"multi-chamber coupling: {self:?}") }
}
impl<E: std::error::Error+'static, F: std::error::Error+'static> std::error::Error for NetworkError<E,F> {}

/// Original chamber endpoint and per-contact sources, in stable ID order.
pub struct ChamberBoundary<I> {
    /// Original physical chamber identity.
    pub id: u64,
    /// Exact shared-midpoint gas proposal. Each source also retains its cell map.
    pub boundary: ThermalPortsCandidate<(Contact,I)>,
}
/// One matching material evaluation and all gas candidates for the same interval.
pub struct NetworkSources<S,I> {
    /// Original material-network source payload.
    pub material: S,
    /// All gas candidates in stable chamber-ID order.
    pub chambers: Vec<ChamberBoundary<I>>,
    /// Actual signed heat debit handed to the thermal source map [J].
    pub transferred_heat_j: f64,
}
/// Combined physical balance and work, not a time-accuracy certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkStep {
    /// Network boundary/source heat excluding internal electrical/contact heat [J].
    pub external_network_heat_j: f64,
    /// Signed evaluated heat removed from all network cells [J].
    pub transferred_heat_j: f64,
    /// Circuit + cells + all gas/kinetic storage and exterior work minus supply [J].
    pub total_energy_residual_j: f64,
    /// Every chamber solve, including discarded numerical trials.
    pub boundary_solves: usize,
    /// Every nested scalar heat proposal across all chamber solves.
    pub exchange_trials: usize,
    /// Every nested gas-material proposal across all chamber solves.
    pub material_trials: usize,
}
/// No accepted circuit, thermal cell or gas is changed by constructing this.
pub struct NetworkCandidate<'m,S,I> {
    network: CircuitEnthalpyCandidate<'m,NetworkSources<S,I>>,
    step: NetworkStep,
}
impl<'m,S,I> NetworkCandidate<'m,S,I> {
    /// Proposed thermal history, with original masses, laws and event cursor.
    pub fn thermal(&self) -> &DepositionSimulation<'m> { self.network.thermal() }
    /// Proposed conserved electrical charge/flux history.
    pub fn circuit(&self) -> &MaterialCircuit { self.network.circuit() }
    /// Proposed gases and all original contact evaluations in chamber-ID order.
    pub fn chambers(&self) -> &[ChamberBoundary<I>] { &self.network.report().source.chambers }
    /// Original material-network evaluation at the accepted thermal endpoint.
    pub fn material(&self) -> &S { &self.network.report().source.material }
    /// Actual independent all-owner energy check and bounded-work counts.
    pub const fn step(&self) -> &NetworkStep { &self.step }
    /// Publish the network and all chamber candidates together. Contact and
    /// electrical heat are ALREADY included. Gas endpoints are retained in the
    /// report source alongside their original interface records; move them with
    /// `boundary.into_parts()` rather than advancing or depositing heat again.
    #[allow(clippy::type_complexity)]
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialCircuit,
        CircuitEnthalpyReport<NetworkSources<S,I>>, NetworkStep) {
        let (thermal,circuit,report) = self.network.into_parts();
        (thermal,circuit,report,self.step)
    }
}

/// Solve all chambers, contacts, resistor bodies and finite electrical storage
/// over one interval. Each gas has ONE summed heat/momentum solve; each contact
/// has its OWN signed receiving-cell heat, including equal/opposite transfers.
///
/// Circuit and thermal geometry, C/L, gas masses and compositions are fixed.
/// Chambers have independent walls, not a shared crank or mechanical linkage.
/// Thermal endpoint coefficients with gas midpoint motion give first-order
/// accuracy in general. Source/event discontinuities require interval boundaries.
/// The caller owns surface coverage for custom contact laws; the material adapter
/// checks the sum of all patches. Callbacks must not mutate live physical history.
/// Repeated chamber IDs, aliased gas references, bad clocks/maps and budget,
/// material, convergence or cancellation failures return no partial endpoint.
#[allow(clippy::too_many_arguments,clippy::type_complexity)]
pub fn step_candidate<'m,S,I,E,F>(
    thermal: &DepositionSimulation<'m>, target_s: f64, electrical: CircuitEnthalpyInput<'_>,
    chambers: &[Chamber<'_>], contacts: &[Contact], controls: NetworkSolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>,E>,
    resolve_port: impl Fn(Contact,&DepositionSimulation<'m>,&S,&ResolvedThermochemicalGasState)
        -> Result<ConductanceSample<I>,F>,
    cancelled: impl FnMut() -> bool,
) -> Result<NetworkCandidate<'m,S,I>,CircuitEnthalpyError<NetworkError<E,F>>> {
    let cancel = RefCell::new(cancelled);
    let check = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
    check()?;
    let dt = target_s-thermal.time_s();
    if !positive(dt) || chambers.len()>controls.maximum_chambers || contacts.len()>controls.maximum_contacts
        || (!chambers.is_empty() && controls.maximum_boundary_solves==0)
        || !positive(controls.total_energy_tolerance_j) {
        return Err(CircuitEnthalpyError::Invalid("invalid shared clock, chamber/contact bounds or energy tolerance"));
    }
    let mut ordered = BTreeMap::new();
    for (index,chamber) in chambers.iter().enumerate() {
        check()?;
        if chamber.input.duration_s!=dt || ordered.insert(chamber.id,chamber).is_some()
            || chambers[..index].iter().any(|old| std::ptr::eq(old.piston,chamber.piston)) {
            return Err(CircuitEnthalpyError::Invalid("each chamber needs one unique gas owner, ID and shared clock"));
        }
    }
    let mut mapped = BTreeMap::new();
    for contact in contacts {
        check()?;
        if !ordered.contains_key(&contact.chamber) || !thermal.cells().contains_key(&contact.cell)
            || mapped.insert(contact.id,*contact).is_some() {
            return Err(CircuitEnthalpyError::Invalid("duplicate contact or unknown chamber/thermal cell"));
        }
    }
    let mut groups: BTreeMap<u64,Vec<Contact>> = ordered.keys().map(|&id|(id,Vec::new())).collect();
    for &contact in mapped.values() { groups.get_mut(&contact.chamber).expect("admitted chamber").push(contact); }
    let calls=Cell::new(0usize); let heat_trials=Cell::new(0usize); let material_trials=Cell::new(0usize);
    let candidate = circuit_enthalpy::step_candidate(thermal,target_s,electrical,controls.network,guess_j,|state| {
        let mut coefficients=resolve(state).map_err(NetworkError::Material)?;
        let mut boundaries=Vec::with_capacity(ordered.len());
        let mut cell_heat=BTreeMap::<CellId,f64>::new();
        let mut evaluated_heat=0.0;
        for (&id,chamber) in &ordered {
            check().map_err(NetworkError::Thermal)?;
            if calls.get()>=controls.maximum_boundary_solves { return Err(NetworkError::Budget); }
            let ports: Vec<_> = groups[&id].iter().map(|contact| state.temperature_k(contact.cell)
                .map(|temperature_k| ThermalPort {id:contact.id,temperature_k}))
                .collect::<Result<_,_>>().map_err(NetworkError::Thermal)?;
            calls.set(calls.get()+1);
            let boundary=chamber.piston.thermal_ports_candidate(chamber.input,chamber.piston_solve,
                chamber.exchange,&ports,controls.maximum_contacts,|port,gas| {
                    let contact=mapped[&port.id];
                    resolve_port(contact,state,&coefficients.source,gas).map(|sample| ConductanceSample {
                        conductance_w_k:sample.conductance_w_k,source:(contact,sample.source) })
                },||(*cancel.borrow_mut())())
                .map_err(|error|NetworkError::Boundary {chamber:id,error})?;
            heat_trials.set(heat_trials.get().checked_add(boundary.step().exchange_trials).ok_or(NetworkError::Budget)?);
            material_trials.set(material_trials.get().checked_add(boundary.step().material_trials).ok_or(NetworkError::Budget)?);
            for sample in boundary.source() {
                let sum=cell_heat.entry(sample.source.0.cell).or_insert(0.0);
                *sum+=sample.heat_j; evaluated_heat+=sample.heat_j;
                if !sum.is_finite() || !evaluated_heat.is_finite() { return Err(NetworkError::Transfer); }
            }
            boundaries.push(ChamberBoundary {id,boundary});
        }
        // Aggregate by the true receiving cell, not by gas or electrical branch.
        // Keep individual signed contacts above even when a cell's net is zero.
        let mut deposited=0.0;
        for (cell,q) in cell_heat {
            check().map_err(NetworkError::Thermal)?;
            let watts=q/dt;
            let source=coefficients.sources_w.entry(cell).or_insert(0.0);
            let previous=*source; let next=previous-watts;
            if !watts.is_finite() || !next.is_finite() || (q!=0.0 && (watts==0.0 || next==previous)) {
                return Err(NetworkError::Transfer);
            }
            *source=next; deposited+=dt*watts;
        }
        if !deposited.is_finite() || (deposited-evaluated_heat).abs()>controls.total_energy_tolerance_j {
            return Err(NetworkError::Transfer);
        }
        Ok(CircuitThermalCoefficients {resistors:coefficients.resistors,sources_w:coefficients.sources_w,
            conductances_w_k:coefficients.conductances_w_k,source:NetworkSources {material:coefficients.source,chambers:boundaries,transferred_heat_j:deposited}})
    },||(*cancel.borrow_mut())());
    let candidate=match candidate {
        Err(CircuitEnthalpyError::Source(NetworkError::Boundary {error:ThermalPistonError::Cancelled,..})) =>
            return Err(CircuitEnthalpyError::Cancelled),
        Err(CircuitEnthalpyError::Source(NetworkError::Thermal(DepositionError::Cancelled))) =>
            return Err(CircuitEnthalpyError::Cancelled),
        other=>other?,
    };
    let report=candidate.report();
    let external=report.non_electrical_heat_j+report.source.transferred_heat_j;
    let mut total=(report.electrical.stored_after_j-report.electrical.stored_before_j)
        +report.thermal.stored_energy_change_j-report.electrical.native.supplied-external;
    for chamber in &report.source.chambers {
        let p=&chamber.boundary.step().piston;
        total+=p.gas.stored_energy_change_j+p.kinetic_energy_change_j+p.exterior_pressure_work_j
            -p.external_force_work_j-ordered[&chamber.id].input.heat_j;
    }
    if !external.is_finite() || !total.is_finite() || total.abs()>controls.total_energy_tolerance_j {
        return Err(CircuitEnthalpyError::Invalid("all chambers, thermal bodies and circuit fail combined energy balance"));
    }
    check()?;
    let step=NetworkStep {external_network_heat_j:external,transferred_heat_j:report.source.transferred_heat_j,
        total_energy_residual_j:total,boundary_solves:calls.get(),exchange_trials:heat_trials.get(),
        material_trials:material_trials.get()};
    Ok(NetworkCandidate {network:candidate,step})
}

#[cfg(test)]
mod tests;
