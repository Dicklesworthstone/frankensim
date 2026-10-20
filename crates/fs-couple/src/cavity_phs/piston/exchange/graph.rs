//! Original-state multiport gas exchange coupled to independent inertial walls.
//!
//! All restrictions see one accepted gas snapshot. The existing conservative
//! network resolves every inlet/outlet together; each existing piston then
//! advances its post-transfer gas, heat and momentum over the SAME interval.
//! This transfer-then-motion split is first order, not a simultaneous implicit
//! flow solve. No second EOS, mixing law, pressure-work law or wall integrator.
//!
//! The enclosing balance uses actual final gas and kinetic storage. Nozzle
//! enthalpy cancels internally and is never added again as prescribed heat.
//! Independent walls have independent momentum; a shared crank or opposing
//! faces require their dedicated mechanics owner, not duplicated pistons.

use std::collections::BTreeMap;
use super::{ClosedGasParcel, GasEnergyError, GasPiston, PistonError, PistonInput,
    PistonSolve, PistonStep, positive, stored_change};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::{
    GasSnapshotSolve, MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::{
    self, SnapshotNozzleError, SnapshotNozzleLink, SnapshotNozzleStep,
};

/// Original piston input and admission. Its duration must equal the graph clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingGasDrive {
    /// External heat/force/pressure only; transferred enthalpy is already stored.
    pub input: PistonInput,
    /// Original volume, gas-source, impulse and work bounds.
    pub solve: PistonSolve,
}
/// Global acceptance supplements every original material and mechanical gate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingGasSolve {
    /// Shared-state mass/species/energy and gross original-inventory bounds.
    pub exchange: GasSnapshotSolve,
    /// Absolute all-gas/all-wall energy tolerance [J].
    pub total_energy_tolerance_j: f64,
}
/// Failure returns no graph endpoint, including already solved trial walls.
#[derive(Clone, Debug, PartialEq)]
pub enum MovingGasError {
    /// Missing/overlapping identities, mismatched clocks or invalid arithmetic.
    Invalid(&'static str),
    /// Original common-snapshot flow and inventory owner refused the proposal.
    Flow(SnapshotNozzleError),
    /// Original material/momentum owner refused this moving chamber.
    Piston { id: u64, error: PistonError },
    /// Cancellation before publication of all physical owners.
    Cancelled,
}
impl core::fmt::Display for MovingGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "moving material gas network: {self:?}")
    }
}
impl std::error::Error for MovingGasError {}

/// Stage accounts and an independent enclosing physical balance.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingGasStep {
    /// Original nozzle/throat, species and extensive-transfer accounts.
    pub flow: SnapshotNozzleStep,
    /// Original heat/work/momentum results in moving-chamber id order.
    pub pistons: BTreeMap<u64, PistonStep>,
    /// Actual total gas internal-energy change, including changed inventories [J].
    pub gas_energy_change_j: f64,
    /// Actual sum of the walls' kinetic-energy changes [J].
    pub kinetic_energy_change_j: f64,
    /// Work delivered against the independently declared exterior pressures [J].
    pub exterior_pressure_work_j: f64,
    /// Additional applied forces times actual wall displacements [J].
    pub external_force_work_j: f64,
    /// External gas heating only, excluding internal transported enthalpy [J].
    pub supplied_heat_j: f64,
    /// Fresh net mass discrepancy over the final physical endpoints [kg].
    pub total_mass_residual_kg: f64,
    /// dU_gases+dK_walls+W_exterior-W_applied-Q_external [J].
    pub total_energy_residual_j: f64,
}
/// Complete proposed physical owners. No separate, potentially stale gas cache.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingGasCandidate<P = GasPiston> {
    pistons: BTreeMap<u64, P>,
    reservoirs: BTreeMap<u64, ClosedGasParcel>,
    step: MovingGasStep,
}
impl<P> MovingGasCandidate<P> {
    /// Complete moving chambers, each owning its gas, geometry and momentum.
    #[must_use]
    pub fn pistons(&self) -> &BTreeMap<u64, P> { &self.pistons }
    /// Finite fixed-volume chambers; they are not infinite pressure sources.
    #[must_use]
    pub fn reservoirs(&self) -> &BTreeMap<u64, ClosedGasParcel> { &self.reservoirs }
    /// All original stage results and the independently checked whole balance.
    #[must_use]
    pub const fn step(&self) -> &MovingGasStep { &self.step }
    /// Publish all gas inventories and momenta together, without redepositing heat.
    #[must_use]
    pub fn into_parts(self) -> (BTreeMap<u64, P>, BTreeMap<u64, ClosedGasParcel>, MovingGasStep) {
        (self.pistons, self.reservoirs, self.step)
    }
    /// Preserve an existing source-binding wrapper around each numerical owner.
    pub(crate) fn map_pistons<Q>(self, mut wrap: impl FnMut(u64, P) -> Q) -> MovingGasCandidate<Q> {
        MovingGasCandidate { pistons: self.pistons.into_iter().map(|(id, piston)|
            (id, wrap(id, piston))).collect(), reservoirs: self.reservoirs, step: self.step }
    }
}
use MovingGasError as Error;
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn piston_error(id: u64, error: PistonError) -> Error {
    if error == PistonError::Cancelled { Error::Cancelled } else { Error::Piston { id, error } }
}
// Fixed-order sums avoid subtracting large formation-energy totals. They are
// only numerical accumulation, not another source of physical storage.
#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), Error> {
        let adjusted = term - self.correction;
        let next = self.value + adjusted;
        let correction = (next - self.value) - adjusted;
        if ![term, adjusted, next, correction].iter().all(|x| x.is_finite()) {
            return Err(Error::Invalid("moving gas network extensive sum is unrepresentable"));
        }
        self.value = next; self.correction = correction;
        Ok(())
    }
}

// The existing wall/heat port schema has no stream port. Keep its real transfers
// visible but never certify the open chamber from its closed two-port audit.
fn record_open_account(
    old: &GasPiston, next: &mut GasPiston, drive: MovingGasDrive,
    report: &PistonStep, gas_change: f64, tolerance_j: f64,
) -> Result<(), PistonError> {
    let initial_s = old.gas.state().specific_entropy_j_kg_k().map_err(GasEnergyError::from)?;
    let final_s = next.gas.state().specific_entropy_j_kg_k().map_err(GasEnergyError::from)?;
    let entropy = next.gas.mass_kg() * (final_s - initial_s)
        + (next.gas.mass_kg() - old.gas.mass_kg()) * initial_s;
    let storage = gas_change + report.kinetic_energy_change_j;
    if !entropy.is_finite() || !storage.is_finite() {
        return Err(PistonError::Invalid("open moving-chamber storage or entropy is unrepresentable"));
    }
    next.ports.accept(crate::port_adoption::StepAccount {
        tick: old.ports.next_tick(), duration_s: drive.input.duration_s,
        storage_change_j: Some(storage), loss_j: Some(0.0), included_work_j: 0.0,
        entropy_change_j_k: entropy, tolerance_j: Some(tolerance_j),
        unavailable: Some("multiport gas exchange needs stream ports; use the enclosing moving-network balance"),
    }, [crate::port_adoption::PortTransfer::work(
        report.external_force_work_j - report.exterior_pressure_work_j, true),
        crate::port_adoption::PortTransfer::heat(drive.input.heat_j, report.midpoint.state().temperature)]);
    Ok(())
}

/// Advance arbitrary branching/cyclic gas plumbing and independent moving walls.
///
/// IDs must be disjoint across the two owner maps. `endpoints` covers their union
/// exactly; `drives` covers the moving map exactly, with the SAME `duration_s`.
/// Every original flow sees the accepted snapshot, not another wall's updated
/// pressure. The flow-stage withdrawal/equilibrium checks remain unchanged.
/// All wall trials then start from the fully reconstructed post-transfer state.
///
/// This first-order operator split advances time once. Motion can reverse a
/// pressure ordering after the flow stage; the next interval selects flow again.
/// No flow event locator, jet thrust, implicit-network stability or exact entropy
/// conservation is claimed. Split at prescribed events and refine time.
/// Fixed chambers receive no external heating in this call. Walls are independent
/// constant-area translations, not opposite faces of one shared inertial body.
///
/// Flow, all walls, port history and full balances must succeed before return.
/// Zero-flow paths retain ordinary closed-piston behavior and its port account.
/// Any nonzero gross inlet OR outlet marks that chamber's two-port-only audit
/// unavailable, even when net transported mass/enthalpy happens to be zero.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    pistons: &BTreeMap<u64, GasPiston>, reservoirs: &BTreeMap<u64, ClosedGasParcel>,
    links: &[SnapshotNozzleLink], duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>, drives: &BTreeMap<u64, MovingGasDrive>,
    controls: MovingGasSolve, mut cancelled: impl FnMut() -> bool,
) -> Result<MovingGasCandidate, Error> {
    check(&mut cancelled)?;
    let count = pistons.len().checked_add(reservoirs.len())
        .ok_or(Error::Invalid("gas network size overflow"))?;
    if count == 0 || count > MAX_SNAPSHOT_PARCELS || links.len() > MAX_SNAPSHOT_TRANSFERS
        || !positive(duration_s) || !positive(controls.total_energy_tolerance_j)
        || !pistons.keys().eq(drives.keys()) || endpoints.len() != count
        || pistons.keys().any(|id| reservoirs.contains_key(id)) {
        return Err(Error::Invalid("bounded disjoint gas owners, exact drive/endpoint maps and positive clock/tolerance required"));
    }
    for drive in drives.values() {
        check(&mut cancelled)?;
        if drive.input.duration_s != duration_s {
            return Err(Error::Invalid("all moving chambers must use the same physical interval"));
        }
    }
    let mut initial = BTreeMap::new();
    for (&id, gas) in reservoirs {
        check(&mut cancelled)?;
        initial.insert(id, gas.clone());
    }
    for (&id, piston) in pistons {
        check(&mut cancelled)?;
        initial.insert(id, piston.gas.clone());
    }
    let candidate = snapshot::step_candidate(&initial, links, duration_s, endpoints,
        controls.exchange, &mut cancelled).map_err(|error| {
            if error == SnapshotNozzleError::Cancelled { Error::Cancelled } else { Error::Flow(error) }
        })?;
    let (mut remaining, flow) = candidate.into_parts();
    let open: std::collections::BTreeSet<_> = flow.exchange.parcels.iter()
        .filter(|p| p.incoming_mass_kg > 0.0 || p.outgoing_mass_kg > 0.0).map(|p| p.id).collect();
    let mut next = BTreeMap::new();
    let mut reports = BTreeMap::new();
    let (mut gas_energy, mut mass, mut kinetic, mut exterior, mut applied, mut heat) =
        (Sum::default(), Sum::default(), Sum::default(), Sum::default(), Sum::default(), Sum::default());
    for (&id, old) in pistons {
        check(&mut cancelled)?;
        let mut wall = old.clone();
        wall.gas = remaining.remove(&id).ok_or(Error::Invalid("flow proposal omitted a moving chamber"))?;
        let transferred_mass = wall.gas.mass_kg();
        let drive = drives[&id];
        let report = wall.step(drive.input, drive.solve, &mut cancelled)
            .map_err(|error| piston_error(id, error))?;
        if wall.gas.mass_kg().to_bits() != transferred_mass.to_bits() {
            return Err(Error::Invalid("closed mechanical stage changed transported gas mass"));
        }
        let du = stored_change(&old.gas, &wall.gas);
        gas_energy.add(du)?;
        mass.add(wall.gas.mass_kg() - old.gas.mass_kg())?;
        // Read actual endpoint motion independently of the per-wall receipts.
        kinetic.add((0.5 * old.mass_kg * (wall.velocity_m_s - old.velocity_m_s))
            * (wall.velocity_m_s + old.velocity_m_s))?;
        let dv = wall.gas.volume_m3() - old.gas.volume_m3();
        exterior.add(drive.input.exterior_pressure_pa * dv)?;
        applied.add(drive.input.external_force_n * (dv / old.area_m2))?;
        heat.add(drive.input.heat_j)?;
        if open.contains(&id) {
            record_open_account(old, &mut wall, drive, &report, du, controls.total_energy_tolerance_j)
                .map_err(|error| piston_error(id, error))?;
        }
        next.insert(id, wall); reports.insert(id, report);
    }
    if !remaining.keys().eq(reservoirs.keys()) {
        return Err(Error::Invalid("flow proposal changed the fixed-chamber identity set"));
    }
    for (&id, old) in reservoirs {
        check(&mut cancelled)?;
        gas_energy.add(stored_change(old, &remaining[&id]))?;
        mass.add(remaining[&id].mass_kg() - old.mass_kg())?;
    }
    let mut balance = Sum::default();
    for value in [gas_energy.value, kinetic.value, exterior.value, -applied.value, -heat.value] {
        balance.add(value)?;
    }
    if mass.value.abs() > controls.exchange.mass_tolerance_kg
        || balance.value.abs() > controls.total_energy_tolerance_j {
        return Err(Error::Invalid("final gas/wall endpoints do not close total mass and energy"));
    }
    check(&mut cancelled)?;
    Ok(MovingGasCandidate { pistons: next, reservoirs: remaining,
        step: MovingGasStep { flow, pistons: reports, gas_energy_change_j: gas_energy.value,
            kinetic_energy_change_j: kinetic.value, exterior_pressure_work_j: exterior.value,
            external_force_work_j: applied.value, supplied_heat_j: heat.value,
            total_mass_residual_kg: mass.value, total_energy_residual_j: balance.value } })
}

#[cfg(test)]
mod tests;
