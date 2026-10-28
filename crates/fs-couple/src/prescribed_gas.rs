//! Finite gas networks driven by prescribed chamber volumes.
//!
//! A crank, compressor drive or moving geometry supplies the volume endpoint;
//! the original gas material supplies pressure and the work returned to that
//! drive. Every nozzle first sees the same accepted network. The original
//! midpoint heat/volume owner then advances each post-transfer parcel. This
//! flow-then-motion split is first order overall and advances one physical
//! interval. There is no independent piston inertia or duplicated gas law.
//!
//! Transported upstream total enthalpy already includes flow work and the jet
//! energy thermalized by the receiving volume. It is never supplied again as
//! heat. Formation energy retains its original species reference. Gross gas
//! work, exterior-pressure work and net mechanical drive work stay separate.

use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::midpoint::GasHeatVolumeStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::{
    self, SnapshotNozzleError, SnapshotNozzleLink, SnapshotNozzleStep,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::{
    GasSnapshotSolve, MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS,
};
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use std::collections::BTreeMap;

/// Declared chemistry in open chambers with prescribed mechanical drives.
pub mod reactive;

/// Explicit volume and external heat for one chamber over the common interval.
/// Held reservoirs use their current volume. Zero heat is adiabatic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedGasDrive {
    /// Positive endpoint volume [m3], supplied by the mechanical owner.
    pub final_volume_m3: f64,
    /// Signed external heat INTO the gas [J], excluding transported enthalpy.
    pub heat_j: f64,
    /// Absolute exterior pressure held over this interval [Pa]; zero is vacuum.
    /// It affects the net drive work, never the chamber's internal-energy law.
    pub exterior_pressure_pa: f64,
}

/// Original flow admission and independent enclosing energy acceptance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedGasSolve {
    /// Original simultaneous mass/species/enthalpy and gross-withdrawal bounds.
    pub exchange: GasSnapshotSolve,
    /// Absolute energy tolerance [J], applied locally and to the whole network.
    pub total_energy_tolerance_j: f64,
}

/// No physical owner is returned or modified on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum PrescribedGasError {
    /// Missing identities, invalid geometry, controls or extensive arithmetic.
    Invalid(&'static str),
    /// The original same-snapshot nozzle/transport owner refused the interval.
    Flow(SnapshotNozzleError),
    /// A selected chamber's original chemical owner refused its mixed state.
    Reaction {
        id: u64,
        error: crate::cavity_phs::piston::reactive::GasReactionError,
    },
    /// Actual whole-network species change disagrees with chemical stoichiometry.
    SpeciesBalance { maximum_residual_mol: f64 },
    /// The original caloric/EOS/transport or heat-volume owner refused a chamber.
    Parcel { id: u64, error: GasEnergyError },
    /// An actual endpoint fails its complete open-volume energy balance.
    EnergyBalance { id: u64, residual_j: f64 },
    /// The original finite inventories and external exchanges do not close.
    NetworkBalance {
        mass_residual_kg: f64,
        energy_residual_j: f64,
    },
    /// Cancelled before returning all endpoints.
    Cancelled,
}
impl core::fmt::Display for PrescribedGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "prescribed-volume gas network: {self:?}")
    }
}
impl std::error::Error for PrescribedGasError {}

/// One open-volume balance with the actual work-conjugate material pressure.
#[derive(Clone, Debug, PartialEq)]
pub struct PrescribedGasParcelStep {
    /// Original simultaneous heat/volume update, including gross gas work [J].
    pub material: GasHeatVolumeStep,
    /// Original source-bearing midpoint state whose pressure performed the work.
    pub midpoint: ResolvedThermochemicalGasState,
    /// Actual full-interval mass change, including all inlets/outlets [kg].
    pub mass_change_kg: f64,
    /// Actual full-interval internal-energy change, including transported energy [J].
    pub energy_change_j: f64,
    /// Incoming minus outgoing upstream total enthalpy [J], counted once.
    pub net_enthalpy_into_gas_j: f64,
    /// Work delivered against the declared exterior pressure [J].
    pub exterior_pressure_work_j: f64,
    /// Gross gas work minus exterior-pressure work [J]; positive reaches the drive.
    /// For a crank interval dtheta, work-conjugate average torque is this / dtheta.
    pub net_drive_work_j: f64,
    /// Actual dU + gross gas work - external heat - net incoming enthalpy [J].
    pub energy_residual_j: f64,
}

/// Complete same-interval flow, thermodynamics and mechanical-energy response.
#[derive(Clone, Debug, PartialEq)]
pub struct PrescribedGasStep {
    /// Physical interval used by every original nozzle [s].
    pub duration_s: f64,
    /// Original shared-inventory transfer and source-resolved throat reports.
    pub flow: SnapshotNozzleStep,
    /// Local heat/work and actual open-volume balances in chamber-id order.
    pub parcels: BTreeMap<u64, PrescribedGasParcelStep>,
    /// Actual sum of final minus initial gas internal energies [J].
    pub gas_energy_change_j: f64,
    /// Sum of external heat into all volumes [J].
    pub supplied_heat_j: f64,
    /// Sum of midpoint gas pressure times actual volume change [J].
    pub gross_gas_work_j: f64,
    /// Sum of exterior pressure times actual volume change [J].
    pub exterior_pressure_work_j: f64,
    /// Sum of work delivered to the prescribed mechanical drives [J].
    pub net_drive_work_j: f64,
    /// Independently measured net mass change across final physical owners [kg].
    pub total_mass_residual_kg: f64,
    /// Original shared-flow species residual [kg]. The heat/volume operation
    /// changes no species. The opt-in reaction path separately audits its
    /// stoichiometric changes over the complete interval.
    pub maximum_species_mass_residual_kg: f64,
    /// Actual dU_gases + gross gas work - external heat [J].
    pub total_energy_residual_j: f64,
}

/// All proposed gas owners together. The caller publishes only after accepting
/// the returned pressure work at its prescribed mechanical boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct PrescribedGasCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: PrescribedGasStep,
}
impl PrescribedGasCandidate {
    /// Complete proposed media, preserving original caloric and transport sources.
    #[must_use]
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    /// Actual physical balances and source evaluations.
    #[must_use]
    pub const fn step(&self) -> &PrescribedGasStep {
        &self.step
    }
    /// Accept every gas endpoint together without a second material query.
    #[must_use]
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, PrescribedGasStep) {
        (self.parcels, self.step)
    }
}

use PrescribedGasError as Error;
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn material_error(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Parcel { id, error }
    }
}

// Numerical accumulation only; no alternate thermodynamic state or energy law.
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), Error> {
        let adjusted = value - self.correction;
        let next = self.value + adjusted;
        let correction = (next - self.value) - adjusted;
        if ![value, adjusted, next, correction]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(Error::Invalid(
                "gas network extensive sum is unrepresentable",
            ));
        }
        self.value = next;
        self.correction = correction;
        Ok(())
    }
}

/// Drive finite, well-mixed chambers through prescribed volume/heat endpoints.
///
/// `endpoints` and `drives` cover every parcel exactly, including held reservoirs.
/// Every restriction uses the SAME accepted material snapshot. Its conservative
/// update combines concurrent inlets and outlets before reconstructing a volume.
/// Each resulting parcel then uses its original nonlinear midpoint heat/volume
/// law. Its actual p_mid * dV is work OUT of the gas, including compression's
/// negative work. The drive must absorb or supply that work; no fictitious free
/// piston or mechanical inertia is introduced to prescribe the motion.
///
/// The whole network obeys dU + W_gross - Q_external = 0. Transported total
/// enthalpy cancels internally. Exterior pressure partitions W_gross into work
/// on the exterior and net work to the drive, without changing the gas state.
/// No flow work, jet kinetic energy or reaction heat is added a second time.
///
/// The combined open-system update is FIRST ORDER. Source domains, donor
/// withdrawal and pre-motion pressure-crossing checks remain with their original
/// owners. Motion may reverse a pressure ordering; the next interval reselects
/// upwind flow. Refine time and split prescribed valve/volume/heat events. No
/// event locator, stiff implicit flow, exact isentropy or combustion is claimed.
/// All trials, including late caloric failures, preserve every accepted input.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[SnapshotNozzleLink],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    drives: &BTreeMap<u64, PrescribedGasDrive>,
    controls: PrescribedGasSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedGasCandidate, Error> {
    staged_candidate(
        parcels,
        links,
        duration_s,
        endpoints,
        drives,
        controls,
        |_, _| Ok(()),
        cancelled,
    )
    .map(|(candidate, ())| candidate)
}

// The optional chemical stage shares every original flow, heat/volume and
// original-to-final balance. It cannot advance a second mechanical clock.
#[allow(clippy::too_many_arguments)]
fn staged_candidate<R, C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[SnapshotNozzleLink],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    drives: &BTreeMap<u64, PrescribedGasDrive>,
    controls: PrescribedGasSolve,
    stage: impl FnOnce(&mut BTreeMap<u64, ClosedGasParcel>, &mut C) -> Result<R, Error>,
    mut cancelled: C,
) -> Result<(PrescribedGasCandidate, R), Error> {
    poll(&mut cancelled)?;
    if parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
        || drives.len() != parcels.len()
        || endpoints.len() != parcels.len()
        || !positive(duration_s)
        || !positive(controls.total_energy_tolerance_j)
        || !parcels.keys().eq(drives.keys())
        || !parcels.keys().eq(endpoints.keys())
    {
        return Err(Error::Invalid(
            "bounded nonempty gas graph, positive interval/tolerance and exact chamber drive/solver maps required",
        ));
    }
    for drive in drives.values() {
        poll(&mut cancelled)?;
        if !positive(drive.final_volume_m3)
            || !drive.heat_j.is_finite()
            || !drive.exterior_pressure_pa.is_finite()
            || drive.exterior_pressure_pa < 0.0
        {
            return Err(Error::Invalid(
                "positive volume, finite heat and nonnegative absolute exterior pressure required",
            ));
        }
    }
    let flow = snapshot::step_candidate(
        parcels,
        links,
        duration_s,
        endpoints,
        controls.exchange,
        &mut cancelled,
    )
    .map_err(|error| {
        if error == SnapshotNozzleError::Cancelled {
            Error::Cancelled
        } else {
            Error::Flow(error)
        }
    })?;
    let (mut post_flow, flow) = flow.into_parts();
    let staged = stage(&mut post_flow, &mut cancelled)?;
    let mut enthalpy: BTreeMap<u64, Sum> = parcels.keys().map(|&id| (id, Sum::default())).collect();
    for flux in &flow.exchange.fluxes {
        poll(&mut cancelled)?;
        enthalpy
            .get_mut(&flux.transfer.donor)
            .expect("original flow donor")
            .add(-flux.enthalpy_j)?;
        enthalpy
            .get_mut(&flux.transfer.receiver)
            .expect("original flow receiver")
            .add(flux.enthalpy_j)?;
    }
    let mut next = BTreeMap::new();
    let mut reports = BTreeMap::new();
    let (mut energy, mut mass, mut heat) = (Sum::default(), Sum::default(), Sum::default());
    let (mut gross, mut outside, mut net) = (Sum::default(), Sum::default(), Sum::default());
    for (&id, old) in parcels {
        poll(&mut cancelled)?;
        let drive = drives[&id];
        let trial = post_flow[&id]
            .midpoint_heat_volume(
                drive.final_volume_m3,
                drive.heat_j,
                endpoints[&id],
                &mut cancelled,
            )
            .map_err(|error| material_error(id, error))?;
        let midpoint = trial.midpoint().clone();
        let (gas, material) = trial.into_parts();
        let dm = gas.mass_kg() - old.mass_kg();
        // Subtract specific storage before multiplying, retaining formation
        // references even when the chamber's total inventory changes.
        let du = gas.mass_kg()
            * (gas.state().specific_internal_energy_j_kg()
                - old.state().specific_internal_energy_j_kg())
            + dm * old.state().specific_internal_energy_j_kg();
        let exterior_work = drive.exterior_pressure_pa * (drive.final_volume_m3 - old.volume_m3());
        let net_work = material.work_by_gas_j - exterior_work;
        let incoming = enthalpy[&id].value;
        let mut balance = Sum::default();
        for value in [du, material.work_by_gas_j, -drive.heat_j, -incoming] {
            balance.add(value)?;
        }
        if balance.value.abs() > controls.total_energy_tolerance_j {
            return Err(Error::EnergyBalance {
                id,
                residual_j: balance.value,
            });
        }
        energy.add(du)?;
        mass.add(dm)?;
        heat.add(drive.heat_j)?;
        gross.add(material.work_by_gas_j)?;
        outside.add(exterior_work)?;
        net.add(net_work)?;
        reports.insert(
            id,
            PrescribedGasParcelStep {
                material,
                midpoint,
                mass_change_kg: dm,
                energy_change_j: du,
                net_enthalpy_into_gas_j: incoming,
                exterior_pressure_work_j: exterior_work,
                net_drive_work_j: net_work,
                energy_residual_j: balance.value,
            },
        );
        next.insert(id, gas);
    }
    let mut balance = Sum::default();
    for value in [energy.value, gross.value, -heat.value] {
        balance.add(value)?;
    }
    if mass.value.abs() > controls.exchange.mass_tolerance_kg
        || balance.value.abs() > controls.total_energy_tolerance_j
    {
        return Err(Error::NetworkBalance {
            mass_residual_kg: mass.value,
            energy_residual_j: balance.value,
        });
    }
    let maximum_species_mass_residual_kg = flow.exchange.maximum_species_mass_residual_kg;
    poll(&mut cancelled)?;
    Ok((
        PrescribedGasCandidate {
            parcels: next,
            step: PrescribedGasStep {
                duration_s,
                flow,
                parcels: reports,
                gas_energy_change_j: energy.value,
                supplied_heat_j: heat.value,
                gross_gas_work_j: gross.value,
                exterior_pressure_work_j: outside.value,
                net_drive_work_j: net.value,
                total_mass_residual_kg: mass.value,
                maximum_species_mass_residual_kg,
                total_energy_residual_j: balance.value,
            },
        },
        staged,
    ))
}

#[cfg(test)]
mod tests;
