//! Conservative species transport through shock-carrying branches and finite plenums.
//!
//! Each branch and compact volume retains the same original species bank and
//! per-species energy references, while its actual composition evolves. One
//! characteristic port flux supplies opposite exchanges to its two owners.
//! Network plenums are adiabatic and perfectly mixed during a network transport
//! step. They may remain rigid or follow prescribed volume rates, with pressure
//! work integrated in the same stages as port transport. The `moving` owner can
//! also construct a quasistatic adiabatic volume-change candidate between accepted steps;
//! it accounts boundary work without disguising it as heat or a pressure source.
//! Incoming kinetic energy remains in plenum internal-energy balance. Stationary
//! supports receive each port's separately reported axial impulse.
//!
//! This estimated inviscid model does not resolve three-dimensional junction
//! mixing, molecular diffusion, wall losses or chemical reaction. All branches
//! and plenums share one atomic SSP RK2 step and accepted clock.

use super::evolve::{Stage, totals};
use super::{
    Cell, Conserved, Error, FaceFlux, MAX_COMPONENTS, MAX_SPECIES, MixtureEuler1d, MixtureGas,
    Primitive, StepReport, Totals, allocated, poll, positive, sum,
};
pub use crate::dynamics::network::{
    BranchEnd, MAX_NETWORK_BRANCHES, MAX_NETWORK_PLENUMS, PlenumConnection,
};
use crate::dynamics::{Boundaries, MAX_EULER_CELLS};
use fs_qty::Time;

mod coupled;
mod step;

/// Reversible adiabatic compact-volume changes with explicit boundary work.
pub mod moving;

/// Controllable, conservative valve links between finite rigid plenums.
pub mod valve;

/// Conservative molecular species, heat and viscous transport through finite ports.
pub mod dissipation;

/// Atomic molecular/advection/molecular composition of a complete network.
pub mod dissipating;

/// Atomic finite-rate chemistry in every branch cell and compact gas volume.
pub mod reactive;

/// Clock-preserving radial drag and conservative heat exchange with wall storage.
pub mod wall_exchange;

/// Finite, perfectly mixed gas volume with evolving species and zero bulk momentum.
#[derive(Clone, Debug, PartialEq)]
pub struct Plenum {
    volume_m3: f64,
    state: Conserved,
}
impl Plenum {
    /// Retain a finite volume and zero-momentum state in canonical species order.
    /// The network admits its caloric state against the exact shared source bank.
    pub fn new(volume_m3: f64, state: Conserved) -> Result<Self, Error> {
        if !positive(volume_m3)
            || !(1..=MAX_SPECIES).contains(&state.partial_densities.len())
            || state
                .partial_densities
                .iter()
                .any(|&x| !x.is_finite() || x < 0.0)
            || state.momentum != 0.0
            || !state.energy.is_finite()
            || !positive(volume_m3 * state.as_cell().density())
            || !(volume_m3 * state.energy).is_finite()
            || state
                .partial_densities
                .iter()
                .any(|&x| !(volume_m3 * x).is_finite())
        {
            return Err(Error::Invalid(
                "mixture plenum requires finite positive volume and stationary gas",
            ));
        }
        Ok(Self { volume_m3, state })
    }
    /// Current physical volume [m³].
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Borrow current species densities and referenced internal-energy density.
    pub fn state(&self) -> Cell<'_> {
        self.state.as_cell()
    }
}

/// Extensive species and energy storage, without adding unrelated branch axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkTotals {
    /// Canonical species masses [kg]; inactive entries are exactly zero.
    pub species_masses: [f64; MAX_SPECIES],
    /// Active species count.
    pub species_count: usize,
    /// Branch total energy plus compact internal energy in species gauges [J].
    pub energy: f64,
}
impl NetworkTotals {
    /// Total mass derived from the retained species inventories [kg].
    pub fn mass(&self) -> f64 {
        sum(self.species_masses[..self.species_count].iter().copied())
    }
    /// Physical energy including the original per-species formation references [J].
    pub fn absolute_energy(&self, gas: &MixtureGas) -> Result<f64, Error> {
        if self.species_count != gas.species_count() {
            return Err(Error::Invalid("network total/source species axis"));
        }
        let energy = self.energy
            + sum(self.species_masses[..self.species_count]
                .iter()
                .zip(gas.energy_references())
                .map(|(a, b)| a * b));
        if !energy.is_finite() {
            return Err(Error::Invalid("non-finite physical network energy"));
        }
        Ok(energy)
    }
}

/// Accepted exchange into one plenum, in the retained connection order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PortExchange {
    /// Mass of each canonical species delivered to the plenum [kg].
    pub species_masses_to_plenum: [f64; MAX_SPECIES],
    /// Active species count; inactive mass entries are exactly zero.
    pub species_count: usize,
    /// Impulse on the fixed support in this branch's axial coordinate [kg m/s].
    pub momentum_to_plenum: f64,
    /// Referenced total energy delivered to the plenum, including incoming kinetic energy [J].
    pub energy_to_plenum: f64,
}
impl PortExchange {
    /// Total mass delivered to the plenum [kg]; negative denotes reverse flow.
    pub fn mass_to_plenum(&self) -> f64 {
        sum(self.species_masses_to_plenum[..self.species_count]
            .iter()
            .copied())
    }
    /// Physical energy exchange with each transported species reference included once.
    pub fn absolute_energy_to_plenum(&self, gas: &MixtureGas) -> Result<f64, Error> {
        NetworkTotals {
            species_masses: self.species_masses_to_plenum,
            species_count: self.species_count,
            energy: self.energy_to_plenum,
        }
        .absolute_energy(gas)
    }
}

/// Accepted mechanical boundary exchange for one compact plenum.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PlenumBoundaryWork {
    /// Physical volume at the beginning of the accepted step [m³].
    pub initial_volume_m3: f64,
    /// Physical volume at the end of the accepted step [m³].
    pub final_volume_m3: f64,
    /// Mechanical work done ON the gas [J], positive under compression.
    /// Uses the same pressure quadrature as the extensive-energy update.
    pub work_on_gas_j: f64,
}

/// Accepted common step; internal species and energy exchanges cancel globally.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkStepReport {
    /// Actual elapsed duration, no greater than the requested maximum.
    pub elapsed: Time,
    /// Maximum branch/port acoustic or compact-volume fractional-change CFL over both stages.
    pub maximum_cfl: f64,
    /// Number of globally rejected and halved trials.
    pub halvings: usize,
    /// Accepted interfaces using HLL fallback, counted once per stage.
    pub hll_fallback_faces: usize,
    /// Accepted branch cells with reduced reconstruction slopes.
    pub limited_cells: usize,
    /// Whether a numerical refusal switched requested MUSCL reconstruction to first order.
    pub first_order_retry: bool,
    /// Original combined storage.
    pub before: NetworkTotals,
    /// Accepted combined storage.
    pub after: NetworkTotals,
    /// Accepted inflow through external branch ends only.
    pub boundary_change: NetworkTotals,
    /// Mechanical work done ON all plenums [J], separate from transported energy.
    pub boundary_work_on_gas_j: f64,
    /// `after - before - boundary_change - boundary_work`, per species and in the
    /// original energy gauges. Mechanical work changes energy only.
    pub conservation_residual: NetworkTotals,
}

struct NetworkWorkspace {
    plenum_primitive: Vec<Primitive>,
    plenum_stage: Vec<f64>,
    plenum_candidate: Vec<f64>,
    plenum_rhs: Vec<f64>,
    plenum_compensation: Vec<f64>,
    plenum_signal_rates: Vec<f64>,
    endpoint_values: Vec<f64>,
    endpoint_data: Vec<[Option<FaceFlux>; 2]>,
    stages_a: Vec<Stage>,
    stages_b: Vec<Stage>,
    branch_before: Vec<Totals>,
    branch_report_candidate: Vec<StepReport>,
    port_exchange_candidate: Vec<PortExchange>,
    plenum_work_candidate: Vec<PlenumBoundaryWork>,
}

/// Shared-clock species-bearing network of original Euler branches and finite plenums.
///
/// Sources, canonical species order, caloric intervals, per-species energy gauges
/// and accepted clocks must match exactly. Local composition and temperature may
/// differ. Immutable branch access prevents separately advancing a network owner.
/// Successful numerical steps reuse construction-time buffers without allocation.
pub struct MixtureNetwork1d {
    branches: Vec<MixtureEuler1d>,
    plenums: Vec<Plenum>,
    connections: Vec<PlenumConnection>,
    endpoint_connections: Vec<[Option<usize>; 2]>,
    work: NetworkWorkspace,
    port_exchanges: Vec<PortExchange>,
    plenum_work: Vec<PlenumBoundaryWork>,
    branch_reports: Vec<StepReport>,
    last_report: Option<NetworkStepReport>,
}
impl MixtureNetwork1d {
    /// Connect physical branch ends to finite, perfectly mixed gas volumes.
    ///
    /// The combined branch-cell/plenum count is bounded by `MAX_EULER_CELLS`.
    /// Periodic branches may remain isolated but cannot have connected ends.
    /// Each endpoint may appear once; an unconnected rigid plenum retains its state.
    pub fn new(
        branches: Vec<MixtureEuler1d>,
        plenums: Vec<Plenum>,
        connections: Vec<PlenumConnection>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !(1..=MAX_NETWORK_BRANCHES).contains(&branches.len())
            || plenums.len() > MAX_NETWORK_PLENUMS
            || connections.len() > 2 * branches.len()
        {
            return Err(Error::Invalid(
                "mixture network branch, plenum or connection count",
            ));
        }
        let gas = branches[0].gas();
        let n = gas.species_count();
        let m = n + 2;
        let clock = branches[0].time().value();
        let mut cells = plenums.len();
        let mut storage = Accumulator::new(n);
        for branch in &branches {
            poll(&mut cancelled)?;
            if branch.gas() != gas {
                return Err(Error::Invalid(
                    "mixture network branches must retain the identical species source bank",
                ));
            }
            if branch.time().value().to_bits() != clock.to_bits() {
                return Err(Error::Invalid("mixture network branch clocks differ"));
            }
            cells = cells
                .checked_add(branch.cell_count())
                .filter(|&x| x <= MAX_EULER_CELLS)
                .ok_or(Error::Invalid(
                    "mixture network aggregate storage exceeds MAX_EULER_CELLS",
                ))?;
            let value = totals(&branch.states, branch.mesh.volumes(), n, &mut cancelled)?;
            storage.add(&value.species_masses[..n], value.energy)?;
        }
        let mut endpoint_connections = allocated(branches.len(), [None; 2])?;
        for (port, connection) in connections.iter().enumerate() {
            poll(&mut cancelled)?;
            if connection.branch >= branches.len() || connection.plenum >= plenums.len() {
                return Err(Error::Invalid(
                    "mixture network connection index outside topology",
                ));
            }
            if branches[connection.branch].boundaries == Boundaries::Periodic {
                return Err(Error::Invalid(
                    "a periodic mixture branch has no connectable physical endpoints",
                ));
            }
            let slot = &mut endpoint_connections[connection.branch][connection.end.index()];
            if slot.is_some() {
                return Err(Error::Invalid(
                    "a mixture branch endpoint cannot connect more than once",
                ));
            }
            *slot = Some(port);
        }
        let mut plenum_primitive = allocated(plenums.len(), Primitive::default())?;
        for (i, plenum) in plenums.iter().enumerate() {
            plenum_primitive[i] = gas.primitive(plenum.state(), &mut cancelled)?;
            storage.add_cell(plenum.state(), plenum.volume_m3)?;
        }
        let plenum_size = plenums.len().checked_mul(m).ok_or(Error::Allocation)?;
        let endpoint_size = branches.len().checked_mul(2 * m).ok_or(Error::Allocation)?;
        let work = NetworkWorkspace {
            plenum_primitive,
            plenum_stage: allocated(plenum_size, 0.0)?,
            plenum_candidate: allocated(plenum_size, 0.0)?,
            plenum_rhs: allocated(plenum_size, 0.0)?,
            plenum_compensation: allocated(plenum_size, 0.0)?,
            plenum_signal_rates: allocated(plenums.len(), 0.0)?,
            endpoint_values: allocated(endpoint_size, 0.0)?,
            endpoint_data: allocated(branches.len(), [None; 2])?,
            stages_a: allocated(branches.len(), empty_stage())?,
            stages_b: allocated(branches.len(), empty_stage())?,
            branch_before: allocated(branches.len(), empty_totals(n))?,
            branch_report_candidate: allocated(branches.len(), empty_report(n))?,
            port_exchange_candidate: allocated(connections.len(), empty_exchange(n))?,
            plenum_work_candidate: allocated(plenums.len(), PlenumBoundaryWork::default())?,
        };
        let port_exchanges = allocated(connections.len(), empty_exchange(n))?;
        let plenum_work = allocated(plenums.len(), PlenumBoundaryWork::default())?;
        let branch_reports = allocated(branches.len(), empty_report(n))?;
        poll(&mut cancelled)?;
        Ok(Self {
            branches,
            plenums,
            connections,
            endpoint_connections,
            work,
            port_exchanges,
            plenum_work,
            branch_reports,
            last_report: None,
        })
    }
    /// Original common source bank, species axis, caloric path and energy gauges.
    pub fn gas(&self) -> &MixtureGas {
        self.branches[0].gas()
    }
    /// Common accepted physical time.
    pub fn time(&self) -> Time {
        self.branches[0].time()
    }
    /// Original branch order and accepted physical fields.
    pub fn branches(&self) -> &[MixtureEuler1d] {
        &self.branches
    }
    /// Original compact-volume order and current storage.
    pub fn plenums(&self) -> &[Plenum] {
        &self.plenums
    }
    /// Deterministic retained connection order for port reductions and reports.
    pub fn connections(&self) -> &[PlenumConnection] {
        &self.connections
    }
    /// Most recent accepted report; refusal/cancellation leaves it unchanged.
    pub const fn last_report(&self) -> Option<NetworkStepReport> {
        self.last_report
    }
    /// Accepted exchanges in original connection order, or None before stepping.
    pub fn accepted_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.map(|_| self.port_exchanges.as_slice())
    }
    /// Accepted volumes and mechanical work in original plenum order.
    /// Rigid-volume steps report identical endpoint volumes and exactly zero work.
    pub fn accepted_plenum_boundary_work(&self) -> Option<&[PlenumBoundaryWork]> {
        self.last_report.map(|_| self.plenum_work.as_slice())
    }
    /// Complete branch species/momentum/energy balances for the accepted step.
    pub fn accepted_branch_reports(&self) -> Option<&[StepReport]> {
        self.last_report.map(|_| self.branch_reports.as_slice())
    }
}

struct Accumulator {
    values: [f64; MAX_SPECIES + 1],
    corrections: [f64; MAX_SPECIES + 1],
    n: usize,
}
impl Accumulator {
    fn new(n: usize) -> Self {
        Self {
            values: [0.0; MAX_SPECIES + 1],
            corrections: [0.0; MAX_SPECIES + 1],
            n,
        }
    }
    fn add(&mut self, species: &[f64], energy: f64) -> Result<(), Error> {
        for (k, value) in species
            .iter()
            .copied()
            .chain(core::iter::once(energy))
            .enumerate()
        {
            let y = value - self.corrections[k];
            let next = self.values[k] + y;
            if !value.is_finite() || !next.is_finite() {
                return Err(Error::Invalid(
                    "non-finite extensive mixture network balance",
                ));
            }
            self.corrections[k] = (next - self.values[k]) - y;
            self.values[k] = next;
        }
        Ok(())
    }
    fn add_cell(&mut self, state: Cell<'_>, volume: f64) -> Result<(), Error> {
        let mut species = [0.0; MAX_SPECIES];
        for (value, &density) in species.iter_mut().zip(state.partial_densities) {
            *value = volume * density;
        }
        self.add(&species[..self.n], volume * state.energy)
    }
    fn totals(&self) -> NetworkTotals {
        let mut species_masses = [0.0; MAX_SPECIES];
        species_masses[..self.n].copy_from_slice(&self.values[..self.n]);
        NetworkTotals {
            species_masses,
            species_count: self.n,
            energy: self.values[self.n],
        }
    }
}
fn empty_totals(n: usize) -> Totals {
    Totals::from_values(&[0.0; MAX_COMPONENTS][..n + 2], n)
}
fn empty_stage() -> Stage {
    Stage {
        left: [0.0; MAX_COMPONENTS],
        right: [0.0; MAX_COMPONENTS],
        signal_rate: 0.0,
        area_force: 0.0,
        fallbacks: 0,
        limited: 0,
    }
}
fn empty_exchange(n: usize) -> PortExchange {
    PortExchange {
        species_masses_to_plenum: [0.0; MAX_SPECIES],
        species_count: n,
        momentum_to_plenum: 0.0,
        energy_to_plenum: 0.0,
    }
}
fn empty_report(n: usize) -> StepReport {
    StepReport {
        elapsed: Time::new(0.0),
        maximum_cfl: 0.0,
        halvings: 0,
        hll_fallback_faces: 0,
        limited_cells: 0,
        first_order_retry: false,
        before: empty_totals(n),
        after: empty_totals(n),
        boundary_change: empty_totals(n),
        left_boundary_change: empty_totals(n),
        right_boundary_change: empty_totals(n),
        area_momentum_impulse: 0.0,
        conservation_residual: empty_totals(n),
    }
}

#[cfg(test)]
mod tests;
