//! Shock-carrying passages coupled through finite, rigid, perfectly mixed plenums.
//!
//! Every port uses one source-owned characteristic gas flux for both owners.
//! Branches resolve axial momentum and total energy; a compact plenum stores
//! mass and internal energy and has no resolved bulk momentum. Incoming kinetic
//! energy remains in its energy balance. The stationary plenum support receives
//! the reported port impulse without doing mechanical work.
//!
//! All branches and plenums advance in one transactional SSP RK2 step. This is
//! an inviscid, frozen-composition engineering model: a finite plenum is not a
//! zero-volume junction, resolved three-dimensional mixing region, or a proof
//! of entropy stability. Different branch axes are deliberately not added into
//! a fictitious scalar network momentum.

use super::evolve::Stage;
use super::{
    Boundaries, Conserved, Error, Euler1d, FaceFlux, FrozenGas, MAX_EULER_CELLS, Primitive,
    StepReport, Totals, allocated, poll, positive,
};
use fs_qty::Time;

mod step;

/// Explicit topology bound, in addition to the aggregate cell-storage bound.
pub const MAX_NETWORK_BRANCHES: usize = 4096;
/// Explicit compact-volume bound, in addition to the aggregate storage bound.
pub const MAX_NETWORK_PLENUMS: usize = 4096;

/// One physical endpoint in a branch's increasing axial coordinate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BranchEnd {
    /// Face at the smallest axial position; outward normal is -1.
    Left,
    /// Face at the largest axial position; outward normal is +1.
    Right,
}
impl BranchEnd {
    pub(crate) const fn index(self) -> usize {
        match self {
            Self::Left => 0,
            Self::Right => 1,
        }
    }
    pub(crate) const fn normal(self) -> f64 {
        match self {
            Self::Left => -1.0,
            Self::Right => 1.0,
        }
    }
}

/// A branch endpoint connected to one compact plenum.
///
/// Indices refer to the constructor's retained branch and plenum order. Each
/// branch endpoint may occur at most once. This connection replaces that end's
/// original physical boundary; unconnected ends keep their original condition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlenumConnection {
    /// Index of the original branch.
    pub branch: usize,
    /// Left or right physical face.
    pub end: BranchEnd,
    /// Index of the receiving compact volume.
    pub plenum: usize,
}

/// Finite, stationary, adiabatic, well-mixed gas volume.
///
/// The state uses the network's common frozen mixture and energy reference.
/// Its momentum is exactly zero; `state.energy` is referenced internal-energy
/// density. This model assumes that the plenum is acoustically compact on the
/// time and length scales of interest. It does not resolve its internal waves.
#[derive(Clone, Debug, PartialEq)]
pub struct Plenum {
    volume_m3: f64,
    state: Conserved,
}
impl Plenum {
    /// Construct a compact volume in m³ with a zero-velocity conservative state.
    ///
    /// Obtain the state from the intended network gas, for example
    /// `gas.conserved(density, 0.0, pressure)`. The network constructor performs
    /// complete caloric admission with that exact source before any evolution.
    pub fn new(volume_m3: f64, state: Conserved) -> Result<Self, Error> {
        if !positive(volume_m3)
            || !state.finite()
            || !positive(state.density)
            || state.momentum != 0.0
            || !positive(volume_m3 * state.density)
            || !(volume_m3 * state.energy).is_finite()
        {
            return Err(Error::Invalid(
                "plenum requires finite positive volume and stationary gas",
            ));
        }
        Ok(Self { volume_m3, state })
    }
    /// Original rigid volume [m³].
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Current density and referenced internal-energy density; momentum is zero.
    pub const fn state(&self) -> Conserved {
        self.state
    }
}

/// Extensive network storage. All branches share one exact energy reference.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkTotals {
    /// Total branch plus plenum mass [kg].
    pub mass: f64,
    /// Branch total energy plus plenum internal energy in the common gauge [J].
    pub energy: f64,
}
impl NetworkTotals {
    /// Recover physical energy using the original signed formation reference.
    pub fn absolute_energy(self, gas: &FrozenGas) -> f64 {
        self.energy + gas.energy_reference() * self.mass
    }
}

/// One accepted, time-integrated exchange in the retained connection order.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PortExchange {
    /// Mass delivered from the branch to its plenum [kg]; negative is reverse flow.
    pub mass_to_plenum: f64,
    /// Impulse delivered to the fixed plenum support [kg m/s], in this branch's
    /// increasing axial coordinate. The branch receives the exact opposite term.
    pub momentum_to_plenum: f64,
    /// Referenced total energy delivered to the plenum [J], including kinetic
    /// energy of incoming flow. Physical energy adds `e_ref * mass_to_plenum`.
    pub energy_to_plenum: f64,
}

/// Accepted common step. Internal port exchanges cancel from the global balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkStepReport {
    /// Actual common elapsed duration, no greater than the requested maximum.
    pub elapsed: Time,
    /// Largest branch or plenum face-speed CFL over both accepted stages.
    pub maximum_cfl: f64,
    /// Number of rejected, halved trials before acceptance.
    pub halvings: usize,
    /// Accepted branch interfaces using the HLL fallback, counted once per stage.
    pub hll_fallback_faces: usize,
    /// Accepted branch cells with source-admission slope reduction.
    pub limited_cells: usize,
    /// Whether a numerical refusal changed the requested MUSCL step to first order.
    pub first_order_retry: bool,
    /// Original combined storage.
    pub before: NetworkTotals,
    /// Accepted combined storage.
    pub after: NetworkTotals,
    /// Accepted inflow through external branch ends only. Connected ports are internal.
    pub boundary_change: NetworkTotals,
    /// `after - before - boundary_change`, in the common energy gauge.
    pub conservation_residual: NetworkTotals,
}

struct NetworkWorkspace {
    plenum_primitive: Vec<Primitive>,
    plenum_stage: Vec<Conserved>,
    plenum_candidate: Vec<Conserved>,
    plenum_rhs: Vec<[f64; 2]>,
    plenum_compensation: Vec<[f64; 2]>,
    plenum_signal_rates: Vec<f64>,
    endpoint_fluxes: Vec<[Option<FaceFlux>; 2]>,
    stages_a: Vec<Stage>,
    stages_b: Vec<Stage>,
    branch_before: Vec<Totals>,
    branch_report_candidate: Vec<StepReport>,
    port_exchange_candidate: Vec<PortExchange>,
}

/// A shared-clock network of existing Euler branches and finite rigid plenums.
///
/// Original source identity, composition, caloric interval and energy gauge must
/// match exactly across branches. Accepted branch clocks must also match. Each
/// successful call advances every owner at the same two RK stages and commits
/// all states, clocks and reports together. There is no independently mutable
/// branch handle while the network owns the field.
pub struct GasNetwork1d {
    branches: Vec<Euler1d>,
    plenums: Vec<Plenum>,
    connections: Vec<PlenumConnection>,
    endpoint_connections: Vec<[Option<usize>; 2]>,
    work: NetworkWorkspace,
    port_exchanges: Vec<PortExchange>,
    branch_reports: Vec<StepReport>,
    last_report: Option<NetworkStepReport>,
}
impl GasNetwork1d {
    /// Couple existing stationary passages through a declared finite-volume topology.
    ///
    /// The sum of branch cells and plenums is bounded by `MAX_EULER_CELLS`.
    /// Periodic branches may remain isolated but cannot have connected endpoints.
    /// An unconnected plenum retains its initial state. Every plenum state belongs
    /// to the common frozen source and is admitted before construction succeeds.
    pub fn new(
        branches: Vec<Euler1d>,
        plenums: Vec<Plenum>,
        connections: Vec<PlenumConnection>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !(1..=MAX_NETWORK_BRANCHES).contains(&branches.len())
            || plenums.len() > MAX_NETWORK_PLENUMS
            || connections.len() > 2 * branches.len()
        {
            return Err(Error::Invalid("network branch, plenum or connection count"));
        }
        let gas = branches[0].gas();
        let clock = branches[0].time().value();
        let mut cells = plenums.len();
        for branch in &branches {
            poll(&mut cancelled)?;
            if branch.atmosphere().is_some() {
                return Err(Error::Invalid(
                    "network gravity needs a shared potential datum",
                ));
            }
            if branch.gas() != gas {
                return Err(Error::Invalid(
                    "network branches must retain the identical frozen gas source",
                ));
            }
            if branch.time().value().to_bits() != clock.to_bits() {
                return Err(Error::Invalid("network branch clocks differ"));
            }
            cells = cells
                .checked_add(branch.states().len())
                .filter(|&n| n <= MAX_EULER_CELLS)
                .ok_or(Error::Invalid(
                    "network aggregate storage exceeds MAX_EULER_CELLS",
                ))?;
        }
        let mut endpoint_connections = allocated(branches.len(), [None; 2])?;
        for (port, connection) in connections.iter().enumerate() {
            poll(&mut cancelled)?;
            if connection.branch >= branches.len() || connection.plenum >= plenums.len() {
                return Err(Error::Invalid("network connection index outside topology"));
            }
            if branches[connection.branch].boundaries == Boundaries::Periodic {
                return Err(Error::Invalid(
                    "a periodic branch has no connectable physical endpoints",
                ));
            }
            let slot = &mut endpoint_connections[connection.branch][connection.end.index()];
            if slot.is_some() {
                return Err(Error::Invalid(
                    "a branch endpoint cannot connect more than once",
                ));
            }
            *slot = Some(port);
        }
        let mut plenum_primitive = allocated(plenums.len(), Primitive::default())?;
        let mut storage = Accumulator::default();
        for branch in &branches {
            let value =
                super::evolve::totals(branch.states(), branch.mesh().volumes(), &mut cancelled)?;
            storage.add(value.mass, value.energy)?;
        }
        for (i, plenum) in plenums.iter().enumerate() {
            plenum_primitive[i] = gas.decode(plenum.state, None, &mut cancelled)?;
            storage.add(
                plenum.volume_m3 * plenum.state.density,
                plenum.volume_m3 * plenum.state.energy,
            )?;
        }
        let work = NetworkWorkspace {
            plenum_primitive,
            plenum_stage: allocated(plenums.len(), Conserved::default())?,
            plenum_candidate: allocated(plenums.len(), Conserved::default())?,
            plenum_rhs: allocated(plenums.len(), [0.0; 2])?,
            plenum_compensation: allocated(plenums.len(), [0.0; 2])?,
            plenum_signal_rates: allocated(plenums.len(), 0.0)?,
            endpoint_fluxes: allocated(branches.len(), [None; 2])?,
            stages_a: allocated(branches.len(), empty_stage())?,
            stages_b: allocated(branches.len(), empty_stage())?,
            branch_before: allocated(branches.len(), Totals::default())?,
            branch_report_candidate: allocated(branches.len(), empty_report())?,
            port_exchange_candidate: allocated(connections.len(), PortExchange::default())?,
        };
        let port_exchanges = allocated(connections.len(), PortExchange::default())?;
        let branch_reports = allocated(branches.len(), empty_report())?;
        poll(&mut cancelled)?;
        Ok(Self {
            branches,
            plenums,
            connections,
            endpoint_connections,
            work,
            port_exchanges,
            branch_reports,
            last_report: None,
        })
    }
    /// Exact original gas source, composition, caloric path and energy gauge.
    pub fn gas(&self) -> &FrozenGas {
        self.branches[0].gas()
    }
    /// Accepted common physical time.
    pub fn time(&self) -> Time {
        self.branches[0].time()
    }
    /// Immutable original branch order and current fields.
    pub fn branches(&self) -> &[Euler1d] {
        &self.branches
    }
    /// Immutable original compact-volume order and current storage.
    pub fn plenums(&self) -> &[Plenum] {
        &self.plenums
    }
    /// Retained deterministic connection order used for all port reductions.
    pub fn connections(&self) -> &[PlenumConnection] {
        &self.connections
    }
    /// Most recent accepted common report, or None before the first step.
    pub const fn last_report(&self) -> Option<NetworkStepReport> {
        self.last_report
    }
    /// Accepted exchanges in the retained connection order, or None before stepping.
    /// A failed or cancelled step leaves this slice's values unchanged.
    pub fn accepted_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.map(|_| self.port_exchanges.as_slice())
    }
    /// Each branch's complete mass/momentum/energy report for the accepted common step.
    /// A branch's boundary change includes its connected and external ends.
    pub fn accepted_branch_reports(&self) -> Option<&[StepReport]> {
        self.last_report.map(|_| self.branch_reports.as_slice())
    }
}

#[derive(Default)]
struct Accumulator {
    sum: [f64; 2],
    compensation: [f64; 2],
}
impl Accumulator {
    fn add(&mut self, mass: f64, energy: f64) -> Result<(), Error> {
        for (k, value) in [mass, energy].into_iter().enumerate() {
            let y = value - self.compensation[k];
            let next = self.sum[k] + y;
            if !value.is_finite() || !next.is_finite() {
                return Err(Error::Invalid("non-finite extensive network balance"));
            }
            self.compensation[k] = (next - self.sum[k]) - y;
            self.sum[k] = next;
        }
        Ok(())
    }
    fn totals(self) -> NetworkTotals {
        NetworkTotals {
            mass: self.sum[0],
            energy: self.sum[1],
        }
    }
}
fn empty_stage() -> Stage {
    Stage {
        left: FaceFlux::default(),
        right: FaceFlux::default(),
        signal_rate: 0.0,
        area_force: 0.0,
        gravity_force: 0.0,
        gravity_power: 0.0,
        fallbacks: 0,
        limited: 0,
    }
}
fn empty_report() -> StepReport {
    StepReport {
        elapsed: Time::new(0.0),
        maximum_cfl: 0.0,
        halvings: 0,
        hll_fallback_faces: 0,
        limited_cells: 0,
        first_order_retry: false,
        before: Totals::default(),
        after: Totals::default(),
        boundary_change: Totals::default(),
        left_boundary_change: Totals::default(),
        right_boundary_change: Totals::default(),
        area_momentum_impulse: 0.0,
        gravity_momentum_impulse: 0.0,
        gravity_work: 0.0,
        conservation_residual: Totals::default(),
    }
}

#[cfg(test)]
#[path = "network/physics_tests.rs"]
mod physics_tests;
#[cfg(test)]
mod tests;
