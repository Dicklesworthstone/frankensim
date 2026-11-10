//! The enthalpy timeline as a composition-kernel Component.
//!
//! [`ThermalTimeline`] puts [`DepositionSimulation`] behind
//! `fs_ports::component::Component` (contract v0), so the kernel in
//! `fs-system` can deliver heat to its cells from any other component and
//! route what leaves it, with no code written for the combination. Nothing
//! is integrated here: a step is the owner's own
//! [`DepositionSimulation::advance_observed`], with the heat the kernel
//! delivers given to the owner as its `step_sources`.
//!
//! Ports, in this order:
//!
//! - `cell<id>`, one per declared cell, a heat port receiving flow: the heat
//!   rate into that cell [W], held over the step. It returns the mean of the
//!   cell's temperatures at the two ends of the step [K].
//! - `environment`, only if declared, a heat port receiving effort: it
//!   returns the mean rate [W] at which heat entered the timeline through
//!   its own authored boundaries (ambient films, reservoirs, authored
//!   sources), which the owner integrates itself. The timeline does not use
//!   the temperature the port receives; the kernel needs it positive for its
//!   own account of the wire.
//! - `temperature<id>`, one per cell named to
//!   [`ThermalTimeline::publishing`], a signal output: that cell's
//!   temperature at the start of the kernel's step [K]. It does not depend
//!   on the heat the step delivers, so a consumer reads a lagged
//!   temperature. A hand-written coupling that refreshes its feedback at
//!   every thermal substep is reproduced when the kernel takes one substep
//!   per step.
//!
//! Boundary heat is on that wire or it does not exist: without the port, a
//! timeline in which any active cell has a film coefficient, a reservoir or
//! an authored source is refused ([`TimelineRefusal::Boundary`]) at
//! [`ThermalTimeline::start`] and at every step, and so is a step in which
//! an event gives a cell such a boundary. The test is on the declared
//! boundary, not on the heat it happens to exchange.
//!
//! Interval and energy: the owner's clock is one `f64`. A step of `dt` from
//! time `t` ends at the representable `t + dt`, and the owner integrates
//! over its own substep durations, whose sum is not exactly `dt` once `t`
//! is large beside `dt`. The kernel's energy is kept, not its rate: each
//! port's `rate * dt` is handed to the owner substep by substep, in
//! proportion to the clock interval each substep covers, and whatever is
//! left goes into the substep that ends the step. A step that does not move
//! the owner's clock is refused ([`TimelineRefusal::Interval`]). A step of
//! exactly the distance [`Component::next_event_s`] reported ends on that
//! event, although `t + (event - t)` may round to its neighbour.
//!
//! Step size: [`Component::limit_dt_s`] returns the smaller of the proposal
//! and the owner's explicit bound for a substep from that state (which
//! includes its `max_dt_s`), read from the owner without taking a step. A
//! kernel that honours it normally takes one thermal substep per step, so
//! every other component advances on the thermal clock, as the hand-written
//! couplings of this timeline do; a longer step is the owner's several
//! substeps. [`Component::next_event_s`] is the time to the next scheduled
//! event.
//!
//! Receipt: the stored change is the owner's ([`EnergyAudit`]). Numerical
//! loss is what entered through the ports less what the owner stored: minus
//! the owner's residual, and the rounding of handing the energy over. The
//! environment heat is the owner's external heat less the port heat it was
//! handed, so beside port heat it carries that heat's rounding (about 1e-16
//! of it); it is exactly zero while no boundary exchanges.
//!
//! Implicit conduction: [`ThermalTimeline::implicit`] steps the same
//! timeline with the owner's backward-Euler or SDIRK2 solve on a prepared
//! [`ImplicitConductionGraph`]
//! (`advance_observed_to` / `advance_sdirk2_observed_to`), with the same
//! ports and the same hand-over of energy. No conduction bound limits the
//! step: [`Component::limit_dt_s`] is the smaller of the proposal and the
//! controls' `maximum_dt_s`. The graph can own heat sources and boundary
//! laws that the cells' boundaries do not show, so an implicit timeline
//! always has the `environment` port, and its heat is the owner's external
//! heat less the port heat whether or not a cell's boundary exchanges.
//!
//! When the owner's nonlinear solve rejects an interval (no convergence, an
//! entropy decrease, a singular or unconverged linear solve), the step is
//! not taken: the Component answers `Proposal::Shorten` with half the step,
//! and the kernel retries every component over the shorter interval, as the
//! hand-written couplings of this solve re-run their mechanics over a halved
//! interval. The owner's own halving is switched off for the step, so the
//! controls' `maximum_halvings` is not used and the kernel's limit on
//! shortenings applies; a step whose half would be under the controls'
//! `minimum_dt_s` returns the owner's refusal. The receipt carries the
//! solve's largest cell residual and the tolerance it met. One step
//! allocates the solver's workspace. With scheduled events the graph admits
//! only what leaves its bindings as they are: a boundary change is taken, a
//! birth, a contact or a separation is the owner's refusal.
//!
//! Class: `Smooth`. With the explicit rule the returned temperatures are as
//! differentiable in the heat rates as the laws are in enthalpy; a law with
//! a latent interval has kinks at the interval's ends, where a kernel that
//! differences the port map gets a one-sided slope. With the implicit solve
//! the map is smooth to the solve's tolerance only: its value moves by up to
//! that tolerance where the iteration count changes, and a step is refused
//! outright on one side of a rejection.
//!
//! The kernel keeps component states as `'static` values, so a timeline
//! added to a system must borrow laws that live that long.
//!
//! No claim: the entropy generated inside the timeline (conduction between
//! cells at different temperatures, exchange with its boundaries) is
//! reported as zero. That is a lower bound for the conduction; for the
//! boundaries it is one only if the environment wire is at the authored
//! ambient temperature, which is not checked. A birth after the start is
//! refused, because mass and its enthalpy would enter with no port to carry
//! them; the port map has no analytic tangent; the owner's adaptive
//! (step-doubling) and staged variants are not wrapped; the owner's own
//! accuracy statements are unchanged by the wrapping.
//!
//! [`EnergyAudit`]: crate::enthalpy::EnergyAudit

use core::fmt;
use core::num::NonZeroUsize;
use std::cell::Cell;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::rc::Rc;

use fs_exec::Cx;
use fs_la::tridiagonal::TridiagonalError;
use fs_ports::StableId;
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, Proposal, Receipt,
};

use crate::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionError, ImplicitConductionGraph,
    ImplicitConductionObservedError, ImplicitConductionWorkspace,
};
use crate::enthalpy::{
    Boundary, CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
    StepConfig,
};

/// Why a [`ThermalTimeline`] could not be built.
#[derive(Debug, Clone, PartialEq)]
pub enum TimelineBuildError {
    /// No cell port and no environment port: nothing could be connected.
    NoPorts,
    /// A cell was named twice.
    Duplicate {
        /// The repeated cell.
        cell: u64,
    },
    /// A port identifier was refused.
    Identifier(fs_ports::CoupleError),
    /// A port specification was refused.
    Port(ContractViolation),
    /// The implicit controls' longest interval is not a positive time.
    Interval {
        /// The value given [s].
        maximum_dt_s: f64,
    },
}

impl fmt::Display for TimelineBuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoPorts => write!(f, "a thermal timeline needs a cell or an environment port"),
            Self::Duplicate { cell } => write!(f, "cell {cell} is named twice"),
            Self::Identifier(error) => write!(f, "port identifier refused: {error:?}"),
            Self::Port(violation) => write!(f, "port refused: {violation}"),
            Self::Interval { maximum_dt_s } => {
                write!(f, "a longest thermal interval of {maximum_dt_s} s is not a positive time")
            }
        }
    }
}

impl std::error::Error for TimelineBuildError {}

/// Why a [`ThermalTimeline`] refused a call.
#[derive(Debug, Clone, PartialEq)]
pub enum TimelineRefusal {
    /// The owner refused; its error is unchanged.
    Thermal(DepositionError),
    /// The owner's implicit solve refused; its error is unchanged.
    Implicit(ImplicitConductionError),
    /// An event is scheduled at the timeline's own time. A conduction graph
    /// is prepared against the timeline as it stands, so the events a run
    /// starts with are applied before the graph is built.
    PendingEvent {
        /// The timeline's time [s].
        time_s: f64,
    },
    /// A cell with a port is not active in the timeline.
    InactiveCell {
        /// The cell.
        cell: u64,
    },
    /// A birth fell inside the step: mass and enthalpy would enter with no
    /// port to carry them.
    Birth {
        /// The first cell born.
        cell: u64,
        /// Energy all the births of the step carried [J].
        energy_j: f64,
    },
    /// A cell has a boundary that exchanges heat (a film coefficient, a
    /// reservoir or an authored source), and no `environment` port was
    /// declared to carry that heat.
    Boundary {
        /// The first such cell.
        cell: u64,
    },
    /// A heat rate that is not finite.
    Input {
        /// The value received.
        value: f64,
    },
    /// A step that does not move the timeline's clock.
    Interval {
        /// The step asked for [s].
        dt_s: f64,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl fmt::Display for TimelineRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Thermal(error) => write!(f, "thermal timeline refused: {error}"),
            Self::Implicit(error) => write!(f, "thermal timeline refused: {error}"),
            Self::PendingEvent { time_s } => write!(
                f,
                "an event is still scheduled at the timeline's time, {time_s} s; apply it before \
                 the conduction graph is prepared"
            ),
            Self::InactiveCell { cell } => write!(f, "cell {cell} has a port and is not active"),
            Self::Birth { cell, energy_j } => write!(
                f,
                "cell {cell} was born inside the step ({energy_j} J of births); births after \
                 the start need a stream port, which this component does not have"
            ),
            Self::Boundary { cell } => write!(
                f,
                "cell {cell} has a boundary that exchanges heat and no environment port was \
                 declared"
            ),
            Self::Input { value } => write!(f, "heat rate {value} is not finite"),
            Self::Interval { dt_s } => {
                write!(f, "a step of {dt_s} s does not move the timeline's clock")
            }
            Self::Cancelled => write!(f, "cancelled during the thermal step"),
        }
    }
}

impl std::error::Error for TimelineRefusal {}

impl From<DepositionError> for TimelineRefusal {
    fn from(error: DepositionError) -> Self {
        match error {
            DepositionError::Cancelled => Self::Cancelled,
            error => Self::Thermal(error),
        }
    }
}

impl From<ImplicitConductionError> for TimelineRefusal {
    fn from(error: ImplicitConductionError) -> Self {
        match error {
            ImplicitConductionError::Thermal(error) => error.into(),
            error => Self::Implicit(error),
        }
    }
}

/// Whether the owner's solve rejected the interval, as distinct from
/// refusing the problem: the errors on which the owner itself retries a
/// shorter interval.
fn rejected(error: &ImplicitConductionError) -> bool {
    matches!(
        error,
        ImplicitConductionError::NotConverged { .. }
            | ImplicitConductionError::EntropyDecrease { .. }
            | ImplicitConductionError::KrylovNotConverged { .. }
            | ImplicitConductionError::Linear(TridiagonalError::Singular { .. })
    )
}

/// The owner's implicit thermal quadrature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImplicitScheme {
    /// Backward Euler, `advance_observed_to`.
    BackwardEuler,
    /// The owner's two-stage SDIRK, `advance_sdirk2_observed_to`.
    Sdirk2,
}

/// How the owner advances the timeline over a step.
#[derive(Clone)]
enum Driver<'a> {
    /// The explicit rule, `DepositionSimulation::advance_observed`.
    Explicit,
    /// An implicit solve on a prepared conduction graph.
    Implicit {
        graph: Rc<ImplicitConductionGraph<'a>>,
        controls: ImplicitConductionControls,
        scheme: ImplicitScheme,
    },
}

impl fmt::Debug for Driver<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Explicit => write!(f, "Explicit"),
            Self::Implicit {
                graph,
                controls,
                scheme,
            } => f
                .debug_struct("Implicit")
                .field("cells", &graph.cell_ids().len())
                .field("controls", controls)
                .field("scheme", scheme)
                .finish(),
        }
    }
}

/// A timeline inside the kernel: the owner's own checkpoint, admitted by
/// [`ThermalTimeline::start`] or produced by a step.
///
/// Such a state has no event pending at its own time: the owner applies the
/// events at the end time of an advance.
///
/// Two states are equal when their times, active cells (with their
/// enthalpies and boundaries), links, next scheduled event time and budgets
/// are equal. Equality does not see the schedule past the next event or the
/// contact coverage of a link.
#[derive(Clone)]
pub struct TimelineState<'a> {
    simulation: DepositionSimulation<'a>,
}

impl<'a> TimelineState<'a> {
    /// The owner's timeline at this state.
    #[must_use]
    pub fn simulation(&self) -> &DepositionSimulation<'a> {
        &self.simulation
    }

    /// The owner's timeline, to carry on outside the kernel.
    #[must_use]
    pub fn into_simulation(self) -> DepositionSimulation<'a> {
        self.simulation
    }
}

impl PartialEq for TimelineState<'_> {
    fn eq(&self, other: &Self) -> bool {
        let (a, b) = (&self.simulation, &other.simulation);
        a.time_s().to_bits() == b.time_s().to_bits()
            && a.cells() == b.cells()
            && a.links() == b.links()
            && a.next_event_time_s().map(f64::to_bits) == b.next_event_time_s().map(f64::to_bits)
            && a.budget() == b.budget()
    }
}

impl fmt::Debug for TimelineState<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TimelineState")
            .field("time_s", &self.simulation.time_s())
            .field("cells", &self.simulation.cells())
            .field("links", &self.simulation.links().len())
            .finish()
    }
}

/// Whether a boundary can exchange heat with the cell it is on.
fn exchanges(boundary: &Boundary) -> bool {
    boundary.htc_w_m2_k != 0.0 || boundary.reservoir.is_some() || boundary.source_w != 0.0
}

/// What the wrapper tells the owner and reads back from it over one advance.
#[derive(Clone)]
struct Held {
    /// Energy still to hand to each port cell [J], and the time the advance
    /// ends [s].
    remaining_j: BTreeMap<CellId, f64>,
    target_s: f64,
    /// The owner's explicit bound as it was first proposed. Shared, so that
    /// it can be read after an advance that ends in an error.
    first_bound_s: Rc<Cell<Option<f64>>>,
    /// Stop at that bound, before any enthalpy changes.
    probe: bool,
    /// The shares prepared for the interval now on trial [J]. The owner may
    /// reject a trial and prepare a shorter one; a share is handed over only
    /// when the owner accepts its interval.
    pending_j: BTreeMap<CellId, f64>,
    /// Port heat handed to the owner [J], and the sum of its terms' sizes.
    took_j: f64,
    took_terms_j: f64,
    /// The first cell born inside the advance.
    born: Option<u64>,
    /// The first cell an event gave a boundary that exchanges heat.
    exchanging: Option<u64>,
}

impl Held {
    /// An observer that hands `energy_j` to the owner over an advance to
    /// `target_s`, or, as a `probe`, stops the owner at its first bound.
    fn new(energy_j: BTreeMap<CellId, f64>, target_s: f64, probe: bool) -> Self {
        Self {
            remaining_j: energy_j,
            target_s,
            first_bound_s: Rc::new(Cell::new(None)),
            probe,
            pending_j: BTreeMap::new(),
            took_j: 0.0,
            took_terms_j: 0.0,
            born: None,
            exchanging: None,
        }
    }
}

impl DepositionObserver for Held {
    type Error = DepositionError;

    fn requires_step_interval(&self) -> bool {
        true
    }

    /// The mean rate into each port cell over one trial interval: that
    /// interval's share of the energy still to hand over. The interval that
    /// ends the advance takes all that is left, so the cell receives the
    /// kernel's `rate * dt` whatever the owner's durations sum to. Nothing
    /// is booked here: the owner may reject the trial and ask again.
    fn step_sources(
        &mut self,
        _thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        // Bit for bit: the owner ends its last substep on the target itself.
        let last = end_s.to_bits() == self.target_s.to_bits();
        let fraction = duration_s / (self.target_s - start_s);
        let mut rates_w = BTreeMap::new();
        self.pending_j.clear();
        for (cell, remaining_j) in &self.remaining_j {
            let share_j = if last {
                *remaining_j
            } else {
                *remaining_j * fraction
            };
            self.pending_j.insert(*cell, share_j);
            rates_w.insert(*cell, share_j / duration_s);
        }
        Ok(rates_w)
    }

    fn limit_timestep_s(
        &mut self,
        _thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, DepositionError> {
        if self.first_bound_s.get().is_none() {
            self.first_bound_s.set(Some(proposed_s));
        }
        if self.probe {
            // The bound is all a probe is for. The owner discards the
            // candidate on any error; this one is recognised by the bound
            // having been recorded.
            return Err(DepositionError::Cancelled);
        }
        Ok(proposed_s)
    }

    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), DepositionError> {
        match event {
            Event::Birth(cell) => {
                self.born.get_or_insert(cell.id.0);
            }
            Event::Boundary(cell, boundary) if exchanges(boundary) => {
                self.exchanging.get_or_insert(cell.0);
            }
            Event::Boundary(..) | Event::Contact(_) | Event::Separation(_) => {}
        }
        Ok(())
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        // The owner accepted the interval: its shares are handed over.
        for (cell, share_j) in core::mem::take(&mut self.pending_j) {
            if let Some(remaining_j) = self.remaining_j.get_mut(&cell) {
                *remaining_j -= share_j;
            }
            self.took_j += share_j;
            self.took_terms_j += share_j.abs();
        }
        Ok(())
    }
}

/// The enthalpy timeline of [`DepositionSimulation`], stepped by its owner:
/// by the explicit rule ([`ThermalTimeline::new`]) or by an implicit solve
/// on a conduction graph ([`ThermalTimeline::implicit`]). `'a` is the
/// lifetime of the thermal laws the timeline borrows.
#[derive(Debug, Clone)]
pub struct ThermalTimeline<'a> {
    cells: Vec<CellId>,
    environment: bool,
    published: Vec<CellId>,
    config: StepConfig,
    driver: Driver<'a>,
    ports: Vec<PortSpec>,
    laws: PhantomData<&'a ()>,
}

impl<'a> ThermalTimeline<'a> {
    /// A timeline with a heat port for each of `cells`, in that order, and
    /// an `environment` port if `environment` is set, stepped with the
    /// owner's `config`.
    ///
    /// # Errors
    /// [`TimelineBuildError`] for no port at all or a cell named twice.
    pub fn new(
        cells: Vec<CellId>,
        environment: bool,
        config: StepConfig,
    ) -> Result<Self, TimelineBuildError> {
        if cells.is_empty() && !environment {
            return Err(TimelineBuildError::NoPorts);
        }
        let mut ports = Vec::with_capacity(cells.len() + usize::from(environment));
        for (index, cell) in cells.iter().enumerate() {
            if cells[..index].contains(cell) {
                return Err(TimelineBuildError::Duplicate { cell: cell.0 });
            }
            let id = StableId::new(format!("cell{}", cell.0))
                .map_err(TimelineBuildError::Identifier)?;
            ports.push(PortSpec::heat(id, PortRole::FlowIn).map_err(TimelineBuildError::Port)?);
        }
        if environment {
            let id = StableId::new("environment").map_err(TimelineBuildError::Identifier)?;
            ports.push(PortSpec::heat(id, PortRole::EffortIn).map_err(TimelineBuildError::Port)?);
        }
        Ok(Self {
            cells,
            environment,
            published: Vec::new(),
            config,
            driver: Driver::Explicit,
            ports,
            laws: PhantomData,
        })
    }

    /// A timeline stepped by the owner's implicit solve on `graph`, with a
    /// heat port for each of `cells`, in that order, and the `environment`
    /// port. `controls.maximum_dt_s` is the longest thermal interval and the
    /// step limit; `controls.total_tolerance_j` gates each step's balance.
    ///
    /// # Errors
    /// [`TimelineBuildError::Duplicate`] for a cell named twice;
    /// [`TimelineBuildError::Interval`] for a longest interval that is not
    /// positive and finite. The owner checks the other controls at a step.
    pub fn implicit(
        cells: Vec<CellId>,
        graph: ImplicitConductionGraph<'a>,
        controls: ImplicitConductionControls,
        scheme: ImplicitScheme,
    ) -> Result<Self, TimelineBuildError> {
        if !(controls.maximum_dt_s > 0.0 && controls.maximum_dt_s.is_finite()) {
            return Err(TimelineBuildError::Interval {
                maximum_dt_s: controls.maximum_dt_s,
            });
        }
        let config = StepConfig {
            max_dt_s: controls.maximum_dt_s,
            energy_tolerance_j: controls.total_tolerance_j,
        };
        let mut timeline = Self::new(cells, true, config)?;
        timeline.driver = Driver::Implicit {
            graph: Rc::new(graph),
            controls,
            scheme,
        };
        Ok(timeline)
    }

    /// Also publish the temperature of each of `cells` at the start of the
    /// step, as signal outputs `temperature<id>` [K] after the ports already
    /// declared. A published cell needs no heat port of its own; the
    /// timeline still needs a heat or an environment port.
    ///
    /// # Errors
    /// [`TimelineBuildError::Duplicate`] for a cell published twice.
    pub fn publishing(mut self, cells: Vec<CellId>) -> Result<Self, TimelineBuildError> {
        for cell in cells {
            if self.published.contains(&cell) {
                return Err(TimelineBuildError::Duplicate { cell: cell.0 });
            }
            let id = StableId::new(format!("temperature{}", cell.0))
                .map_err(TimelineBuildError::Identifier)?;
            self.ports.push(PortSpec::signal_out(id, NonZeroUsize::MIN));
            self.published.push(cell);
        }
        Ok(self)
    }

    /// The state to hand to the kernel: `simulation` with the events
    /// scheduled at its current time applied (the births a run starts
    /// with).
    ///
    /// An implicit timeline applies no events here: its graph was prepared
    /// against `simulation` as it stands.
    ///
    /// # Errors
    /// The owner's refusal; a cell with a port or a published temperature
    /// that is not active; without an `environment` port, a cell whose
    /// boundary exchanges heat; for an implicit timeline, an event still
    /// scheduled at the simulation's time.
    pub fn start(
        &self,
        simulation: &DepositionSimulation<'a>,
    ) -> Result<TimelineState<'a>, TimelineRefusal> {
        let time_s = simulation.time_s();
        let started = match self.driver {
            Driver::Explicit => simulation.advance_to(time_s, self.config, || Ok(()))?.0,
            Driver::Implicit { .. } => {
                if simulation
                    .next_event_time_s()
                    .is_some_and(|event_s| event_s <= time_s)
                {
                    return Err(TimelineRefusal::PendingEvent { time_s });
                }
                simulation.clone()
            }
        };
        if let Some(cell) = self
            .cells
            .iter()
            .chain(&self.published)
            .find(|cell| !started.cells().contains_key(cell))
        {
            return Err(TimelineRefusal::InactiveCell { cell: cell.0 });
        }
        if let (false, Some(cell)) = (self.environment, exchanging_cell(&started)) {
            return Err(TimelineRefusal::Boundary { cell });
        }
        Ok(TimelineState {
            simulation: started,
        })
    }
}

/// The first active cell whose boundary exchanges heat.
fn exchanging_cell(simulation: &DepositionSimulation<'_>) -> Option<u64> {
    simulation
        .cells()
        .values()
        .find(|cell| exchanges(&cell.boundary))
        .map(|cell| cell.id.0)
}

/// `sum |m e|` over the active cells: the size of the terms of the stored
/// energy, which may carry a large datum and cancel.
fn stored_terms_j(simulation: &DepositionSimulation<'_>) -> f64 {
    simulation
        .cells()
        .values()
        .map(|cell| (cell.mass_kg * cell.specific_h_j_kg).abs())
        .sum()
}

impl<'a> Component for ThermalTimeline<'a> {
    type State = TimelineState<'a>;
    type Error = TimelineRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &TimelineState<'a>) -> Result<f64, TimelineRefusal> {
        let energy_j = state.simulation.energy_j();
        if energy_j.is_finite() {
            Ok(energy_j)
        } else {
            Err(TimelineRefusal::Thermal(DepositionError::Arithmetic))
        }
    }

    fn limit_dt_s(
        &self,
        state: &TimelineState<'a>,
        proposed_s: f64,
    ) -> Result<f64, TimelineRefusal> {
        if let Driver::Implicit { controls, .. } = &self.driver {
            // No conduction bound: the longest interval the controls admit.
            return Ok(proposed_s.min(controls.maximum_dt_s));
        }
        // The owner proposes its explicit bound to its observer before it
        // changes any enthalpy; the probe stops it there. The bound comes
        // from the start state's capacities and conductances, not from the
        // heat delivered, so no rates are needed and no trial step is taken.
        let target_s = state.simulation.time_s() + proposed_s;
        let probe = Held::new(BTreeMap::new(), target_s, true);
        let outcome =
            state
                .simulation
                .advance_observed_limited(&probe, target_s, self.config, 1, || Ok(()));
        match (probe.first_bound_s.get(), outcome) {
            (Some(bound_s), _) => Ok(proposed_s.min(bound_s)),
            // No substep was needed: the proposal does not move the clock,
            // or nothing is active.
            (None, Ok(_)) => Ok(proposed_s),
            (None, Err(error)) => Err(error.into()),
        }
    }

    fn next_event_s(&self, state: &TimelineState<'a>) -> Option<f64> {
        // A state has no event pending at its own time (see
        // `TimelineState`), so this is the first event ahead of it.
        state
            .simulation
            .next_event_time_s()
            .map(|event_s| event_s - state.simulation.time_s())
            .filter(|ahead_s| *ahead_s > 0.0)
    }

    fn propose(
        &self,
        old: &TimelineState<'a>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut TimelineState<'a>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, TimelineRefusal> {
        cx.checkpoint().map_err(|_| TimelineRefusal::Cancelled)?;
        let exchanging = exchanging_cell(&old.simulation);
        if let (false, Some(cell)) = (self.environment, exchanging) {
            return Err(TimelineRefusal::Boundary { cell });
        }
        let start_s = old.simulation.time_s();
        let mut target_s = start_s + dt_s;
        // A step of exactly the reported distance to the next event ends on
        // it: `start + (event - start)` can round to a neighbour of the
        // event's time, which would leave a substep, or a step, one ulp long.
        if let Some(event_s) = old.simulation.next_event_time_s()
            && dt_s.to_bits() == (event_s - start_s).to_bits()
        {
            target_s = event_s;
        }
        if !(dt_s > 0.0 && target_s > start_s && target_s.is_finite()) {
            return Err(TimelineRefusal::Interval { dt_s });
        }
        // The energy the kernel books for each port, for the owner to store.
        let mut energy_j = BTreeMap::new();
        let (mut delivered_j, mut delivered_terms_j) = (0.0_f64, 0.0_f64);
        for (cell, &rate_w) in self.cells.iter().zip(inputs) {
            if !rate_w.is_finite() {
                return Err(TimelineRefusal::Input { value: rate_w });
            }
            energy_j.insert(*cell, dt_s * rate_w);
            delivered_j += dt_s * rate_w;
            delivered_terms_j += (dt_s * rate_w).abs();
        }
        let observer = Held::new(energy_j, target_s, false);
        // The implicit solve's largest cell residual and its tolerance.
        let mut solve = (0.0_f64, 0.0_f64);
        let (next, held, audit) = match &self.driver {
            Driver::Explicit => {
                old.simulation
                    .advance_observed(&observer, target_s, self.config, || {
                        cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                    })?
            }
            Driver::Implicit {
                graph,
                controls,
                scheme,
            } => {
                let mut workspace = ImplicitConductionWorkspace::new(graph.cell_ids().len())?;
                let checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
                // A rejected interval is the kernel's to shorten, with every
                // other component: the owner does not halve it here.
                let single = ImplicitConductionControls {
                    maximum_halvings: 0,
                    ..*controls
                };
                let advanced = match scheme {
                    ImplicitScheme::BackwardEuler => graph.advance_observed_to(
                        &old.simulation,
                        &observer,
                        target_s,
                        single,
                        &mut workspace,
                        checkpoint,
                    ),
                    ImplicitScheme::Sdirk2 => graph.advance_sdirk2_observed_to(
                        &old.simulation,
                        &observer,
                        target_s,
                        single,
                        &mut workspace,
                        checkpoint,
                    ),
                };
                let (next, held, report) = match advanced {
                    Ok(advanced) => advanced,
                    Err(ImplicitConductionObservedError::Implicit(error))
                        if rejected(&error) && 0.5 * dt_s >= controls.minimum_dt_s =>
                    {
                        return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
                    }
                    Err(ImplicitConductionObservedError::Implicit(error)) => {
                        return Err(error.into());
                    }
                    Err(ImplicitConductionObservedError::Observer(error)) => {
                        return Err(error.into());
                    }
                };
                solve = (report.maximum_cell_residual_j, controls.cell_tolerance_j);
                (next, held, report.energy)
            }
        };
        if let Some(cell) = held.born {
            return Err(TimelineRefusal::Birth {
                cell,
                energy_j: audit.incoming_energy_j,
            });
        }
        if let (false, Some(cell)) = (self.environment, held.exchanging) {
            return Err(TimelineRefusal::Boundary { cell });
        }

        for (port, cell) in self.cells.iter().enumerate() {
            let before_k = old.simulation.temperature_k(*cell)?;
            let after_k = next.temperature_k(*cell)?;
            outputs[port] = 0.5 * (before_k + after_k);
        }
        let signals = self.cells.len() + usize::from(self.environment);
        for (port, cell) in self.published.iter().enumerate() {
            outputs[signals + port] = old.simulation.temperature_k(*cell)?;
        }
        // The owner's external heat is the port heat it was handed plus what
        // its own boundaries exchanged. The second part is on the
        // environment wire; while no boundary exchanges it is nothing, not
        // the rounding between the two sums. A conduction graph can own
        // exchanges no cell's boundary shows, so there it is always read.
        let owned = matches!(self.driver, Driver::Implicit { .. });
        let environment_j = if owned || exchanging.is_some() || held.exchanging.is_some() {
            audit.external_heat_j - held.took_j
        } else {
            0.0
        };
        // What the kernel books as entering through the ports.
        let booked_j = if self.environment {
            outputs[self.cells.len()] = environment_j / dt_s;
            delivered_j + environment_j
        } else {
            delivered_j
        };

        // The tolerance covers the kernel's port sums, the re-measurement of
        // the stored energy at both ends (a plain sum of terms that carry
        // the laws' datum, whose rounding grows with the number of cells),
        // and the owner's rounded totals.
        let scale_j = delivered_terms_j
            + held.took_terms_j
            + environment_j.abs()
            + audit.external_heat_j.abs()
            + stored_terms_j(&old.simulation)
            + stored_terms_j(&next)
            + audit.stored_energy_change_j.abs()
            + audit.residual_j.abs();
        #[allow(clippy::cast_precision_loss)]
        let units = 64.0 + 2.0 * next.cells().len() as f64;
        new.simulation = next;
        Ok(Proposal::Step(Receipt {
            stored_change_j: audit.stored_energy_change_j,
            dissipated_j: 0.0,
            // What entered through the ports less what the owner stored:
            // minus the owner's residual, and the rounding of handing the
            // energy over.
            numerical_loss_j: booked_j - audit.stored_energy_change_j,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: units * f64::EPSILON * scale_j + f64::MIN_POSITIVE,
            residual: solve.0,
            tolerance: solve.1,
        }))
    }
}
