//! The planar nonlinear winding owner as one Component of the composition
//! contract (fs-ports).
//!
//! [`PlanarWinding`] steps [`NonlinearWindingCircuit`] once per proposal, from
//! the accepted state, with conductors fixed at construction
//! ([`PlanarWinding::new`]) or given by laws at the temperatures its heat
//! ports receive ([`PlanarWinding::at_port_temperatures`]). Ports, in order:
//! `source<p>` for each original external input `p` of the native graph
//! (`EffortIn` for a voltage source, `FlowIn` for a current source),
//! returning the native step-average conjugate `StepRecord::y[p]`; then
//! `heat<k>` for each declared [`HeatReceiver`], receiving temperature and
//! returning minus the mean heat rate of the stranded cells, core cells and
//! native graph resistors mapped to it. Every dissipating element maps to
//! exactly one receiver. The owner is backward Euler with held inputs, so
//! every port answers [`EndOfStep::Untold`].
//!
//! The receipt takes the owner's stored change, its numerical (backward
//! Euler) dissipation and the heat that leaves through the heat ports. A
//! step with magnetic parameter work or mechanical work refuses: the receipt
//! has no term for either. The temperatures the heat ports receive enter
//! the generated entropy and, with laws, the conductors and graph resistors
//! of the step: a [`ConductorLaw`] per receiver for its triangles, and a
//! [`ResistorLaw`] per graph resistor through the owner's
//! `with_graph_resistors`, or none, keeping the circuit's own resistors. A
//! law with no answer refuses the step; the Component never clamps a
//! temperature.
//!
//! The Component borrows the owner's problem for `'m` and leaks nothing. A
//! system of the kernel holds `'static` members, so a caller that adds one
//! to a system gives it a problem that lives for `'static`. The state has no
//! codec: a system holding it cannot be written to a durable snapshot.

use core::fmt;
use core::marker::PhantomData;
use core::num::NonZeroUsize;
use std::collections::BTreeMap;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::conductor::{ResolvedConductor, ResolvedConductorState};
use fs_phs::arithmetic_rounding_floor;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::component::{
    Class, Component, ContractViolation, EndOfStep, PortLayout, PortRole, PortSpec, Proposal,
    Receipt,
};
use fs_ports::{CoupleError, PortKind, StableId};

use super::nonlinear::{
    NonlinearWindingCircuit, NonlinearWindingControls, NonlinearWindingError, NonlinearWindingStep,
    resistors::GraphResistorMaterial,
};
use super::ohmic::WindingConductor;
use crate::planar::eddy::EddyConductor;

/// One thermal receiver: the dissipating elements whose heat leaves through
/// one heat port.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct HeatReceiver {
    /// Mesh triangles of stranded winding or core conductors.
    pub triangles: Vec<usize>,
    /// Branch indices of native graph resistors.
    pub resistors: Vec<usize>,
}

/// The conductors of one receiver's triangles at the temperature its heat
/// port received.
pub trait ConductorLaw: fmt::Debug + Send + Sync {
    /// The source-bound state of `triangle`'s conductor at `kelvin`. The
    /// owner refuses a state whose card, basis or claims are not those of the
    /// circuit's own conductor there.
    ///
    /// # Errors
    /// Why the law has no answer there, for example a temperature outside the
    /// validity of the card behind it. The step then refuses with
    /// [`PlanarWindingError::Law`]; the temperature is never clamped.
    fn conductor(&self, triangle: usize, kelvin: f64) -> Result<ResolvedConductorState, String>;

    /// Words that identify the law. They join the Component's identity.
    fn identity_words(&self, words: &mut Vec<u64>);
}

/// One native graph resistor at the temperature its heat port received.
pub trait ResistorLaw: fmt::Debug + Send + Sync {
    /// The source-bound resistor at `kelvin`, a static Ohmic scalar as the
    /// owner's `with_graph_resistors` requires.
    ///
    /// # Errors
    /// As [`ConductorLaw::conductor`].
    fn conductor(&self, kelvin: f64) -> Result<ResolvedConductor, String>;

    /// Words that identify the law. They join the Component's identity.
    fn identity_words(&self, words: &mut Vec<u64>);
}

/// The element whose law had no answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LawElement {
    /// A stranded or core conductor's mesh triangle.
    Triangle(usize),
    /// A native graph resistor's branch.
    Resistor(usize),
}

/// Refusal of the winding owner or of its Component boundary.
#[derive(Debug)]
pub enum PlanarWindingError {
    /// The owner's refusal, unchanged.
    Owner(NonlinearWindingError),
    /// Invalid or repeated port metadata.
    Port(ContractViolation),
    /// An invalid port identifier.
    Identifier(CoupleError),
    /// A shape, a value or a construction that is refused.
    Invalid {
        /// What was refused.
        what: &'static str,
    },
    /// The step did magnetic parameter work or mechanical work, which the
    /// receipt has no term for.
    UnaccountedWork {
        /// Signed magnetic parameter work \[J\].
        parameter_j: f64,
        /// Signed mechanical work \[J\].
        mechanical_j: f64,
    },
    /// Supplied work, stored change, numerical dissipation and the heat that
    /// leaves do not close within the rounding floor and the declared
    /// tolerance.
    Balance {
        /// Supplied work minus stored change minus heat out minus the
        /// owner's numerical dissipation \[J\].
        defect_j: f64,
        /// The allowance \[J\].
        tolerance_j: f64,
    },
    /// A law had no answer at the temperature its port received.
    Law {
        /// The heat receiver whose port gave the temperature.
        receiver: usize,
        /// The element asked for.
        element: LawElement,
        /// The temperature \[K\].
        kelvin: f64,
        /// The law's reason.
        message: String,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl fmt::Display for PlanarWindingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Owner(error) => write!(f, "winding owner refused: {error}"),
            Self::Port(violation) => write!(f, "winding port refused: {violation}"),
            Self::Identifier(error) => write!(f, "winding port identifier refused: {error:?}"),
            Self::Invalid { what } => write!(f, "winding component refused: {what}"),
            Self::UnaccountedWork {
                parameter_j,
                mechanical_j,
            } => write!(
                f,
                "winding step did parameter work {parameter_j} J and mechanical work \
                 {mechanical_j} J, which the receipt has no term for"
            ),
            Self::Balance {
                defect_j,
                tolerance_j,
            } => write!(
                f,
                "winding energy account does not close: defect {defect_j} J beyond {tolerance_j} J"
            ),
            Self::Law {
                receiver,
                element,
                kelvin,
                message,
            } => write!(
                f,
                "winding law of {element:?} (receiver {receiver}) has no answer at {kelvin} K: \
                 {message}"
            ),
            Self::Cancelled => f.write_str("winding component cancelled"),
        }
    }
}

impl std::error::Error for PlanarWindingError {}

impl From<NonlinearWindingError> for PlanarWindingError {
    fn from(error: NonlinearWindingError) -> Self {
        Self::Owner(error)
    }
}

const fn invalid(what: &'static str) -> PlanarWindingError {
    PlanarWindingError::Invalid { what }
}

/// The accepted state: the owner's circuit, which carries the field history.
#[derive(Clone)]
pub struct PlanarWindingState<'m>(NonlinearWindingCircuit<'m>);

impl<'m> PlanarWindingState<'m> {
    /// The state of an accepted circuit.
    #[must_use]
    pub const fn new(circuit: NonlinearWindingCircuit<'m>) -> Self {
        Self(circuit)
    }

    /// The accepted circuit.
    #[must_use]
    pub const fn circuit(&self) -> &NonlinearWindingCircuit<'m> {
        &self.0
    }

    /// The accepted circuit, by value.
    #[must_use]
    pub fn into_circuit(self) -> NonlinearWindingCircuit<'m> {
        self.0
    }
}

fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

/// Equal states share the owner's problem (the same allocation) and graph,
/// and hold the same bits in every accepted coordinate: clock, circuit
/// state, currents, linkages, resistances, conductors and field potential,
/// energy and clock. The graph is compared by value: a step with resistor
/// laws builds a new one.
impl PartialEq for PlanarWindingState<'_> {
    fn eq(&self, other: &Self) -> bool {
        let (a, b) = (&self.0, &other.0);
        a.graph() == b.graph()
            && core::ptr::eq(a.field().problem(), b.field().problem())
            && a.time_s().to_bits() == b.time_s().to_bits()
            && same_bits(a.state(), b.state())
            && same_bits(a.currents_a(), b.currents_a())
            && same_bits(a.flux_linkages_wb(), b.flux_linkages_wb())
            && same_bits(a.resistance_ohm(), b.resistance_ohm())
            && a.conductors() == b.conductors()
            && a.field().conductors() == b.field().conductors()
            && same_bits(a.field().potential_wb_m(), b.field().potential_wb_m())
            && a.field().energy_j().to_bits() == b.field().energy_j().to_bits()
            && a.field().time_s().to_bits() == b.field().time_s().to_bits()
    }
}

impl fmt::Debug for PlanarWindingState<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PlanarWindingState")
            .field("time_s", &self.0.time_s())
            .field("currents_a", &self.0.currents_a())
            .field("state", &self.0.state())
            .finish_non_exhaustive()
    }
}

/// A step's circuit, when resistor laws rebuilt it, and its core and
/// stranded conductors.
type Supplied<'m> = (
    Option<NonlinearWindingCircuit<'m>>,
    Vec<EddyConductor>,
    Vec<WindingConductor>,
);

/// Where a step's conductors and graph resistances come from.
enum Supply {
    /// Given at construction and used at every step.
    Fixed {
        core: Vec<EddyConductor>,
        windings: Vec<WindingConductor>,
    },
    /// Laws asked at the temperatures the heat ports receive.
    Laws {
        /// Every conductor triangle, with its core region or `None` when
        /// stranded, as the construction's circuit holds them.
        cells: BTreeMap<usize, Option<u32>>,
        /// One per heat receiver.
        conductors: Vec<Arc<dyn ConductorLaw>>,
        /// One per native graph resistor, by branch, or none: the
        /// circuit's own resistors stay.
        resistors: BTreeMap<usize, Arc<dyn ResistorLaw>>,
    },
}

/// The planar nonlinear winding owner behind the composition contract's
/// source and heat ports. See the module documentation.
pub struct PlanarWinding<'m> {
    graph: CircuitGraph,
    sources: usize,
    supply: Supply,
    receivers: Vec<HeatReceiver>,
    triangle_port: BTreeMap<usize, usize>,
    resistor_port: BTreeMap<usize, usize>,
    controls: NonlinearWindingControls,
    energy_tolerance_j: f64,
    ports: Vec<PortSpec>,
    problem: PhantomData<&'m ()>,
}

impl fmt::Debug for PlanarWinding<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PlanarWinding")
            .field("sources", &self.sources)
            .field("receivers", &self.receivers)
            .field("energy_tolerance_j", &self.energy_tolerance_j)
            .finish_non_exhaustive()
    }
}

fn place(
    map: &mut BTreeMap<usize, usize>,
    element: usize,
    port: usize,
) -> Result<(), PlanarWindingError> {
    if map.insert(element, port).is_some() {
        return Err(invalid("a dissipating element is mapped to two receivers"));
    }
    Ok(())
}

impl<'m> PlanarWinding<'m> {
    /// The Component of `circuit`'s topology, stepped with these conductors
    /// and controls, its heat leaving through `receivers`.
    ///
    /// # Errors
    /// The energy tolerance is not positive and finite; an external port
    /// drives no source or two; a receiver is empty; a stranded or core
    /// triangle or a graph resistor is mapped to no receiver or to two; a
    /// receiver names an element that does not dissipate.
    pub fn new(
        circuit: &NonlinearWindingCircuit<'m>,
        core_conductors: Vec<EddyConductor>,
        winding_conductors: Vec<WindingConductor>,
        receivers: Vec<HeatReceiver>,
        controls: NonlinearWindingControls,
        energy_tolerance_j: f64,
    ) -> Result<Self, PlanarWindingError> {
        let triangles: Vec<usize> = core_conductors
            .iter()
            .map(|c| c.triangle)
            .chain(winding_conductors.iter().map(|c| c.triangle))
            .collect();
        let supply = Supply::Fixed {
            core: core_conductors,
            windings: winding_conductors,
        };
        Self::build(
            circuit,
            &triangles,
            supply,
            receivers,
            controls,
            energy_tolerance_j,
        )
    }

    /// The Component of `circuit`'s topology whose conductors and graph
    /// resistors follow the temperatures its heat ports receive: receiver
    /// `k`'s law gives the conductor of each of its triangles at port `k`'s
    /// temperature, and each graph resistor's law gives it at the
    /// temperature of the port its heat leaves through. With no resistor
    /// laws, the graph resistors keep `circuit`'s values. The core and
    /// stranded triangles, and the core regions, are `circuit`'s own.
    ///
    /// # Errors
    /// As [`Self::new`]; the number of conductor laws is not the number of
    /// receivers; some graph resistors have laws and others do not, or a law
    /// names a branch that is not a resistor.
    pub fn at_port_temperatures(
        circuit: &NonlinearWindingCircuit<'m>,
        conductor_laws: Vec<Arc<dyn ConductorLaw>>,
        resistor_laws: BTreeMap<usize, Arc<dyn ResistorLaw>>,
        receivers: Vec<HeatReceiver>,
        controls: NonlinearWindingControls,
        energy_tolerance_j: f64,
    ) -> Result<Self, PlanarWindingError> {
        if conductor_laws.len() != receivers.len() {
            return Err(invalid("one conductor law per heat receiver required"));
        }
        let cells: BTreeMap<usize, Option<u32>> = circuit
            .field()
            .conductors()
            .iter()
            .map(|c| (c.triangle, Some(c.region)))
            .chain(circuit.conductors().iter().map(|c| (c.triangle, None)))
            .collect();
        let triangles: Vec<usize> = cells.keys().copied().collect();
        let supply = Supply::Laws {
            cells,
            conductors: conductor_laws,
            resistors: resistor_laws,
        };
        Self::build(
            circuit,
            &triangles,
            supply,
            receivers,
            controls,
            energy_tolerance_j,
        )
    }

    fn build(
        circuit: &NonlinearWindingCircuit<'m>,
        triangles: &[usize],
        supply: Supply,
        receivers: Vec<HeatReceiver>,
        controls: NonlinearWindingControls,
        energy_tolerance_j: f64,
    ) -> Result<Self, PlanarWindingError> {
        if !(energy_tolerance_j.is_finite() && energy_tolerance_j > 0.0) {
            return Err(invalid("energy tolerance must be positive and finite"));
        }
        let graph = circuit.graph().clone();
        let mut roles: Vec<Option<PortRole>> = Vec::new();
        let mut resistors = Vec::new();
        for (branch, &(_, _, element)) in graph.branches.iter().enumerate() {
            let (port, role) = match element {
                Branch::Resistor { .. } => {
                    resistors.push(branch);
                    continue;
                }
                Branch::VoltageSource { port } => (port, PortRole::EffortIn),
                Branch::CurrentSource { port } => (port, PortRole::FlowIn),
                Branch::Inductor { .. } | Branch::Capacitor { .. } => continue,
            };
            if roles.len() <= port {
                roles.resize(port + 1, None);
            }
            if roles[port].replace(role).is_some() {
                return Err(invalid("an external port drives more than one source"));
            }
        }
        let mut triangle_port = BTreeMap::new();
        let mut resistor_port = BTreeMap::new();
        for (port, receiver) in receivers.iter().enumerate() {
            if receiver.triangles.is_empty() && receiver.resistors.is_empty() {
                return Err(invalid("a heat receiver maps no element"));
            }
            for &triangle in &receiver.triangles {
                if !triangles.contains(&triangle) {
                    return Err(invalid(
                        "a receiver names a triangle that carries no conductor",
                    ));
                }
                place(&mut triangle_port, triangle, port)?;
            }
            for &branch in &receiver.resistors {
                if !resistors.contains(&branch) {
                    return Err(invalid("a receiver names a branch that is not a resistor"));
                }
                place(&mut resistor_port, branch, port)?;
            }
        }
        if triangles
            .iter()
            .any(|triangle| !triangle_port.contains_key(triangle))
            || resistors
                .iter()
                .any(|branch| !resistor_port.contains_key(branch))
        {
            return Err(invalid("a dissipating element is mapped to no receiver"));
        }
        if let Supply::Laws {
            resistors: laws, ..
        } = &supply
        {
            if laws.keys().any(|branch| !resistors.contains(branch)) {
                return Err(invalid(
                    "a resistor law names a branch that is not a resistor",
                ));
            }
            // The owner rebinds every graph resistor or none.
            if !laws.is_empty() && resistors.iter().any(|branch| !laws.contains_key(branch)) {
                return Err(invalid("every graph resistor or none has a law"));
            }
        }
        let mut ports = Vec::with_capacity(roles.len() + receivers.len());
        for (port, role) in roles.iter().enumerate() {
            let role = role.ok_or_else(|| invalid("an external port drives no source"))?;
            let id =
                StableId::new(format!("source{port}")).map_err(PlanarWindingError::Identifier)?;
            ports.push(
                PortSpec::power(
                    id,
                    PortKind::ElectricalVoltageCurrent,
                    role,
                    NonZeroUsize::MIN,
                )
                .map_err(PlanarWindingError::Port)?,
            );
        }
        for k in 0..receivers.len() {
            let id = StableId::new(format!("heat{k}")).map_err(PlanarWindingError::Identifier)?;
            ports.push(PortSpec::heat(id, PortRole::EffortIn).map_err(PlanarWindingError::Port)?);
        }
        PortLayout::new(&ports).map_err(PlanarWindingError::Port)?;
        Ok(Self {
            graph,
            sources: roles.len(),
            supply,
            receivers,
            triangle_port,
            resistor_port,
            controls,
            energy_tolerance_j,
            ports,
            problem: PhantomData,
        })
    }

    /// The owner's clock at the end of a step of `dt_s` from `old`, when the
    /// step, its ports and the state's topology are this Component's.
    fn admit(
        &self,
        old: &PlanarWindingState<'m>,
        dt_s: f64,
        inputs: &[f64],
        outputs: usize,
    ) -> Result<f64, PlanarWindingError> {
        let ports = self.sources + self.receivers.len();
        if !(dt_s.is_finite() && dt_s > 0.0)
            || inputs.len() != ports
            || outputs != ports
            || inputs.iter().any(|value| !value.is_finite())
            || inputs[self.sources..].iter().any(|&kelvin| kelvin <= 0.0)
        {
            return Err(invalid(
                "positive duration, complete finite ports and positive temperatures required",
            ));
        }
        let same = match &self.supply {
            Supply::Laws { resistors, .. } if !resistors.is_empty() => {
                same_topology(old.0.graph(), &self.graph)
            }
            Supply::Fixed { .. } | Supply::Laws { .. } => old.0.graph() == &self.graph,
        };
        if !same {
            return Err(invalid(
                "the state's circuit is not this Component's topology",
            ));
        }
        let target_s = old.0.time_s() + dt_s;
        if !(target_s.is_finite() && target_s > old.0.time_s()) {
            return Err(invalid("the step does not advance the owner's clock"));
        }
        Ok(target_s)
    }

    /// The circuit to step, when resistor laws rebuild it, and the core and
    /// stranded conductors of a step whose heat ports received `kelvin`.
    fn supply_at(
        &self,
        old: &NonlinearWindingCircuit<'m>,
        kelvin: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Supplied<'m>, PlanarWindingError> {
        let (cells, conductors, resistors) = match &self.supply {
            Supply::Fixed { core, windings } => return Ok((None, core.clone(), windings.clone())),
            Supply::Laws {
                cells,
                conductors,
                resistors,
            } => (cells, conductors, resistors),
        };
        let mut core = Vec::new();
        let mut windings = Vec::new();
        for (&triangle, &region) in cells {
            let receiver = self.triangle_port[&triangle];
            let kelvin = kelvin[receiver];
            let conductor =
                conductors[receiver]
                    .conductor(triangle, kelvin)
                    .map_err(|message| PlanarWindingError::Law {
                        receiver,
                        element: LawElement::Triangle(triangle),
                        kelvin,
                        message,
                    })?;
            match region {
                Some(region) => core.push(EddyConductor {
                    triangle,
                    region,
                    conductor,
                }),
                None => windings.push(WindingConductor {
                    triangle,
                    conductor,
                }),
            }
        }
        if resistors.is_empty() {
            return Ok((None, core, windings));
        }
        let mut sources = Vec::with_capacity(resistors.len());
        for (&branch, law) in resistors {
            let receiver = self.resistor_port[&branch];
            let kelvin = kelvin[receiver];
            let conductor = law
                .conductor(kelvin)
                .map_err(|message| PlanarWindingError::Law {
                    receiver,
                    element: LawElement::Resistor(branch),
                    kelvin,
                    message,
                })?;
            sources.push(GraphResistorMaterial { branch, conductor });
        }
        let driven = old.with_graph_resistors(sources, cancelled)?;
        Ok((Some(driven), core, windings))
    }

    /// The heat of each receiver over a step [J], from the owner's records
    /// of its stranded cells, core cells and graph resistors, in that order.
    fn heat_by_receiver(
        &self,
        step: &NonlinearWindingStep,
    ) -> Result<Vec<f64>, PlanarWindingError> {
        let mut heat = vec![0.0; self.receivers.len()];
        let triangle = |t: usize| {
            self.triangle_port
                .get(&t)
                .copied()
                .ok_or_else(|| invalid("the owner reported heat in a triangle no receiver maps"))
        };
        for cell in &step.cells {
            heat[triangle(cell.source.triangle)?] += cell.heat_j;
        }
        for cell in &step.field.cells {
            heat[triangle(cell.source.triangle)?] += cell.heat_j;
        }
        for resistor in &step.graph_resistors {
            let port = self
                .resistor_port
                .get(&resistor.branch)
                .copied()
                .ok_or_else(|| invalid("the owner reported heat in a resistor no receiver maps"))?;
            heat[port] += resistor.heat_j;
        }
        if heat.iter().any(|q| !(q.is_finite() && *q >= 0.0)) {
            return Err(invalid("a receiver's heat is not finite and nonnegative"));
        }
        Ok(heat)
    }

    /// The declared heat receivers, in heat-port order.
    #[must_use]
    pub fn receivers(&self) -> &[HeatReceiver] {
        &self.receivers
    }
}

/// Circuit, field and end-mode storage of an accepted state, summed as the
/// owner sums `stored_before_j`.
fn energy(state: &PlanarWindingState<'_>) -> Result<f64, PlanarWindingError> {
    let circuit = &state.0;
    let end = circuit
        .end_inductance()
        .stored_energy_j(circuit.currents_a())
        .map_err(NonlinearWindingError::from)?;
    let total =
        circuit.dae().system.hamiltonian(circuit.state()) + circuit.field().energy_j() + end;
    if !total.is_finite() {
        return Err(invalid("the winding's stored energy is not finite"));
    }
    Ok(total)
}

fn conductor_words(
    conductor: &fs_material::conductor::ResolvedConductorState,
    words: &mut Vec<u64>,
) {
    words.push(conductor.resistivity_ohm_m().to_bits());
    for chunk in conductor.material().card_identity().0.chunks(8) {
        let mut bytes = [0_u8; 8];
        bytes.copy_from_slice(chunk);
        words.push(u64::from_le_bytes(bytes));
    }
}

/// A law's words, preceded by their count.
fn law_words(words: &mut Vec<u64>, law: impl FnOnce(&mut Vec<u64>)) {
    let at = words.len();
    words.push(0);
    law(words);
    words[at] = (words.len() - at - 1) as u64;
}

/// Equal graphs but for the values of their resistors, which resistor laws
/// set at every step.
fn same_topology(a: &CircuitGraph, b: &CircuitGraph) -> bool {
    a.node_count == b.node_count
        && a.transformers == b.transformers
        && a.branches.len() == b.branches.len()
        && a.branches.iter().zip(&b.branches).all(|(x, y)| {
            x.0 == y.0
                && x.1 == y.1
                && (matches!(
                    (x.2, y.2),
                    (Branch::Resistor { .. }, Branch::Resistor { .. })
                ) || x.2 == y.2)
        })
}

fn branch_words(element: Branch, words: &mut Vec<u64>) {
    match element {
        Branch::Resistor { ohms } => words.extend([0, ohms.to_bits()]),
        Branch::Inductor { henries } => words.extend([1, henries.to_bits()]),
        Branch::Capacitor { farads } => words.extend([2, farads.to_bits()]),
        Branch::VoltageSource { port } => words.extend([3, port as u64]),
        Branch::CurrentSource { port } => words.extend([4, port as u64]),
    }
}

impl<'m> Component for PlanarWinding<'m> {
    type State = PlanarWindingState<'m>;
    type Error = PlanarWindingError;

    /// The graph; the fixed conductors' triangles, regions, resistivity bits
    /// and card identities, or the laws' cells and words; the receiver map;
    /// the controls by their `Debug` text (which prints each `f64` so that it
    /// parses back to the same bits) and the energy tolerance. Each list is
    /// preceded by its length.
    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            self.graph.node_count as u64,
            self.graph.branches.len() as u64,
        ]);
        for &(a, b, element) in &self.graph.branches {
            words.extend([a as u64, b as u64]);
            branch_words(element, words);
        }
        words.push(self.graph.transformers.len() as u64);
        for link in &self.graph.transformers {
            words.extend([
                link.primary.0 as u64,
                link.primary.1 as u64,
                link.secondary.0 as u64,
                link.secondary.1 as u64,
                link.ratio.to_bits(),
            ]);
        }
        match &self.supply {
            Supply::Fixed { core, windings } => {
                words.extend([0, core.len() as u64]);
                for core in core {
                    words.extend([core.triangle as u64, u64::from(core.region)]);
                    conductor_words(&core.conductor, words);
                }
                words.push(windings.len() as u64);
                for winding in windings {
                    words.push(winding.triangle as u64);
                    conductor_words(&winding.conductor, words);
                }
            }
            Supply::Laws {
                cells,
                conductors,
                resistors,
            } => {
                words.extend([1, cells.len() as u64]);
                for (&triangle, &region) in cells {
                    words.extend([triangle as u64, region.map_or(u64::MAX, u64::from)]);
                }
                words.push(conductors.len() as u64);
                for law in conductors {
                    law_words(words, |words| law.identity_words(words));
                }
                words.push(resistors.len() as u64);
                for (&branch, law) in resistors {
                    words.push(branch as u64);
                    law_words(words, |words| law.identity_words(words));
                }
            }
        }
        words.push(self.receivers.len() as u64);
        for receiver in &self.receivers {
            words.push(receiver.triangles.len() as u64);
            words.extend(receiver.triangles.iter().map(|&t| t as u64));
            words.push(receiver.resistors.len() as u64);
            words.extend(receiver.resistors.iter().map(|&b| b as u64));
        }
        let controls = format!("{:?}", self.controls);
        words.push(controls.len() as u64);
        for chunk in controls.as_bytes().chunks(8) {
            let mut bytes = [0_u8; 8];
            bytes[..chunk.len()].copy_from_slice(chunk);
            words.push(u64::from_le_bytes(bytes));
        }
        words.push(self.energy_tolerance_j.to_bits());
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        // Backward Euler with held inputs: no flow ends its step by the
        // midpoint relation, and none is told.
        EndOfStep::Untold
    }

    fn stored_energy_j(&self, state: &PlanarWindingState<'m>) -> Result<f64, PlanarWindingError> {
        energy(state)
    }

    fn propose(
        &self,
        old: &PlanarWindingState<'m>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut PlanarWindingState<'m>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, PlanarWindingError> {
        cx.checkpoint().map_err(|_| PlanarWindingError::Cancelled)?;
        let target_s = self.admit(old, dt_s, inputs, outputs.len())?;
        let (sources, heats) = (self.sources, self.receivers.len());
        let (driven, core, windings) =
            self.supply_at(&old.0, &inputs[sources..], &mut || cx.checkpoint().is_err())?;
        let candidate = driven.as_ref().unwrap_or(&old.0).step_to_candidate(
            target_s,
            core,
            windings,
            &inputs[..sources],
            self.controls,
            &mut || cx.checkpoint().is_err(),
        )?;
        let step = candidate.step();
        if step.parameter_work_j != 0.0 || step.mechanical_work_j != 0.0 {
            return Err(PlanarWindingError::UnaccountedWork {
                parameter_j: step.parameter_work_j,
                mechanical_j: step.mechanical_work_j,
            });
        }
        let heat = self.heat_by_receiver(step)?;
        let heat_j: f64 = heat.iter().sum();
        let entropy: f64 = heat
            .iter()
            .zip(&inputs[sources..])
            .map(|(q, t)| q / t)
            .sum();
        let numerical_loss_j = step.supplied_work_j - step.stored_change_j - heat_j;
        let magnitude = step.stored_before_j.abs()
            + step.stored_after_j.abs()
            + step.supplied_work_j.abs()
            + step.numerical_dissipation_j.abs()
            + heat_j;
        let elements = step.cells.len() + step.field.cells.len() + step.graph_resistors.len();
        let rounding = arithmetic_rounding_floor(
            64_usize
                .saturating_add(8_usize.saturating_mul(old.0.state().len()))
                .saturating_add(8_usize.saturating_mul(old.0.field().potential_wb_m().len()))
                .saturating_add(8_usize.saturating_mul(elements + sources + heats)),
            magnitude,
        );
        let allowance = rounding + self.energy_tolerance_j;
        let defect_j = numerical_loss_j - step.numerical_dissipation_j;
        if ![rounding, numerical_loss_j, entropy, defect_j]
            .iter()
            .all(|value| value.is_finite())
        {
            return Err(invalid("the winding's energy account is not finite"));
        }
        if defect_j.abs() > allowance {
            return Err(PlanarWindingError::Balance {
                defect_j,
                tolerance_j: allowance,
            });
        }
        let report = &step.report;
        let initial = report
            .history
            .first()
            .map_or(report.residual_norm, |iteration| iteration.residual_before);
        let tolerance = self
            .controls
            .newton
            .absolute_tolerance
            .max(self.controls.newton.relative_tolerance * initial);
        cx.checkpoint().map_err(|_| PlanarWindingError::Cancelled)?;
        outputs[..sources].copy_from_slice(&step.native.y[..sources]);
        for (rate, q) in outputs[sources..].iter_mut().zip(&heat) {
            *rate = -q / dt_s;
        }
        let receipt = Receipt {
            stored_change_j: step.stored_change_j,
            dissipated_j: heat_j,
            numerical_loss_j,
            entropy_generated_j_per_k: entropy,
            energy_tolerance_j: rounding,
            residual: report.residual_norm,
            tolerance,
        };
        *new = PlanarWindingState(candidate.into_parts().0);
        Ok(Proposal::Step(receipt))
    }
}
