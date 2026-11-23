//! A Kirchhoff circuit of fs-phs on the composition kernel, with one heat
//! port per resistor.
//!
//! [`KirchhoffCircuit`] puts the descriptor port-Hamiltonian circuit of
//! `fs_phs::circuit` behind `fs_ports::component::Component`. Any netlist
//! that owner admits (resistors, inductors with mutual coupling, capacitors,
//! voltage and current sources, ideal transformers) is one Component. Its
//! resistors deliver their Joule heat to whatever the kernel joins their heat
//! ports to, and a resistor may follow the temperature its heat port receives
//! through a [`ResistanceLaw`]. Nothing is integrated here: a step is the
//! owner's own `CircuitDae::step_audited` on the circuit assembled with the
//! resistances of that step.
//!
//! Ports, in this order:
//!
//! - `source<p>`, one per external port `p` of the circuit, in port order. A
//!   port that drives a voltage source receives that voltage (effort, \[V\])
//!   and returns the current the source delivers into the circuit, the
//!   owner's step-mean conjugate output (flow, \[A\]). A port that drives a
//!   current source receives the current and returns the voltage. The power
//!   into the circuit is their product, the owner's `u^T y`.
//! - `heat<b>`, one per resistor, named by the resistor's index `b` in the
//!   graph's branch list, in branch order. It receives the temperature \[K\]
//!   at which the resistor's heat leaves and returns the heat rate into the
//!   resistor \[W\], which is minus its Joule power.
//!
//! State: the differential coordinates the owner keeps as history, the
//! inductors' flux linkages in the owner's flux order and then the
//! capacitors' charges in its charge order. Node potentials and constraint
//! currents are algebraic multipliers of one step and are not state: each
//! step starts from them at zero, as `fs_couple::material_circuit` does.
//!
//! Heat: the heat of resistor `b` over a step is `dt v i` with `v` the
//! voltage across it at the midpoint of the step (the mean of the step's two
//! ends of its terminals' potentials, the algebraic midpoint the owner's
//! dissipation uses) and `i = v / R`, `R` its resistance in that step. The
//! sum over resistors equals the owner's dissipation to rounding; the
//! difference, and the owner's own balance defect, are booked as numerical
//! loss and refused beyond the declared energy tolerance.
//!
//! Resistance: a resistor with no law keeps the resistance in the graph. A
//! resistor following a law is given `law.ohms(T)` at the temperature `T`
//! its heat port receives in that step. With a backward-Euler thermal
//! timeline that returns end-of-step temperatures, that is the endpoint
//! material resistance of the hand-written circuit/enthalpy couplings; with
//! a heat capacity that returns its step-mean temperature, it is the
//! midpoint's. A law with no answer refuses the step; it is never clamped.
//!
//! Class: [`Class::Smooth`]. The heat is quadratic in the sources and a law
//! makes the map nonlinear in temperature; no tangent is claimed, so a
//! kernel differences the port map.
//!
//! Cost: every proposal assembles the circuit, because a law can change any
//! resistance. The owner's admission (floating nodes, shorted or parallel
//! sources, bad values) therefore runs at every step as well as at
//! construction.

use core::fmt;
use core::num::NonZeroUsize;
use std::sync::Arc;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_phs::arithmetic_rounding_floor;
use fs_phs::circuit::{
    Branch, CircuitDae, CircuitError, CircuitGraph, MutualInductance,
    assemble_circuit_with_mutual_inductance,
};
use fs_ports::component::{
    Class, Component, ContractViolation, PortLayout, PortRole, PortSpec, Proposal, Receipt,
};
use fs_ports::{CoupleError, PortKind, StableId};

/// A resistance that follows the temperature where the resistor's heat
/// leaves.
pub trait ResistanceLaw: fmt::Debug + Send + Sync {
    /// The resistance \[ohm\] at `kelvin`.
    ///
    /// # Errors
    /// Why the law has no answer there, for example a temperature outside the
    /// validity of the card behind it. The circuit then refuses the step with
    /// [`CircuitComponentError::Law`]; it never clamps the temperature or
    /// extrapolates the law.
    fn ohms(&self, kelvin: f64) -> Result<f64, String>;

    /// Words that identify the law. They join the circuit's identity, so two
    /// circuits whose laws differ differ in identity.
    fn identity_words(&self, words: &mut Vec<u64>);
}

/// Refusal of the circuit owner or of its Component boundary.
#[derive(Debug)]
pub enum CircuitComponentError {
    /// The owner's refusal, unchanged: admission of the netlist or a step.
    Circuit(CircuitError),
    /// A resistance law had no answer, or answered a resistance that is not
    /// positive and finite.
    Law {
        /// Index of the resistor in the graph's branch list.
        branch: usize,
        /// The law's reason.
        reason: String,
    },
    /// Invalid or repeated port metadata.
    Port(ContractViolation),
    /// An invalid port identifier.
    Identifier(CoupleError),
    /// A shape, a value or a construction that is refused.
    Invalid {
        /// What was refused.
        what: &'static str,
    },
    /// Work, heat and storage do not close within the declared tolerance and
    /// rounding.
    Balance {
        /// Supplied work minus resistor heat minus stored change \[J\].
        defect_j: f64,
        /// The allowance \[J\].
        tolerance_j: f64,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl fmt::Display for CircuitComponentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Circuit(error) => write!(f, "circuit owner refused: {error:?}"),
            Self::Law { branch, reason } => {
                write!(
                    f,
                    "resistance law of branch {branch} has no answer: {reason}"
                )
            }
            Self::Port(error) => write!(f, "circuit port refused: {error:?}"),
            Self::Identifier(error) => write!(f, "circuit port identifier refused: {error:?}"),
            Self::Invalid { what } => write!(f, "circuit component refused: {what}"),
            Self::Balance {
                defect_j,
                tolerance_j,
            } => write!(
                f,
                "circuit balance defect {defect_j} J exceeds {tolerance_j} J"
            ),
            Self::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for CircuitComponentError {}

impl From<CircuitError> for CircuitComponentError {
    fn from(error: CircuitError) -> Self {
        Self::Circuit(error)
    }
}

fn invalid(what: &'static str) -> CircuitComponentError {
    CircuitComponentError::Invalid { what }
}

/// One resistor of the graph: its branch index, its terminals and, if it
/// follows one, its law.
#[derive(Debug, Clone)]
struct Resistor {
    branch: usize,
    a: usize,
    b: usize,
    law: Option<Arc<dyn ResistanceLaw>>,
}

/// A Kirchhoff circuit of `fs_phs::circuit` as one Component.
///
/// See the module documentation for the ports, the state, the heat and the
/// resistances.
pub struct KirchhoffCircuit {
    graph: CircuitGraph,
    mutual: Vec<MutualInductance>,
    resistors: Vec<Resistor>,
    /// The owner's assembly at the graph's own resistances: index maps and
    /// the stored energy, neither of which depends on a resistance.
    reference: CircuitDae,
    ports: Vec<PortSpec>,
    energy_tolerance_j: f64,
}

impl fmt::Debug for KirchhoffCircuit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("KirchhoffCircuit")
            .field("graph", &self.graph)
            .field("mutual", &self.mutual)
            .field("resistors", &self.resistors)
            .field("energy_tolerance_j", &self.energy_tolerance_j)
            .finish_non_exhaustive()
    }
}

impl KirchhoffCircuit {
    /// The circuit `graph` with the signed mutual inductances `mutual`
    /// between its inductor branches, every resistor at the resistance in
    /// the graph.
    ///
    /// `energy_tolerance_j` \[J\] is the largest balance defect a step may
    /// leave, beyond the rounding of its own arithmetic: the owner's solve
    /// residual times its efforts, which this boundary cannot bound by
    /// itself. A defect within it is booked as numerical loss.
    ///
    /// # Errors
    /// The owner's refusal of the netlist ([`CircuitComponentError::Circuit`]);
    /// a tolerance that is not positive and finite; an external port that
    /// does not drive exactly one source branch.
    pub fn new(
        graph: CircuitGraph,
        mutual: Vec<MutualInductance>,
        energy_tolerance_j: f64,
    ) -> Result<Self, CircuitComponentError> {
        if !(energy_tolerance_j.is_finite() && energy_tolerance_j > 0.0) {
            return Err(invalid("energy tolerance must be positive and finite"));
        }
        let reference = assemble_circuit_with_mutual_inductance(&graph, &mutual)?;
        let mut roles = vec![None; reference.ports];
        let mut resistors = Vec::new();
        for (branch, &(a, b, element)) in graph.branches.iter().enumerate() {
            let (port, role) = match element {
                Branch::Resistor { .. } => {
                    resistors.push(Resistor {
                        branch,
                        a,
                        b,
                        law: None,
                    });
                    continue;
                }
                Branch::VoltageSource { port } => (port, PortRole::EffortIn),
                Branch::CurrentSource { port } => (port, PortRole::FlowIn),
                Branch::Inductor { .. } | Branch::Capacitor { .. } => continue,
            };
            let slot = roles
                .get_mut(port)
                .ok_or_else(|| invalid("a source names a port the circuit does not have"))?;
            if slot.replace(role).is_some() {
                return Err(invalid("an external port drives more than one source"));
            }
        }
        let mut ports = Vec::with_capacity(roles.len() + resistors.len());
        for (port, role) in roles.into_iter().enumerate() {
            let role = role.ok_or_else(|| invalid("an external port drives no source"))?;
            let id = StableId::new(format!("source{port}"))
                .map_err(CircuitComponentError::Identifier)?;
            ports.push(
                PortSpec::power(
                    id,
                    PortKind::ElectricalVoltageCurrent,
                    role,
                    NonZeroUsize::MIN,
                )
                .map_err(CircuitComponentError::Port)?,
            );
        }
        for resistor in &resistors {
            let id = StableId::new(format!("heat{}", resistor.branch))
                .map_err(CircuitComponentError::Identifier)?;
            ports
                .push(PortSpec::heat(id, PortRole::EffortIn).map_err(CircuitComponentError::Port)?);
        }
        PortLayout::new(&ports).map_err(CircuitComponentError::Port)?;
        Ok(Self {
            graph,
            mutual,
            resistors,
            reference,
            ports,
            energy_tolerance_j,
        })
    }

    /// The same circuit with the resistor at `branch` following `law` at the
    /// temperature its heat port receives, in every proposal.
    ///
    /// # Errors
    /// `branch` is not a resistor of the graph, or already follows a law.
    pub fn following(
        mut self,
        branch: usize,
        law: Arc<dyn ResistanceLaw>,
    ) -> Result<Self, CircuitComponentError> {
        let resistor = self
            .resistors
            .iter_mut()
            .find(|resistor| resistor.branch == branch)
            .ok_or_else(|| invalid("a law is given to a branch that is not a resistor"))?;
        if resistor.law.replace(law).is_some() {
            return Err(invalid("a resistor follows one law only"));
        }
        Ok(self)
    }

    /// The circuit's graph, at the resistances it was built with.
    #[must_use]
    pub fn graph(&self) -> &CircuitGraph {
        &self.graph
    }

    /// A state from the inductors' flux linkages \[Wb\] in the owner's flux
    /// order and the capacitors' charges \[C\] in its charge order.
    ///
    /// # Errors
    /// A count that is not the circuit's, or a value that is not finite.
    pub fn state(
        &self,
        fluxes_wb: &[f64],
        charges_c: &[f64],
    ) -> Result<Vec<f64>, CircuitComponentError> {
        if fluxes_wb.len() != self.reference.flux_index.len()
            || charges_c.len() != self.reference.charge_index.len()
            || fluxes_wb
                .iter()
                .chain(charges_c)
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite flux per inductor and one finite charge per capacitor required",
            ));
        }
        Ok(fluxes_wb.iter().chain(charges_c).copied().collect())
    }

    /// The flux linkages \[Wb\] of a state, in the owner's flux order.
    #[must_use]
    pub fn fluxes_wb<'s>(&self, state: &'s [f64]) -> &'s [f64] {
        &state[..self.reference.flux_index.len().min(state.len())]
    }

    /// The charges \[C\] of a state, in the owner's charge order.
    #[must_use]
    pub fn charges_c<'s>(&self, state: &'s [f64]) -> &'s [f64] {
        &state[self.reference.flux_index.len().min(state.len())..]
    }

    fn dimension(&self) -> usize {
        self.reference.flux_index.len() + self.reference.charge_index.len()
    }

    /// The owner's extended state: the differential coordinates in their
    /// places and every algebraic multiplier at zero.
    fn extended(&self, dae: &CircuitDae, state: &[f64]) -> Result<Vec<f64>, CircuitComponentError> {
        if state.len() != self.dimension() || state.iter().any(|value| !value.is_finite()) {
            return Err(invalid("complete finite circuit state required"));
        }
        let mut x = vec![0.0; dae.system.state_dim()];
        let places = dae.flux_index.iter().chain(&dae.charge_index);
        for (&index, &value) in places.zip(state) {
            x[index] = value;
        }
        Ok(x)
    }

    fn energy(&self, state: &[f64]) -> Result<f64, CircuitComponentError> {
        let x = self.extended(&self.reference, state)?;
        let energy = self.reference.system.hamiltonian(&x);
        if energy.is_finite() {
            Ok(energy)
        } else {
            Err(invalid("the circuit's stored energy is not finite"))
        }
    }

    /// The graph with each resistor at the resistance of this step.
    fn graph_at(
        &self,
        temperatures: &[f64],
    ) -> Result<(CircuitGraph, Vec<f64>), CircuitComponentError> {
        let mut graph = self.graph.clone();
        let mut ohms = Vec::with_capacity(self.resistors.len());
        for (resistor, &kelvin) in self.resistors.iter().zip(temperatures) {
            let (_, _, element) = &mut graph.branches[resistor.branch];
            let Branch::Resistor { ohms: value } = element else {
                return Err(invalid("a resistor branch changed its element"));
            };
            if let Some(law) = &resistor.law {
                let answer = law
                    .ohms(kelvin)
                    .map_err(|reason| CircuitComponentError::Law {
                        branch: resistor.branch,
                        reason,
                    })?;
                if !(answer.is_finite() && answer > 0.0) {
                    return Err(CircuitComponentError::Law {
                        branch: resistor.branch,
                        reason: format!("{answer} ohm at {kelvin} K is not positive and finite"),
                    });
                }
                *value = answer;
            }
            ohms.push(*value);
        }
        Ok((graph, ohms))
    }
}

impl KirchhoffCircuit {
    /// Each resistor's Joule heat over the step at the algebraic midpoint of
    /// its node potentials, `ends` being the owner's extended state at the
    /// start and at the end. Writes minus each heat rate \[W\] to `rates` and
    /// returns the heat \[J\] and the entropy \[J/K\] it generates at the
    /// temperatures its ports received.
    fn joule_heat(
        &self,
        dae: &CircuitDae,
        ends: [&[f64]; 2],
        ohms: &[f64],
        dt_s: f64,
        kelvin: &[f64],
        rates: &mut [f64],
    ) -> Result<(f64, f64), CircuitComponentError> {
        let potential = |node: usize| {
            if node == 0 {
                0.0
            } else {
                let index = dae.node_potential_index[node - 1];
                f64::midpoint(ends[0][index], ends[1][index])
            }
        };
        let mut heat_j = 0.0;
        let mut entropy = 0.0;
        for (k, (resistor, &r)) in self.resistors.iter().zip(ohms).enumerate() {
            let volts = potential(resistor.a) - potential(resistor.b);
            let amperes = volts / r;
            let power_w = volts * amperes;
            let joules = dt_s * power_w;
            if !(power_w.is_finite() && joules.is_finite()) || power_w < 0.0 {
                return Err(invalid(
                    "a resistor's Joule heat is not finite and nonnegative",
                ));
            }
            rates[k] = -power_w;
            heat_j += joules;
            entropy += joules / kelvin[k];
        }
        Ok((heat_j, entropy))
    }
}

/// Words of a branch element: a tag and its value's bits.
fn element_words(element: Branch, words: &mut Vec<u64>) {
    match element {
        Branch::Resistor { ohms } => words.extend([0, ohms.to_bits()]),
        Branch::Inductor { henries } => words.extend([1, henries.to_bits()]),
        Branch::Capacitor { farads } => words.extend([2, farads.to_bits()]),
        Branch::VoltageSource { port } => words.extend([3, port as u64]),
        Branch::CurrentSource { port } => words.extend([4, port as u64]),
    }
}

impl Component for KirchhoffCircuit {
    type State = Vec<f64>;
    type Error = CircuitComponentError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::circuit::KirchhoffCircuit",
        state_family: "fs_phs::circuit differential coordinates",
        schema_grammar: "Vec<f64>: inductor flux linkages in flux order, then capacitor charges in charge order",
        codec_grammar: "one word per coordinate, the f64 bits, in schema order",
        codec_version: 1,
    });

    fn encode_state(&self, state: &Vec<f64>, words: &mut Vec<u64>) -> bool {
        words.extend(state.iter().map(|value| value.to_bits()));
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Vec<f64>> {
        let state: Vec<f64> = words.iter().map(|&word| f64::from_bits(word)).collect();
        (state.len() == self.dimension() && state.iter().all(|value| value.is_finite()))
            .then_some(state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.graph.node_count as u64);
        words.push(self.graph.branches.len() as u64);
        for &(a, b, element) in &self.graph.branches {
            words.extend([a as u64, b as u64]);
            element_words(element, words);
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
        words.push(self.mutual.len() as u64);
        for coupling in &self.mutual {
            words.extend([
                coupling.branches.0 as u64,
                coupling.branches.1 as u64,
                coupling.henries.to_bits(),
            ]);
        }
        for resistor in &self.resistors {
            if let Some(law) = &resistor.law {
                words.push(resistor.branch as u64);
                law.identity_words(words);
            }
        }
        words.push(self.energy_tolerance_j.to_bits());
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, CircuitComponentError> {
        self.energy(state)
    }

    fn propose(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Vec<f64>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, CircuitComponentError> {
        cx.checkpoint()
            .map_err(|_| CircuitComponentError::Cancelled)?;
        let sources = self.reference.ports;
        let heats = self.resistors.len();
        if !(dt_s.is_finite() && dt_s > 0.0)
            || inputs.len() != sources + heats
            || outputs.len() != sources + heats
            || inputs.iter().any(|value| !value.is_finite())
            || inputs[sources..].iter().any(|&kelvin| kelvin <= 0.0)
        {
            return Err(invalid(
                "positive duration, complete finite ports and positive temperatures required",
            ));
        }
        let (graph, ohms) = self.graph_at(&inputs[sources..])?;
        let dae = assemble_circuit_with_mutual_inductance(&graph, &self.mutual)?;
        if dae.ports != sources
            || dae.flux_index != self.reference.flux_index
            || dae.charge_index != self.reference.charge_index
        {
            return Err(invalid(
                "the circuit's coordinates changed with its resistances",
            ));
        }
        let x0 = self.extended(&dae, old)?;
        cx.checkpoint()
            .map_err(|_| CircuitComponentError::Cancelled)?;
        let (record, defect) = dae.step_audited(&x0, &inputs[..sources], dt_s)?;
        cx.checkpoint()
            .map_err(|_| CircuitComponentError::Cancelled)?;
        if record.x.len() != x0.len() || record.y.len() != sources {
            return Err(invalid(
                "the owner's step has another shape than the circuit",
            ));
        }

        let (heat_j, entropy) = self.joule_heat(
            &dae,
            [&x0, &record.x],
            &ohms,
            dt_s,
            &inputs[sources..],
            &mut outputs[sources..],
        )?;
        let next: Vec<f64> = dae
            .flux_index
            .iter()
            .chain(&dae.charge_index)
            .map(|&index| record.x[index])
            .collect();
        let energy_old = self.energy(old)?;
        let energy_new = self.energy(&next)?;
        let stored_change_j = energy_new - energy_old;
        let supplied_j = record.supplied;
        let numerical_loss_j = supplied_j - stored_change_j - heat_j;
        let magnitude = energy_old.abs()
            + energy_new.abs()
            + supplied_j.abs()
            + record.dissipated.abs()
            + heat_j;
        let rounding = arithmetic_rounding_floor(
            32_usize
                .saturating_add(8_usize.saturating_mul(x0.len()))
                .saturating_add(8_usize.saturating_mul(sources + heats)),
            magnitude,
        );
        let allowance = rounding + self.energy_tolerance_j;
        if ![rounding, numerical_loss_j, entropy, defect]
            .iter()
            .all(|value| value.is_finite())
        {
            return Err(invalid("the circuit's energy account is not finite"));
        }
        if numerical_loss_j.abs() > allowance {
            return Err(CircuitComponentError::Balance {
                defect_j: numerical_loss_j,
                tolerance_j: allowance,
            });
        }
        let (residual, tolerance) = match &record.solver_certificate {
            Some(certificate) => (certificate.residual, certificate.effective_budget),
            None => (record.solver_residual, record.solver_residual),
        };
        cx.checkpoint()
            .map_err(|_| CircuitComponentError::Cancelled)?;
        outputs[..sources].copy_from_slice(&record.y);
        new.clone_from(&next);
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j: heat_j,
            numerical_loss_j,
            entropy_generated_j_per_k: entropy,
            energy_tolerance_j: rounding,
            residual,
            tolerance,
        }))
    }
}
