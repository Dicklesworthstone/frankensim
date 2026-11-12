//! Executable Component contract, version 0.
//!
//! A Component is one physics owner seen from outside: ports, an accepted
//! state the system owns, and a pure step from an accepted state to a
//! candidate. Coupled behaviour comes from connecting components, never from
//! code written for one pairing.
//!
//! The contract is tiered because most owners offer "accepted state in,
//! candidate and receipt out" and only some offer derivatives:
//!
//! - [`Component`] is mandatory: `propose` and `admit`.
//! - [`PortTangent`] adds the derivative of the port map, for interface Newton
//!   and for adjoints.
//! - [`StageResidual`] adds a residual in the component's own unknowns, for
//!   monolithic solves.
//!
//! # Port values during one step
//!
//! Each port carries one representative pair per step. The side the component
//! receives is held constant over the step. The side it returns is the step
//! average of the conjugate variable. With that definition the energy that
//! entered the owner through a port is computed by the kernel from the two
//! values on the wire ([`PortSpec::energy_in_j`]), so both ends of a lossless
//! connection book the same number and interface energy cancels by
//! construction. A component whose integrator returns something other than
//! the step average must book the difference as numerical loss in its
//! [`Receipt`]; the kernel checks that the receipt explains the defect.
//!
//! Sign: flow, heat rate and a stream's rates are positive into the owner.
//!
//! # What v0 does not cover
//!
//! Entropy carried by streams, multistage port values, and snapshot identity
//! are later beads
//! (`frankensim-rc-1109-root-3gc12.5.2`, `.5.4`). The contract stays v0 until
//! the five kernel slices pass on it.

pub mod conformance;

use core::fmt;
use core::num::NonZeroUsize;
use core::ops::Range;

use fs_exec::Cx;

use crate::{PortKind, StableId};

/// Version of the executable Component contract.
pub const COMPONENT_CONTRACT_VERSION: u16 = 0;

/// What a port carries during one step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Carrier {
    /// A power-conjugate pair of the given kind; power into the owner is the
    /// component-wise product of effort and flow, summed, or `e^T M f` when
    /// the port declares a [`DualityMetric`].
    Power(PortKind),
    /// Heat. Effort is absolute temperature [K]; flow is heat rate [W] into
    /// the owner. Power into the owner is the flow itself, and entropy flow is
    /// flow divided by effort. This is the wire form of the thermal port kind:
    /// owners report heat and temperature, the kernel derives entropy.
    Heat,
    /// A one-way signal. It carries no power and no conserved quantity.
    Signal,
    /// A material stream between a storage and a transport element.
    ///
    /// The side that receives flow (a storage: a vessel, a cavity) returns
    /// what its contents are like: pressure [Pa] at [`STREAM_PRESSURE`],
    /// specific enthalpy [J/kg] at [`STREAM_SPECIFIC_ENTHALPY`], then one
    /// value per unit mass for each carried quantity (a species' mass
    /// fraction). The side that receives those (a transport element: a
    /// restriction, a valve, a pump) returns what crosses the port, positive
    /// into itself: mass rate [kg/s] at [`STREAM_MASS_RATE`], enthalpy rate
    /// [W] at [`STREAM_ENTHALPY_RATE`], then the rate of each carried
    /// quantity.
    ///
    /// Power into the owner is the enthalpy rate itself, as the heat rate is
    /// for [`Carrier::Heat`]; it is not a product of the two sides. What a
    /// stream delivers is set by where it comes from, so a transport element
    /// multiplies its mass rate by the values of the side the mass leaves
    /// (the donor rule). Both directions and zero flow are then one law, and
    /// what leaves one storage is what enters the other.
    ///
    /// Entropy carried by a stream is not on the wire in v0, so entropy
    /// accounts do not include it.
    Stream,
}

/// Index of the pressure [Pa] on the storage side of a stream port.
pub const STREAM_PRESSURE: usize = 0;
/// Index of the specific enthalpy [J/kg] on the storage side of a stream port.
pub const STREAM_SPECIFIC_ENTHALPY: usize = 1;
/// Index of the mass rate [kg/s] on the transport side of a stream port.
pub const STREAM_MASS_RATE: usize = 0;
/// Index of the enthalpy rate [W] on the transport side of a stream port.
pub const STREAM_ENTHALPY_RATE: usize = 1;
/// Index of the first carried quantity on either side of a stream port.
pub const STREAM_CARRIED: usize = 2;

/// Which side of a port the component receives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortRole {
    /// Receives effort, returns flow.
    EffortIn,
    /// Receives flow, returns effort.
    FlowIn,
    /// Receives a signal.
    SignalIn,
    /// Emits a signal.
    SignalOut,
}

/// A symmetric positive-definite duality matrix for field coefficients.
///
/// A field port uses the same number of effort and flow coefficients, with
/// power `e^T M f`. Entries are row-major. Construction validates a finite,
/// exactly symmetric matrix and retains its Cholesky factor for connection
/// admission; multiplication and solves allocate nothing. This covers trace
/// mass matrices in one basis, not rectangular or indefinite dual pairings.
#[derive(Debug, Clone)]
pub struct DualityMetric {
    width: NonZeroUsize,
    entries: Vec<f64>,
    lower: Vec<f64>,
}

impl PartialEq for DualityMetric {
    fn eq(&self, other: &Self) -> bool {
        self.width == other.width && self.entries == other.entries
    }
}

// Construction excludes NaN. Equality concerns the metric, not its cached
// factor, and ordinary floating-point equality equates signed zeros.
impl Eq for DualityMetric {}

impl DualityMetric {
    /// Admit a square field metric and retain its Cholesky factor.
    ///
    /// Positive finite pivots and finite factors are required. This is a
    /// floating-point Cholesky admission, not an interval proof of positive
    /// definiteness or a bound on conditioning.
    ///
    /// # Errors
    /// Refuses the wrong shape, nonfinite coefficients, unequal mirrored
    /// entries, nonpositive pivots, or an unrepresentable factorization.
    #[allow(clippy::float_cmp)] // Symmetry is an exact admission rule.
    pub fn new(width: NonZeroUsize, entries: Vec<f64>) -> Result<Self, ContractViolation> {
        let n = width.get();
        if n.checked_mul(n) != Some(entries.len()) {
            return Err(ContractViolation::MetricShape {
                width: n,
                entries: entries.len(),
            });
        }
        for (index, &value) in entries.iter().enumerate() {
            if !value.is_finite() {
                return Err(ContractViolation::MetricCoefficient { index, value });
            }
        }
        for row in 0..n {
            for column in 0..row {
                if entries[row * n + column] != entries[column * n + row] {
                    return Err(ContractViolation::MetricAsymmetric { row, column });
                }
            }
        }
        let mut lower = vec![0.0; entries.len()];
        for row in 0..n {
            for column in 0..=row {
                let mut value = entries[row * n + column];
                for k in 0..column {
                    value -= lower[row * n + k] * lower[column * n + k];
                }
                if !value.is_finite() {
                    return Err(ContractViolation::MetricFactorization { row, column, value });
                }
                let factor = if row == column {
                    if value <= 0.0 {
                        return Err(ContractViolation::MetricPivot { index: row, value });
                    }
                    value.sqrt()
                } else {
                    value / lower[column * n + column]
                };
                if !factor.is_finite() {
                    return Err(ContractViolation::MetricFactorization {
                        row,
                        column,
                        value: factor,
                    });
                }
                lower[row * n + column] = factor;
            }
        }
        Ok(Self {
            width,
            entries,
            lower,
        })
    }

    /// Number of effort and flow coefficients.
    #[must_use]
    pub const fn width(&self) -> usize {
        self.width.get()
    }

    /// Matrix coefficients in row-major order.
    #[must_use]
    pub fn entries(&self) -> &[f64] {
        &self.entries
    }

    fn check_vectors(&self, input: &[f64], output: &[f64]) -> Result<(), ContractViolation> {
        for (vector, actual) in [
            ("metric input", input.len()),
            ("metric output", output.len()),
        ] {
            if actual != self.width() {
                return Err(ContractViolation::VectorLength {
                    vector,
                    expected: self.width(),
                    actual,
                });
            }
        }
        if !input.iter().all(|value| value.is_finite()) {
            return Err(ContractViolation::NonFinite {
                quantity: "metric input",
            });
        }
        Ok(())
    }

    /// Write `M * input` into a caller-owned buffer without allocation.
    ///
    /// # Errors
    /// Refuses wrong vector lengths, nonfinite input, or a nonfinite result.
    /// The output buffer is unspecified after an arithmetic refusal.
    pub fn multiply(&self, input: &[f64], output: &mut [f64]) -> Result<(), ContractViolation> {
        self.check_vectors(input, output)?;
        let n = self.width();
        for (row, slot) in output.iter_mut().enumerate() {
            let mut value = 0.0;
            for (column, &coefficient) in input.iter().enumerate() {
                value += self.entries[row * n + column] * coefficient;
            }
            if !value.is_finite() {
                return Err(ContractViolation::NonFinite {
                    quantity: "metric product",
                });
            }
            *slot = value;
        }
        Ok(())
    }

    /// Solve `M * solution = rhs` using the retained factor without allocation.
    ///
    /// # Errors
    /// Refuses wrong vector lengths, nonfinite input, or a nonfinite solve.
    /// The solution buffer is unspecified after an arithmetic refusal.
    pub fn solve(&self, rhs: &[f64], solution: &mut [f64]) -> Result<(), ContractViolation> {
        self.check_vectors(rhs, solution)?;
        let n = self.width();
        for (row, &right) in rhs.iter().enumerate() {
            let mut value = right;
            for (column, &solved) in solution[..row].iter().enumerate() {
                value -= self.lower[row * n + column] * solved;
            }
            solution[row] = value / self.lower[row * n + row];
            if !solution[row].is_finite() {
                return Err(ContractViolation::NonFinite {
                    quantity: "metric forward solve",
                });
            }
        }
        for row in (0..n).rev() {
            let mut value = solution[row];
            for (column, &solved) in solution.iter().enumerate().skip(row + 1) {
                value -= self.lower[column * n + row] * solved;
            }
            solution[row] = value / self.lower[row * n + row];
            if !solution[row].is_finite() {
                return Err(ContractViolation::NonFinite {
                    quantity: "metric backward solve",
                });
            }
        }
        Ok(())
    }

    fn pairing(&self, effort: &[f64], flow: &[f64]) -> f64 {
        let n = self.width();
        let mut power = 0.0;
        for (row, &e) in effort.iter().enumerate() {
            let mut dual_flow = 0.0;
            for (column, &f) in flow.iter().enumerate() {
                dual_flow += self.entries[row * n + column] * f;
            }
            power += e * dual_flow;
        }
        power
    }
}

/// Executable description of one port of a component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortSpec {
    id: StableId,
    carrier: Carrier,
    role: PortRole,
    width: NonZeroUsize,
    field_metric: Option<DualityMetric>,
}

impl PortSpec {
    /// A power port of `width` scalar coefficients.
    ///
    /// # Errors
    /// [`ContractViolation::PortRole`] if `role` is a signal role.
    pub fn power(
        id: StableId,
        kind: PortKind,
        role: PortRole,
        width: NonZeroUsize,
    ) -> Result<Self, ContractViolation> {
        if !matches!(role, PortRole::EffortIn | PortRole::FlowIn) {
            return Err(ContractViolation::PortRole {
                port: id.as_str().to_owned(),
            });
        }
        Ok(Self {
            id,
            carrier: Carrier::Power(kind),
            role,
            width,
            field_metric: None,
        })
    }

    /// A field power port with energy `dt * effort^T * metric * flow`.
    ///
    /// Effort and flow have the same width, given by the admitted metric.
    ///
    /// # Errors
    /// [`ContractViolation::PortRole`] if `role` is a signal role.
    pub fn field_power(
        id: StableId,
        kind: PortKind,
        role: PortRole,
        metric: DualityMetric,
    ) -> Result<Self, ContractViolation> {
        let mut port = Self::power(id, kind, role, metric.width)?;
        port.field_metric = Some(metric);
        Ok(port)
    }

    /// A scalar heat port: temperature [K] and heat rate [W] into the owner.
    ///
    /// # Errors
    /// [`ContractViolation::PortRole`] if `role` is a signal role.
    pub fn heat(id: StableId, role: PortRole) -> Result<Self, ContractViolation> {
        if !matches!(role, PortRole::EffortIn | PortRole::FlowIn) {
            return Err(ContractViolation::PortRole {
                port: id.as_str().to_owned(),
            });
        }
        Ok(Self {
            id,
            carrier: Carrier::Heat,
            role,
            width: NonZeroUsize::MIN,
            field_metric: None,
        })
    }

    /// A stream port that carries `carried` quantities besides mass and
    /// enthalpy (see [`Carrier::Stream`]). A storage has `PortRole::FlowIn`,
    /// a transport element `PortRole::EffortIn`.
    ///
    /// # Errors
    /// [`ContractViolation::PortRole`] if `role` is a signal role.
    pub fn stream(id: StableId, role: PortRole, carried: usize) -> Result<Self, ContractViolation> {
        if !matches!(role, PortRole::EffortIn | PortRole::FlowIn) {
            return Err(ContractViolation::PortRole {
                port: id.as_str().to_owned(),
            });
        }
        Ok(Self {
            id,
            carrier: Carrier::Stream,
            role,
            width: NonZeroUsize::MIN.saturating_add(STREAM_CARRIED - 1 + carried),
            field_metric: None,
        })
    }

    /// A signal input of `width` coefficients.
    #[must_use]
    pub fn signal_in(id: StableId, width: NonZeroUsize) -> Self {
        Self {
            id,
            carrier: Carrier::Signal,
            role: PortRole::SignalIn,
            width,
            field_metric: None,
        }
    }

    /// A signal output of `width` coefficients.
    #[must_use]
    pub fn signal_out(id: StableId, width: NonZeroUsize) -> Self {
        Self {
            id,
            carrier: Carrier::Signal,
            role: PortRole::SignalOut,
            width,
            field_metric: None,
        }
    }

    /// Identifier, unique within the component.
    #[must_use]
    pub fn id(&self) -> &StableId {
        &self.id
    }

    /// What the port carries.
    #[must_use]
    pub fn carrier(&self) -> Carrier {
        self.carrier
    }

    /// Which side the component receives.
    #[must_use]
    pub fn role(&self) -> PortRole {
        self.role
    }

    /// Number of scalar coefficients on each side.
    #[must_use]
    pub fn width(&self) -> usize {
        self.width.get()
    }

    /// The explicit field pairing, or `None` for an ordinary port.
    #[must_use]
    pub fn field_metric(&self) -> Option<&DualityMetric> {
        self.field_metric.as_ref()
    }

    /// Coefficients this port contributes to the input vector.
    #[must_use]
    pub fn input_width(&self) -> usize {
        match self.role {
            PortRole::EffortIn | PortRole::FlowIn | PortRole::SignalIn => self.width.get(),
            PortRole::SignalOut => 0,
        }
    }

    /// Coefficients this port contributes to the output vector.
    #[must_use]
    pub fn output_width(&self) -> usize {
        match self.role {
            PortRole::EffortIn | PortRole::FlowIn | PortRole::SignalOut => self.width.get(),
            PortRole::SignalIn => 0,
        }
    }

    /// Energy [J] that entered the owner through this port during a step of
    /// `dt_s`, from the port's own input and output coefficients.
    ///
    /// Signals return zero. The slices must have this port's widths; a length
    /// mismatch is a caller bug and yields NaN, which every audit refuses.
    #[must_use]
    pub fn energy_in_j(&self, dt_s: f64, input: &[f64], output: &[f64]) -> f64 {
        let width = self.width.get();
        match (self.carrier, self.role) {
            (Carrier::Signal, _) | (_, PortRole::SignalIn | PortRole::SignalOut) => 0.0,
            (Carrier::Power(_), _) => {
                if input.len() != width || output.len() != width {
                    return f64::NAN;
                }
                if let Some(metric) = &self.field_metric {
                    let (effort, flow) = if self.role == PortRole::EffortIn {
                        (input, output)
                    } else {
                        (output, input)
                    };
                    return dt_s * metric.pairing(effort, flow);
                }
                let mut sum = 0.0;
                for (a, b) in input.iter().zip(output) {
                    sum += a * b;
                }
                dt_s * sum
            }
            (Carrier::Heat, role) => {
                let flow = if role == PortRole::FlowIn {
                    input
                } else {
                    output
                };
                if flow.len() != width {
                    return f64::NAN;
                }
                dt_s * flow.iter().sum::<f64>()
            }
            (Carrier::Stream, _) => match self.stream_flow(input, output) {
                Some(flow) => dt_s * flow[STREAM_ENTHALPY_RATE],
                None => f64::NAN,
            },
        }
    }

    /// The side of a stream port that holds what crosses it, positive into
    /// the owner; `None` for another carrier or slices of the wrong length.
    fn stream_flow<'a>(&self, input: &'a [f64], output: &'a [f64]) -> Option<&'a [f64]> {
        let flow = match (self.carrier, self.role) {
            (Carrier::Stream, PortRole::FlowIn) => input,
            (Carrier::Stream, PortRole::EffortIn) => output,
            _ => return None,
        };
        (flow.len() == self.width.get()).then_some(flow)
    }

    /// Mass [kg] that entered the owner through this port during a step.
    ///
    /// Only streams carry mass; every other carrier returns zero. A stream
    /// port given slices of the wrong length yields NaN.
    #[must_use]
    pub fn mass_in_kg(&self, dt_s: f64, input: &[f64], output: &[f64]) -> f64 {
        if self.carrier != Carrier::Stream {
            return 0.0;
        }
        self.stream_flow(input, output)
            .map_or(f64::NAN, |flow| dt_s * flow[STREAM_MASS_RATE])
    }

    /// Entropy [J/K] that entered the owner through this port during a step.
    ///
    /// Only heat ports carry entropy on the wire; every other carrier returns
    /// zero. A non-positive temperature yields NaN.
    #[must_use]
    pub fn entropy_in_j_per_k(&self, dt_s: f64, input: &[f64], output: &[f64]) -> f64 {
        if self.carrier != Carrier::Heat {
            return 0.0;
        }
        let (temperature, flow) = match self.role {
            PortRole::EffortIn => (input, output),
            PortRole::FlowIn => (output, input),
            PortRole::SignalIn | PortRole::SignalOut => return 0.0,
        };
        match (temperature.first(), flow.first()) {
            (Some(&t), Some(&q)) if t > 0.0 => dt_s * q / t,
            _ => f64::NAN,
        }
    }
}

/// Offsets of each port inside a component's input and output vectors.
///
/// Ports are laid out in declaration order. Building the layout also admits
/// the port list: identifiers are unique and carriers agree with roles.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortLayout {
    input_offsets: Vec<usize>,
    output_offsets: Vec<usize>,
}

impl PortLayout {
    /// Lay out `ports` in order.
    ///
    /// # Errors
    /// [`ContractViolation::DuplicatePort`] for a repeated identifier and
    /// [`ContractViolation::PortRole`] for a carrier that disagrees with its
    /// role.
    pub fn new(ports: &[PortSpec]) -> Result<Self, ContractViolation> {
        let mut input_offsets = Vec::with_capacity(ports.len() + 1);
        let mut output_offsets = Vec::with_capacity(ports.len() + 1);
        let (mut input, mut output) = (0_usize, 0_usize);
        for (index, port) in ports.iter().enumerate() {
            if ports[..index].iter().any(|earlier| earlier.id == port.id) {
                return Err(ContractViolation::DuplicatePort {
                    port: port.id.as_str().to_owned(),
                });
            }
            let signal_role = matches!(port.role, PortRole::SignalIn | PortRole::SignalOut);
            if signal_role != (port.carrier == Carrier::Signal) {
                return Err(ContractViolation::PortRole {
                    port: port.id.as_str().to_owned(),
                });
            }
            input_offsets.push(input);
            output_offsets.push(output);
            input += port.input_width();
            output += port.output_width();
        }
        input_offsets.push(input);
        output_offsets.push(output);
        Ok(Self {
            input_offsets,
            output_offsets,
        })
    }

    /// Number of ports.
    #[must_use]
    pub fn ports(&self) -> usize {
        self.input_offsets.len() - 1
    }

    /// Length of the input vector.
    #[must_use]
    pub fn input_len(&self) -> usize {
        self.input_offsets[self.input_offsets.len() - 1]
    }

    /// Length of the output vector.
    #[must_use]
    pub fn output_len(&self) -> usize {
        self.output_offsets[self.output_offsets.len() - 1]
    }

    /// Range of `port` inside the input vector (empty for a signal output).
    #[must_use]
    pub fn input_range(&self, port: usize) -> Range<usize> {
        self.input_offsets[port]..self.input_offsets[port + 1]
    }

    /// Range of `port` inside the output vector (empty for a signal input).
    #[must_use]
    pub fn output_range(&self, port: usize) -> Range<usize> {
        self.output_offsets[port]..self.output_offsets[port + 1]
    }
}

/// How a component's port map behaves, which decides how a block is solved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class {
    /// Outputs are an affine function of inputs at a fixed state and step,
    /// and the slopes of that function depend on the step length alone,
    /// not on the state or the inputs. A component whose slopes move with
    /// its state declares [`Class::Smooth`].
    ExactLinear,
    /// Outputs are a differentiable function of inputs.
    Smooth,
    /// A scalar monotone relation; safeguarded bracketing converges.
    Monotone1d,
    /// Set-valued or switching (contact, stick-slip, phase change).
    NonSmooth,
}

/// A component's own account of one proposed step.
///
/// Port energies are not here: the kernel computes them from the wire values.
/// The first law for the step is
/// `sum(port energy in) = stored_change_j + numerical_loss_j`
/// within `energy_tolerance_j`. Physical dissipation is energy turned into
/// heat; it either stays in the component as stored thermal energy or leaves
/// through a heat port, so it is not a separate term of that balance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Receipt {
    /// Change of stored energy [J] as the owner computes it, without
    /// cancellation where the owner can avoid it.
    pub stored_change_j: f64,
    /// Physical dissipation [J] during the step. Never negative.
    pub dissipated_j: f64,
    /// Energy [J] the integrator lost (positive) or created (negative) that
    /// no physical law accounts for.
    pub numerical_loss_j: f64,
    /// Entropy generated inside the component [J/K]. Never negative.
    pub entropy_generated_j_per_k: f64,
    /// Absolute tolerance [J] the owner claims for its energy balance and
    /// for the re-measurement of its endpoints.
    ///
    /// It must cover the roundoff of [`Component::stored_energy_j`] at both
    /// ends of the step. When the stored energy is a sum of terms that cancel
    /// (kinetic plus a negative potential), that roundoff scales with the
    /// size of the terms, not of the total, and only the owner knows it. The
    /// audit trusts this number: a tolerance claimed too wide hides a leak
    /// from the audit, though not from the ledger's accepted-defect account.
    pub energy_tolerance_j: f64,
    /// Norm of the owner's solver residual; zero for a closed-form step.
    pub residual: f64,
    /// Tolerance the residual was required to meet.
    pub tolerance: f64,
}

impl Receipt {
    /// A receipt for a closed-form step with no dissipation and no loss.
    #[must_use]
    pub fn lossless(stored_change_j: f64, energy_tolerance_j: f64) -> Self {
        Self {
            stored_change_j,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j,
            residual: 0.0,
            tolerance: 0.0,
        }
    }
}

/// What a component answers to a proposed step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Proposal {
    /// The candidate state and outputs are filled for the whole step.
    Step(Receipt),
    /// The component will not take this step; retry every component with a
    /// step no longer than `dt_s` (an event was located, or accuracy failed).
    /// The candidate and the outputs are unspecified.
    Shorten {
        /// Longest step the component will take from this state [s].
        dt_s: f64,
    },
}

/// The mandatory tier: a pure step from an accepted state to a candidate.
///
/// The system owns every state. `propose` reads the accepted state and writes
/// a candidate slot; nothing is published until every component of the system
/// has proposed and admitted.
pub trait Component {
    /// Accepted state. Owned by the system, never by the component.
    type State: Clone;
    /// The owner's typed refusal. It is carried to the caller unchanged.
    type Error: fmt::Debug + 'static;

    /// Ports in a stable order; index `i` is the same port in every call.
    fn ports(&self) -> &[PortSpec];

    /// How the port map behaves.
    fn class(&self) -> Class;

    /// Stored energy [J] of a state, measured from the state alone.
    ///
    /// The kernel uses it to re-measure endpoints independently of receipts.
    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error>;

    /// Largest step not above `proposed_s` the component accepts from `state`.
    fn limit_dt_s(&self, _state: &Self::State, proposed_s: f64) -> Result<f64, Self::Error> {
        Ok(proposed_s)
    }

    /// Time [s] from `state` to the component's next scheduled event, if it
    /// has one. The driver never steps past it.
    fn next_event_s(&self, _state: &Self::State) -> Option<f64> {
        None
    }

    /// Propose one step of `dt_s` from `old` under held `inputs`.
    ///
    /// `inputs` and `outputs` follow [`PortLayout`]. `new` is a candidate slot
    /// whose content on entry is unspecified: write all of it, read none of
    /// it. The call must be pure: equal arguments give bitwise-equal results,
    /// and nothing observable changes outside `new` and `outputs`. Poll `cx`
    /// at bounded intervals and return the owner's cancellation error.
    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error>;

    /// Last chance to refuse, called after every partner proposed with the
    /// final interface values and before anything is published.
    fn admit(
        &self,
        _old: &Self::State,
        _new: &Self::State,
        _receipt: &Receipt,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
}

/// Optional tier: derivative of the port map of [`Component::propose`].
///
/// The map is `outputs = f(inputs)` at fixed `old` and `dt_s`. A type cannot
/// offer a tangent without the mandatory tier:
///
/// ```compile_fail
/// use fs_exec::Cx;
/// use fs_ports::component::PortTangent;
/// struct NotAComponent;
/// impl PortTangent for NotAComponent {
///     fn jvp(&self, _: &(), _: f64, _: &[f64], _: &[f64], _: &mut [f64], _: &Cx<'_>)
///         -> Result<(), ()> { Ok(()) }
///     fn vjp(&self, _: &(), _: f64, _: &[f64], _: &[f64], _: &mut [f64], _: &Cx<'_>)
///         -> Result<(), ()> { Ok(()) }
/// }
/// ```
pub trait PortTangent: Component {
    /// `d_outputs = (df/d inputs) d_inputs`.
    fn jvp(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error>;

    /// `bar_inputs = (df/d inputs)^T bar_outputs`.
    fn vjp(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error>;
}

/// Optional tier: a relation at power ports that is not a function of the
/// flow, solved on the load line its neighbours present.
///
/// Dry friction holds any force up to its limit at zero slip, and a closed
/// contact pushes with whatever keeps its gap closed: the effort such a
/// component returns is not a function of the flow it receives, so a solver
/// cannot iterate on that flow. What its neighbours do with an effort is a
/// function, and frozen at an iterate it is affine:
///
/// ```text
/// flow = released - compliance * effort
/// ```
///
/// `released` is the flow the ports would have if the component returned no
/// effort, and `compliance` says how much of it each unit of effort takes
/// away. The component finds the point of that line its own relation
/// allows, and steps there.
///
/// A type cannot offer this tier without the mandatory one:
///
/// ```compile_fail
/// use fs_exec::Cx;
/// use fs_ports::component::{LoadLine, Proposal};
/// struct NotAComponent;
/// impl LoadLine for NotAComponent {
///     fn line_ports(&self) -> &[usize] { &[] }
///     fn propose_on_line(&self, _: &(), _: f64, _: &[f64], _: &[f64], _: &mut (),
///         _: &mut [f64], _: &Cx<'_>) -> Result<Proposal, ()> { unimplemented!() }
/// }
/// ```
pub trait LoadLine: Component {
    /// The ports the line spans, as indices into [`Component::ports`], in
    /// increasing order and at least one. Each is a power port of role
    /// [`PortRole::FlowIn`]: it receives a flow and returns an effort. Their
    /// coefficients, in port order, are the line's flows and efforts.
    fn line_ports(&self) -> &[usize];

    /// [`Component::propose`] on the load line.
    ///
    /// `inputs` is as for `propose`, except that at the line's ports it
    /// holds the released flows. `compliance` is the square matrix of the
    /// line, row-major over the line's coefficients, so that
    /// `flow[i] = released[i] - sum(compliance[i][j] * effort[j])`. It need
    /// not be symmetric; its symmetric part is not negative. A coefficient
    /// whose flow is prescribed has a zero row and column.
    ///
    /// The efforts the component settles on are its outputs at the line's
    /// ports. The candidate, the other outputs and the receipt are those of
    /// the step at the flows on the line, not at the released ones. With
    /// `compliance` all zero the call is `propose`, to the bit. The purity
    /// rules of `propose` hold.
    #[allow(clippy::too_many_arguments)]
    fn propose_on_line(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error>;
}

/// Optional tier: the step as a residual in the component's own unknowns.
///
/// `residual(z) = 0` defines the step for held `inputs`. A monolithic solve
/// assembles these residuals with the connection equations.
///
/// ```compile_fail
/// use fs_exec::Cx;
/// use fs_ports::component::{Receipt, StageResidual};
/// struct NotAComponent;
/// impl StageResidual for NotAComponent {
///     fn unknowns(&self) -> usize { 0 }
///     fn initial_guess(&self, _: &(), _: f64, _: &[f64], _: &mut [f64]) -> Result<(), ()> { Ok(()) }
///     fn residual(&self, _: &(), _: f64, _: &[f64], _: &[f64], _: &mut [f64], _: &mut [f64], _: &Cx<'_>)
///         -> Result<(), ()> { Ok(()) }
///     fn jvp_z(&self, _: &(), _: f64, _: &[f64], _: &[f64], _: &[f64], _: &mut [f64], _: &mut [f64], _: &Cx<'_>)
///         -> Result<(), ()> { Ok(()) }
///     fn vjp_z(&self, _: &(), _: f64, _: &[f64], _: &[f64], _: &[f64], _: &[f64], _: &mut [f64], _: &Cx<'_>)
///         -> Result<(), ()> { Ok(()) }
///     fn accept(&self, _: &(), _: f64, _: &[f64], _: &[f64], _: &mut ()) -> Result<Receipt, ()> {
///         Ok(Receipt::lossless(0.0, 0.0))
///     }
/// }
/// ```
pub trait StageResidual: Component {
    /// Number of unknowns `z`.
    fn unknowns(&self) -> usize;

    /// A starting `z` for the step from `old`.
    fn initial_guess(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), Self::Error>;

    /// Residual `r(z)` and the port outputs implied by `z`.
    #[allow(clippy::too_many_arguments)]
    fn residual(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error>;

    /// Directional derivative of residual and outputs along `dz` at `z`.
    #[allow(clippy::too_many_arguments)]
    fn jvp_z(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error>;

    /// Transposed action: `bar_z = (dr/dz)^T bar_r + (d outputs/dz)^T bar_outputs`.
    #[allow(clippy::too_many_arguments)]
    fn vjp_z(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error>;

    /// Build the candidate state and receipt from converged unknowns.
    fn accept(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        new: &mut Self::State,
    ) -> Result<Receipt, Self::Error>;
}

/// A step that broke the contract, named by the law it broke.
#[derive(Debug, Clone, PartialEq)]
pub enum ContractViolation {
    /// Two ports of one component share an identifier.
    DuplicatePort {
        /// The repeated identifier.
        port: String,
    },
    /// A port's carrier disagrees with its role.
    PortRole {
        /// The offending port.
        port: String,
    },
    /// A field metric does not have one entry per row and column.
    MetricShape {
        /// Number of rows and columns.
        width: usize,
        /// Number of supplied coefficients.
        entries: usize,
    },
    /// A field metric coefficient is not finite.
    MetricCoefficient {
        /// Row-major coefficient index.
        index: usize,
        /// Rejected coefficient.
        value: f64,
    },
    /// Mirrored coefficients of a field metric differ.
    MetricAsymmetric {
        /// Row of the lower-triangular entry.
        row: usize,
        /// Column of the lower-triangular entry.
        column: usize,
    },
    /// A field metric has a nonpositive Cholesky pivot.
    MetricPivot {
        /// Diagonal index.
        index: usize,
        /// Rejected pivot before taking its square root.
        value: f64,
    },
    /// A field metric's Cholesky factorization is not representable.
    MetricFactorization {
        /// Factor row.
        row: usize,
        /// Factor column.
        column: usize,
        /// Nonfinite residual or factor.
        value: f64,
    },
    /// The port list or an input or output vector has the wrong length for
    /// the layout or field metric.
    VectorLength {
        /// Which list or vector has the wrong length.
        vector: &'static str,
        /// Required length.
        expected: usize,
        /// Length supplied.
        actual: usize,
    },
    /// A value that must be finite is not.
    NonFinite {
        /// The quantity.
        quantity: &'static str,
    },
    /// The step is not positive and finite.
    Step {
        /// The rejected step [s].
        dt_s: f64,
    },
    /// Physical dissipation or entropy generation is negative.
    NegativeProduction {
        /// "dissipated_j" or "entropy_generated_j_per_k".
        quantity: &'static str,
        /// The reported value.
        value: f64,
    },
    /// Port energy does not equal stored change plus numerical loss.
    EnergyBalance {
        /// Net energy into the component through its ports [J].
        port_energy_in_j: f64,
        /// Stored change from the receipt [J].
        stored_change_j: f64,
        /// Numerical loss from the receipt [J].
        numerical_loss_j: f64,
        /// Tolerance the receipt claimed [J].
        tolerance_j: f64,
    },
    /// The receipt's stored change disagrees with the re-measured endpoints.
    EndpointEnergy {
        /// Stored energy of the candidate minus that of the accepted state [J].
        measured_change_j: f64,
        /// Stored change from the receipt [J].
        stored_change_j: f64,
        /// Tolerance applied [J].
        tolerance_j: f64,
    },
    /// The owner's residual exceeds its own tolerance.
    Residual {
        /// Reported residual.
        residual: f64,
        /// Reported tolerance.
        tolerance: f64,
    },
}

impl fmt::Display for ContractViolation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicatePort { port } => write!(f, "port identifier `{port}` is used twice"),
            Self::PortRole { port } => {
                write!(
                    f,
                    "port `{port}` has a carrier that disagrees with its role"
                )
            }
            Self::MetricShape { width, entries } => write!(
                f,
                "field metric of width {width} has {entries} entries, expected a square matrix"
            ),
            Self::MetricCoefficient { index, value } => {
                write!(
                    f,
                    "field metric coefficient {index} = {value} is not finite"
                )
            }
            Self::MetricAsymmetric { row, column } => write!(
                f,
                "field metric entries ({row}, {column}) and ({column}, {row}) differ"
            ),
            Self::MetricPivot { index, value } => write!(
                f,
                "field metric Cholesky pivot {index} = {value} is not positive"
            ),
            Self::MetricFactorization { row, column, value } => write!(
                f,
                "field metric Cholesky factor ({row}, {column}) = {value} is not finite"
            ),
            Self::VectorLength {
                vector,
                expected,
                actual,
            } => write!(
                f,
                "{vector} has length {actual}, the port layout needs {expected}"
            ),
            Self::NonFinite { quantity } => write!(f, "{quantity} is not finite"),
            Self::Step { dt_s } => write!(f, "step {dt_s} s is not positive and finite"),
            Self::NegativeProduction { quantity, value } => {
                write!(f, "{quantity} = {value} is negative")
            }
            Self::EnergyBalance {
                port_energy_in_j,
                stored_change_j,
                numerical_loss_j,
                tolerance_j,
            } => write!(
                f,
                "first law: port energy in {port_energy_in_j} J != stored change \
                 {stored_change_j} J + numerical loss {numerical_loss_j} J within {tolerance_j} J"
            ),
            Self::EndpointEnergy {
                measured_change_j,
                stored_change_j,
                tolerance_j,
            } => write!(
                f,
                "endpoint energy: measured change {measured_change_j} J != receipt \
                 {stored_change_j} J within {tolerance_j} J"
            ),
            Self::Residual {
                residual,
                tolerance,
            } => write!(
                f,
                "residual {residual} exceeds the owner's tolerance {tolerance}"
            ),
        }
    }
}

impl std::error::Error for ContractViolation {}

/// Energy and entropy that crossed a component's ports in one step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepAudit {
    /// Net energy into the component through all ports [J].
    pub port_energy_in_j: f64,
    /// Net entropy into the component through its heat ports [J/K].
    pub port_entropy_in_j_per_k: f64,
    /// `port_energy_in_j - stored_change_j - numerical_loss_j` [J].
    pub balance_defect_j: f64,
    /// Measured endpoint change minus the receipt's stored change [J].
    pub endpoint_defect_j: f64,
}

/// Check one proposed step against the contract.
///
/// This is the single place where the first law of a step is tested: the
/// kernel calls it for every component, and so does the conformance kit, so
/// components carry no audit code of their own. `stored_old_j` and
/// `stored_new_j` are [`Component::stored_energy_j`] of the accepted and the
/// candidate state.
///
/// # Errors
/// The first [`ContractViolation`] found, in this order: step, vector
/// lengths, non-finite values, residual, negative production, energy balance,
/// endpoint energy.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn audit_step(
    ports: &[PortSpec],
    layout: &PortLayout,
    dt_s: f64,
    inputs: &[f64],
    outputs: &[f64],
    receipt: &Receipt,
    stored_old_j: f64,
    stored_new_j: f64,
) -> Result<StepAudit, ContractViolation> {
    if !(dt_s.is_finite() && dt_s > 0.0) {
        return Err(ContractViolation::Step { dt_s });
    }
    if ports.len() != layout.ports() {
        return Err(ContractViolation::VectorLength {
            vector: "ports",
            expected: layout.ports(),
            actual: ports.len(),
        });
    }
    if inputs.len() != layout.input_len() {
        return Err(ContractViolation::VectorLength {
            vector: "inputs",
            expected: layout.input_len(),
            actual: inputs.len(),
        });
    }
    if outputs.len() != layout.output_len() {
        return Err(ContractViolation::VectorLength {
            vector: "outputs",
            expected: layout.output_len(),
            actual: outputs.len(),
        });
    }
    let finite = [
        ("inputs", inputs.iter().all(|v| v.is_finite())),
        ("outputs", outputs.iter().all(|v| v.is_finite())),
        ("stored_change_j", receipt.stored_change_j.is_finite()),
        ("dissipated_j", receipt.dissipated_j.is_finite()),
        ("numerical_loss_j", receipt.numerical_loss_j.is_finite()),
        (
            "entropy_generated_j_per_k",
            receipt.entropy_generated_j_per_k.is_finite(),
        ),
        (
            "energy_tolerance_j",
            receipt.energy_tolerance_j.is_finite() && receipt.energy_tolerance_j >= 0.0,
        ),
        ("residual", receipt.residual.is_finite()),
        ("tolerance", receipt.tolerance.is_finite()),
        (
            "stored energy",
            stored_old_j.is_finite() && stored_new_j.is_finite(),
        ),
    ];
    for (quantity, ok) in finite {
        if !ok {
            return Err(ContractViolation::NonFinite { quantity });
        }
    }
    if receipt.residual > receipt.tolerance {
        return Err(ContractViolation::Residual {
            residual: receipt.residual,
            tolerance: receipt.tolerance,
        });
    }
    if receipt.dissipated_j < 0.0 {
        return Err(ContractViolation::NegativeProduction {
            quantity: "dissipated_j",
            value: receipt.dissipated_j,
        });
    }
    if receipt.entropy_generated_j_per_k < 0.0 {
        return Err(ContractViolation::NegativeProduction {
            quantity: "entropy_generated_j_per_k",
            value: receipt.entropy_generated_j_per_k,
        });
    }
    let (mut energy, mut entropy) = (0.0_f64, 0.0_f64);
    for (index, port) in ports.iter().enumerate() {
        let input = &inputs[layout.input_range(index)];
        let output = &outputs[layout.output_range(index)];
        energy += port.energy_in_j(dt_s, input, output);
        entropy += port.entropy_in_j_per_k(dt_s, input, output);
    }
    if !energy.is_finite() {
        return Err(ContractViolation::NonFinite {
            quantity: "port energy",
        });
    }
    if !entropy.is_finite() {
        return Err(ContractViolation::NonFinite {
            quantity: "port entropy",
        });
    }
    let balance_defect_j = energy - receipt.stored_change_j - receipt.numerical_loss_j;
    if !(balance_defect_j.abs() <= receipt.energy_tolerance_j) {
        return Err(ContractViolation::EnergyBalance {
            port_energy_in_j: energy,
            stored_change_j: receipt.stored_change_j,
            numerical_loss_j: receipt.numerical_loss_j,
            tolerance_j: receipt.energy_tolerance_j,
        });
    }
    // The endpoint difference cancels, so it is allowed the roundoff of the
    // two stored energies on top of the owner's tolerance.
    let measured_change_j = stored_new_j - stored_old_j;
    let endpoint_tolerance_j =
        receipt.energy_tolerance_j + 4.0 * f64::EPSILON * (stored_old_j.abs() + stored_new_j.abs());
    let endpoint_defect_j = measured_change_j - receipt.stored_change_j;
    if !(endpoint_defect_j.abs() <= endpoint_tolerance_j) {
        return Err(ContractViolation::EndpointEnergy {
            measured_change_j,
            stored_change_j: receipt.stored_change_j,
            tolerance_j: endpoint_tolerance_j,
        });
    }
    Ok(StepAudit {
        port_energy_in_j: energy,
        port_entropy_in_j_per_k: entropy,
        balance_defect_j,
        endpoint_defect_j,
    })
}
