//! A linear structure in mass-normalized damped modes as one composition-kernel Component.
//!
//! [`ModalStructure`] wraps the exact held-force step of
//! [`fs_modal::time::ModalAcousticTimeModel`]: every mode advances exactly under its held
//! generalized force, the work is that force times the displacement increment the owner added, and
//! the viscous heat is the owner's own integral of `c v(t)^2` along the exact trajectory. Nothing is
//! integrated here. It is the linear case of this crate's modal structures: no quartic term.
//!
//! Ports, in this order:
//!
//! - one power port per [`ModalAttachment`], receiving effort. Row `r` of an attachment carries
//!   weights `s_rk`, one per mode, that turn the row's effort `e_r` into a generalized force on
//!   mode `k` (`f_k = sum_r s_rk e_r`, [N/sqrt(kg)]) and mode `k`'s velocity into the row's flow. The
//!   port returns `sum_k s_rk (q1_k - q0_k) / dt`: the mean flow of the step, which is the
//!   displacement the owner pairs with a held force. A point force on a bar is a row of mode shapes
//!   at the point [1/sqrt(kg)]; a mouth's pressure on an air column is a row of `1/sqrt(M_k)`, with
//!   `M_k` mode `k`'s acoustic mass per unit volume flow [kg/m^4]; an attachment whose weights are the
//!   identity is the modal port itself.
//! - `heat`, a heat port receiving effort: the temperature [K] at which the modal damping releases
//!   its heat. It returns minus the mean heat rate of the step, `-D / dt` [W].
//! - With [`ModalStructure::publishing_ends`], one signal output per attachment after `heat`,
//!   `<id>.end`, of the attachment's width: row `r` returns `sum_k s_rk v1_k`, the flow at the end
//!   of the step.
//!
//! The end of a step. The owner's step is exact for a held force, not the midpoint rule, so the
//! flow does not end its step at twice the mean less its start. An attachment's port says so
//! ([`Component::end_of_step`]): it answers [`EndOfStep::Told`] with its end output when the ends
//! are published, and [`EndOfStep::Untold`] otherwise, so that a system joins the structure to a
//! load line only where the line is told the end.
//!
//! The port map's tangent ([`PortTangent`]). Each mode's step is affine in its held force `f_k`:
//! a unit force adds the owner's unit-force response to its displacement and velocity
//! ([`ModalAcousticTimeModel::unit_force_response_into`], the mobility the owner's held-force
//! solve inverts). The mean flows' rows are `sum_k s_rk s_ck G_q,k / dt` and the ends' rows
//! `sum_k s_rk s_ck G_v,k`, the same at every state. The heat row is minus the derivative of the
//! step's heat over `dt`: per mode, `d(W - dE)/df_k = dq_k + f_k G_q,k - v1_k G_v,k - w_k^2 q1_k
//! G_q,k`, from the step at the iterate (0 for a mode with no damping). That row is a difference
//! of terms the size of the work, and its relative error grows as the damping shrinks: it is
//! measured in `tests/modal_structure.rs`, and nothing beyond that measurement is claimed. The
//! temperature's column is zero. `vjp` is the transpose of the same pieces.
//!
//! Receipt: the stored change is the modal mechanical energy `sum 1/2 (v^2 + w^2 q^2)` measured at
//! both ends by the owner's law ([`ModalAcousticTimeModel::energy_j`]); the dissipation is the
//! owner's viscous heat; the owner's work-heat residual `W - D - (E1 - E0)` is booked as numerical
//! loss (the owner refuses a step whose residual exceeds its own allowance); the entropy generated
//! is `D / T`.
//!
//! No claim: the narrow-band observer transfers of the owner are not used (an observer is a probe of
//! its own), so a model carrying one is refused; the heat row is quadratic in the efforts, so the
//! class is [`Class::Smooth`]; one propose clones the owner's model and allocates its workspace and
//! the force bank, and one tangent call allocates the gains and takes each damped mode's step once
//! more; the owner's statements (exact for a held force, the Nyquist
//! admission of its nominal clock, its budgets) are unchanged by the wrapping. The port flows are
//! endpoint differences `q1 - q0`, while the owner books its work on the increment it added (which
//! it keeps because a subtraction of endpoints can erase a small work); the owner does not publish
//! that increment, so the energy tolerance carries `|f| (|q0| + |q1|)` per mode. Where the
//! displacement is far larger than its change in a step (a statically held load, a free coordinate
//! drifting under a constant force), that tolerance is wide against the step's work, and a wrapper
//! leak below it would show in the ledger's accepted-defect account, not in the audit.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_modal::time::{ModalAcousticTimeError, ModalAcousticWorkspace};
use fs_ports::component::{
    Class, Component, ContractViolation, EndOfStep, PortRole, PortSpec, PortTangent, Proposal,
    Receipt,
};
use fs_ports::{PortKind, StableId};

pub use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};

/// One group of rows through which other components act on the modes.
#[derive(Debug, Clone, PartialEq)]
pub struct ModalAttachment {
    /// Port identifier, unique within the structure; `heat` is taken.
    pub id: StableId,
    /// What the port carries.
    pub kind: PortKind,
    /// Number of rows.
    pub width: NonZeroUsize,
    /// Row-major weights, `width` rows of one entry per mode.
    pub weights: Vec<f64>,
}

/// Why a [`ModalStructure`] could not be built or refused a call.
#[derive(Debug, Clone, PartialEq)]
pub enum ModalStructureRefusal {
    /// The owner refused; its error is unchanged.
    Owner(ModalAcousticTimeError),
    /// The model carries an observer transfer on this mode.
    ObserverTransfer {
        /// Index of the mode.
        mode: usize,
    },
    /// An attachment's weights do not hold `width` rows of one entry per mode.
    Weights {
        /// Index of the attachment in the supplied order.
        attachment: usize,
        /// Entries expected.
        expected: usize,
        /// Entries supplied.
        found: usize,
    },
    /// An attachment carries a weight that is not finite.
    NonFiniteWeight {
        /// Index of the attachment in the supplied order.
        attachment: usize,
    },
    /// A port identifier was refused.
    Identifier(fs_ports::CoupleError),
    /// A port specification or the port list was refused.
    Port(ContractViolation),
    /// The state does not hold one entry per mode.
    States {
        /// Modes of the structure.
        expected: usize,
        /// Entries in the state.
        found: usize,
    },
    /// The heat port received a temperature that is not positive and finite.
    Temperature {
        /// The value received [K].
        kelvin: f64,
    },
    /// A loaded basis was asked for no mode.
    NoModes,
    /// An unloaded angular frequency, or the square of a loaded one, that is not positive and
    /// finite.
    Frequency {
        /// Index of the mode (in the loaded order for a loaded one).
        mode: usize,
        /// The value refused [rad/s, or rad^2/s^2 for a loaded square].
        value: f64,
    },
    /// A matrix does not hold one entry per pair of modes.
    MatrixShape {
        /// Entries expected.
        expected: usize,
        /// Entries supplied.
        found: usize,
    },
    /// A matrix entry is not finite.
    NonFiniteEntry {
        /// Row of the entry.
        row: usize,
        /// Its column.
        column: usize,
    },
    /// A matrix entry differs from its mirror entry (a signed zero equals its mirror).
    Asymmetric {
        /// Row of the entry below the diagonal.
        row: usize,
        /// Its column.
        column: usize,
    },
    /// fs-modal refused the loaded mass `M + M_a`: it is not positive definite.
    LoadedMass(fs_modal::ModalError),
    /// fs-la refused the dense symmetric eigenproblem of this many modes: beyond its work cap.
    LoadedSize(fs_la::eigen::JacobiEighAdmissionError),
    /// Cancellation was requested.
    Cancelled,
}

impl core::fmt::Display for ModalStructureRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Owner(error) => write!(f, "{error}"),
            Self::ObserverTransfer { mode } => write!(
                f,
                "mode {mode} carries an observer transfer; an observer is a probe of its own"
            ),
            Self::Weights {
                attachment,
                expected,
                found,
            } => write!(
                f,
                "attachment {attachment} needs {expected} weights (rows times modes) and has {found}"
            ),
            Self::NonFiniteWeight { attachment } => {
                write!(
                    f,
                    "attachment {attachment} carries a weight that is not finite"
                )
            }
            Self::Identifier(error) => write!(f, "port identifier refused: {error:?}"),
            Self::Port(violation) => write!(f, "port refused: {violation}"),
            Self::States { expected, found } => write!(
                f,
                "the structure has {expected} modes and the state holds {found}"
            ),
            Self::Temperature { kelvin } => write!(
                f,
                "the modal damping releases heat at {kelvin} K, which is not a temperature"
            ),
            Self::NoModes => write!(f, "a loaded basis needs at least one mode"),
            Self::Frequency { mode, value } => write!(
                f,
                "mode {mode} has the frequency value {value}, which is not positive and finite"
            ),
            Self::MatrixShape { expected, found } => write!(
                f,
                "the matrix needs {expected} entries (modes squared) and has {found}"
            ),
            Self::NonFiniteEntry { row, column } => {
                write!(f, "matrix entry ({row}, {column}) is not finite")
            }
            Self::Asymmetric { row, column } => {
                write!(f, "matrix entry ({row}, {column}) differs from its mirror")
            }
            Self::LoadedMass(error) => write!(f, "the loaded mass was refused: {error}"),
            Self::LoadedSize(error) => write!(f, "the loaded basis was refused: {error}"),
            Self::Cancelled => write!(f, "cancelled during the modal step"),
        }
    }
}

impl std::error::Error for ModalStructureRefusal {}

impl From<ModalAcousticTimeError> for ModalStructureRefusal {
    fn from(error: ModalAcousticTimeError) -> Self {
        Self::Owner(error)
    }
}

/// A linear structure in mass-normalized damped modes, stepped by its owner.
#[derive(Debug, Clone)]
pub struct ModalStructure {
    model: ModalAcousticTimeModel,
    /// All attachments' rows, stacked in port order, row-major over modes.
    weights: Vec<f64>,
    rows: usize,
    ports: Vec<PortSpec>,
    /// Each attachment's width, in port order: ports `0..widths.len()` are theirs.
    widths: Vec<NonZeroUsize>,
    /// Whether each attachment's end flow is returned on a signal output after `heat`.
    tells_ends: bool,
}

/// The step's gains in each mode's held force, at one iterate.
struct Gains {
    /// Per mode: the displacement and velocity a unit held force adds to the step.
    unit: Vec<ModalAcousticState>,
    /// Per mode: the derivative of the step's heat in its held force [J sqrt(kg)/N].
    heat: Vec<f64>,
}

impl ModalStructure {
    /// Wrap an admitted owner model with the given attachments. The model's own state is not used:
    /// the system owns the state.
    ///
    /// # Errors
    /// An observer transfer on any mode, weights of the wrong count or not finite, or a refused
    /// port list (repeated identifiers, including `heat`).
    pub fn new(
        model: ModalAcousticTimeModel,
        attachments: Vec<ModalAttachment>,
    ) -> Result<Self, ModalStructureRefusal> {
        let modes = model.modes().len();
        if let Some(mode) = model
            .modes()
            .iter()
            .position(|mode| mode.pressure_per_modal_velocity != fs_math::c64::C64::ZERO)
        {
            return Err(ModalStructureRefusal::ObserverTransfer { mode });
        }
        let mut weights = Vec::new();
        let mut ports = Vec::with_capacity(attachments.len() + 1);
        let mut widths = Vec::with_capacity(attachments.len());
        let mut rows = 0_usize;
        for (index, attachment) in attachments.into_iter().enumerate() {
            let expected = attachment.width.get() * modes;
            if attachment.weights.len() != expected {
                return Err(ModalStructureRefusal::Weights {
                    attachment: index,
                    expected,
                    found: attachment.weights.len(),
                });
            }
            if attachment.weights.iter().any(|w| !w.is_finite()) {
                return Err(ModalStructureRefusal::NonFiniteWeight { attachment: index });
            }
            weights.extend_from_slice(&attachment.weights);
            rows += attachment.width.get();
            widths.push(attachment.width);
            ports.push(
                PortSpec::power(
                    attachment.id,
                    attachment.kind,
                    PortRole::EffortIn,
                    attachment.width,
                )
                .map_err(ModalStructureRefusal::Port)?,
            );
        }
        let heat = StableId::new("heat").map_err(ModalStructureRefusal::Identifier)?;
        ports.push(PortSpec::heat(heat, PortRole::EffortIn).map_err(ModalStructureRefusal::Port)?);
        fs_ports::component::PortLayout::new(&ports).map_err(ModalStructureRefusal::Port)?;
        Ok(Self {
            model,
            weights,
            rows,
            ports,
            widths,
            tells_ends: false,
        })
    }

    /// The same structure telling the end of each attachment's flow: one signal output per
    /// attachment after `heat`, `<id>.end`, of the attachment's width, whose row `r` returns
    /// `sum_k s_rk v1_k`. An attachment's port then answers [`EndOfStep::Told`] with that output.
    /// Asking twice adds them once.
    ///
    /// # Errors
    /// [`ModalStructureRefusal::Port`] when an end's name repeats a port of the structure (an
    /// attachment named `a.end` beside one named `a`).
    pub fn publishing_ends(mut self) -> Result<Self, ModalStructureRefusal> {
        if self.tells_ends {
            return Ok(self);
        }
        for (attachment, &width) in self.widths.iter().enumerate() {
            let id = StableId::new(format!("{}.end", self.ports[attachment].id().as_str()))
                .map_err(ModalStructureRefusal::Identifier)?;
            self.ports.push(PortSpec::signal_out(id, width));
        }
        fs_ports::component::PortLayout::new(&self.ports).map_err(ModalStructureRefusal::Port)?;
        self.tells_ends = true;
        Ok(self)
    }

    /// Admit the owner's model for `modes`, each `(angular frequency [rad/s], viscous damping
    /// ratio)`, with no observer transfer, all at rest, and wrap it. `sample_rate_hz` is the owner's
    /// nominal clock: it admits only modes below its Nyquist guard and caches the full- and
    /// half-sample steps; a step of any other duration is exact as well.
    ///
    /// # Errors
    /// The owner's admission refusals, and every refusal of [`Self::new`].
    pub fn from_modes(
        sample_rate_hz: u32,
        modes: &[(f64, f64)],
        budget: ModalAcousticTimeBudget,
        attachments: Vec<ModalAttachment>,
    ) -> Result<Self, ModalStructureRefusal> {
        let modes = modes
            .iter()
            .map(
                |&(angular_frequency_rad_s, damping_ratio)| ModalAcousticMode {
                    angular_frequency_rad_s,
                    damping_ratio,
                    pressure_per_modal_velocity: fs_math::c64::C64::ZERO,
                },
            )
            .collect();
        Self::new(
            ModalAcousticTimeModel::try_new(sample_rate_hz, modes, budget)?,
            attachments,
        )
    }

    /// Number of modes.
    #[must_use]
    pub fn modes(&self) -> usize {
        self.model.modes().len()
    }

    /// Every mode at rest: the state of a structure nobody has touched.
    #[must_use]
    pub fn at_rest(&self) -> Vec<ModalAcousticState> {
        vec![ModalAcousticState::default(); self.modes()]
    }

    /// The owner's model this structure steps.
    #[must_use]
    pub fn model(&self) -> &ModalAcousticTimeModel {
        &self.model
    }

    fn check_states(&self, states: &[ModalAcousticState]) -> Result<(), ModalStructureRefusal> {
        if states.len() == self.modes() {
            Ok(())
        } else {
            Err(ModalStructureRefusal::States {
                expected: self.modes(),
                found: states.len(),
            })
        }
    }

    /// Each mode's generalized force under the rows' efforts: `f_k = sum_r s_rk e_r`.
    fn forces(&self, efforts: &[f64]) -> Vec<f64> {
        let modes = self.modes();
        let mut force = vec![0.0_f64; modes];
        for (row, &effort) in efforts[..self.rows].iter().enumerate() {
            let weights = &self.weights[row * modes..(row + 1) * modes];
            for (f, &w) in force.iter_mut().zip(weights) {
                *f += w * effort;
            }
        }
        force
    }

    /// The step's gains at the iterate `inputs` from `old` over `dt_s`.
    fn gains(
        &self,
        old: &[ModalAcousticState],
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Gains, ModalStructureRefusal> {
        cx.checkpoint()
            .map_err(|_| ModalStructureRefusal::Cancelled)?;
        self.check_states(old)?;
        let kelvin = inputs[self.rows];
        if !(kelvin.is_finite() && kelvin > 0.0) {
            return Err(ModalStructureRefusal::Temperature { kelvin });
        }
        let modes = self.modes();
        let mut unit = vec![ModalAcousticState::default(); modes];
        self.model.unit_force_response_into(dt_s, &mut unit)?;
        let force = self.forces(inputs);
        let mut heat = vec![0.0_f64; modes];
        for k in 0..modes {
            if self.model.damping_rate_per_s(k) == 0.0 {
                continue;
            }
            // The step itself, without the owner's budgets: those are the proposal's to apply.
            let end = self
                .model
                .advance_exact_zoh_duration_step(k, old[k], force[k], dt_s)?;
            let omega = self.model.modes()[k].angular_frequency_rad_s;
            let (q1, v1) = (end.displacement_m_sqrt_kg, end.velocity_m_sqrt_kg_per_s);
            let (gq, gv) = (
                unit[k].displacement_m_sqrt_kg,
                unit[k].velocity_m_sqrt_kg_per_s,
            );
            // `W = f (q1 - q0)` and `dE = (v1^2 + w^2 q1^2) / 2 - E0`; the exact step's heat is
            // `W - dE`, and each of `q1` and `v1` moves by its gain per unit force.
            let moved = q1 - old[k].displacement_m_sqrt_kg;
            heat[k] = moved + force[k] * gq - v1 * gv - omega * omega * q1 * gq;
        }
        Ok(Gains { unit, heat })
    }
}

impl Component for ModalStructure {
    type State = Vec<ModalAcousticState>;
    type Error = ModalStructureRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(
        &self,
        state: &Vec<ModalAcousticState>,
    ) -> Result<f64, ModalStructureRefusal> {
        self.check_states(state)?;
        Ok(self.model.energy_j(state)?)
    }

    fn end_of_step(&self, port: usize) -> EndOfStep {
        if self.tells_ends && port < self.widths.len() {
            EndOfStep::Told(self.widths.len() + 1 + port)
        } else {
            EndOfStep::Untold
        }
    }

    fn propose(
        &self,
        old: &Vec<ModalAcousticState>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Vec<ModalAcousticState>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ModalStructureRefusal> {
        cx.checkpoint()
            .map_err(|_| ModalStructureRefusal::Cancelled)?;
        self.check_states(old)?;
        let (modes, rows) = (self.modes(), self.rows);
        let kelvin = inputs[rows];
        if !(kelvin.is_finite() && kelvin > 0.0) {
            return Err(ModalStructureRefusal::Temperature { kelvin });
        }
        let force = self.forces(inputs);

        let mut model = self.model.clone();
        model.restore_states(old)?;
        let mut workspace = ModalAcousticWorkspace::new(&model);
        let frame = model.step_duration_into(&force, dt_s, &mut workspace)?;
        // The owner's heat allowance is not part of the tolerance below: the same `D` enters the
        // audit through the heat port and through the numerical loss, so its error cancels.
        let (work_j, heat_j) = (frame.input_work_j, frame.viscous_dissipation_j);
        new.clear();
        new.extend_from_slice(model.states());

        // What the owner's tolerance has to cover, term by term:
        // - the kernel's f64 port products against the owner's work. The owner books
        //   `f_k * increment_k`; the port pairs each row's effort with a sum over modes of
        //   `s_rk (q1_k - q0_k) / dt`, so the sums are reordered and each endpoint difference
        //   rounds at the size of the endpoints, not of the increment;
        // - the heat port's `(-D / dt) * dt` against the owner's `D`;
        // - the subtractions that form the stored change and the numerical loss.
        // The endpoints are re-measured by the very function that forms the stored change, so the
        // endpoint check needs nothing beyond the roundoff the audit itself adds for them.
        let mut paired_j = 0.0_f64;
        for row in 0..rows {
            let weights = &self.weights[row * modes..(row + 1) * modes];
            let (mut flow, mut size) = (0.0_f64, 0.0_f64);
            for k in 0..modes {
                let moved = new[k].displacement_m_sqrt_kg - old[k].displacement_m_sqrt_kg;
                flow += weights[k] * moved;
                size += (weights[k] * moved).abs();
            }
            outputs[row] = flow / dt_s;
            paired_j += inputs[row].abs() * size;
        }
        outputs[rows] = -heat_j / dt_s;
        if self.tells_ends {
            for row in 0..rows {
                let weights = &self.weights[row * modes..(row + 1) * modes];
                let mut end = 0.0_f64;
                for (w, reached) in weights.iter().zip(new.iter()) {
                    end += w * reached.velocity_m_sqrt_kg_per_s;
                }
                outputs[rows + 1 + row] = end;
            }
        }
        let mut endpoints_j = 0.0_f64;
        for k in 0..modes {
            endpoints_j += force[k].abs()
                * (old[k].displacement_m_sqrt_kg.abs() + new[k].displacement_m_sqrt_kg.abs());
        }
        let before_j = self.model.energy_j(old)?;
        let after_j = self.model.energy_j(new)?;
        let stored_change_j = after_j - before_j;
        let scale_j = paired_j + endpoints_j + work_j.abs() + heat_j + stored_change_j.abs();
        // A dot product of n terms rounds by up to about n units in the last place of the sum of
        // the terms' sizes; every sum here has at most `rows + modes` terms.
        let terms = (rows + modes) as f64;
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j: heat_j,
            // `W - D - (E1 - E0)`: positive is energy the step lost that no law accounts for.
            numerical_loss_j: work_j - heat_j - stored_change_j,
            entropy_generated_j_per_k: heat_j / kelvin,
            energy_tolerance_j: 4.0 * (terms + 16.0) * f64::EPSILON * scale_j + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for ModalStructure {
    fn jvp(
        &self,
        old: &Vec<ModalAcousticState>,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ModalStructureRefusal> {
        let gains = self.gains(old, dt_s, inputs, cx)?;
        let (modes, rows) = (self.modes(), self.rows);
        let d_force = self.forces(d_inputs);
        let mut heat = 0.0_f64;
        for (slope, d) in gains.heat.iter().zip(&d_force) {
            heat += slope * d;
        }
        for row in 0..rows {
            let weights = &self.weights[row * modes..(row + 1) * modes];
            let (mut mean, mut end) = (0.0_f64, 0.0_f64);
            for ((w, unit), d) in weights.iter().zip(&gains.unit).zip(&d_force) {
                mean += w * (unit.displacement_m_sqrt_kg * d);
                end += w * (unit.velocity_m_sqrt_kg_per_s * d);
            }
            d_outputs[row] = mean / dt_s;
            if self.tells_ends {
                d_outputs[rows + 1 + row] = end;
            }
        }
        d_outputs[rows] = -heat / dt_s;
        Ok(())
    }

    fn vjp(
        &self,
        old: &Vec<ModalAcousticState>,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ModalStructureRefusal> {
        let gains = self.gains(old, dt_s, inputs, cx)?;
        let (modes, rows) = (self.modes(), self.rows);
        let mut bar_force = vec![0.0_f64; modes];
        for (k, bar) in bar_force.iter_mut().enumerate() {
            let (mut mean, mut end) = (0.0_f64, 0.0_f64);
            for row in 0..rows {
                let weight = self.weights[row * modes + k];
                mean += weight * bar_outputs[row];
                if self.tells_ends {
                    end += weight * bar_outputs[rows + 1 + row];
                }
            }
            *bar = gains.unit[k].displacement_m_sqrt_kg * mean / dt_s
                - gains.heat[k] * bar_outputs[rows] / dt_s
                + gains.unit[k].velocity_m_sqrt_kg_per_s * end;
        }
        for (row, bar_input) in bar_inputs[..rows].iter_mut().enumerate() {
            let weights = &self.weights[row * modes..(row + 1) * modes];
            let mut sum = 0.0_f64;
            for (w, b) in weights.iter().zip(&bar_force) {
                sum += w * b;
            }
            *bar_input = sum;
        }
        bar_inputs[rows] = 0.0;
        Ok(())
    }
}
