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
//!
//! Receipt: the stored change is the modal mechanical energy `sum 1/2 (v^2 + w^2 q^2)` measured at
//! both ends by the owner's law ([`ModalAcousticTimeModel::energy_j`]); the dissipation is the
//! owner's viscous heat; the owner's work-heat residual `W - D - (E1 - E0)` is booked as numerical
//! loss (the owner refuses a step whose residual exceeds its own allowance); the entropy generated
//! is `D / T`.
//!
//! No claim: the narrow-band observer transfers of the owner are not used (an observer is a probe of
//! its own), so a model carrying one is refused; the port map's tangent is not offered, so a kernel
//! loop differences it, although the mechanical rows are affine in the efforts (the heat row is
//! quadratic, so the class is [`Class::Smooth`]); one propose clones the owner's model and allocates
//! its workspace and the force bank; the owner's statements (exact for a held force, the Nyquist
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
    Class, Component, ContractViolation, PortRole, PortSpec, Proposal, Receipt,
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
        })
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
        let mut force = vec![0.0_f64; modes];
        for (row, &effort) in inputs[..rows].iter().enumerate() {
            let weights = &self.weights[row * modes..(row + 1) * modes];
            for (f, &w) in force.iter_mut().zip(weights) {
                *f += w * effort;
            }
        }

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
