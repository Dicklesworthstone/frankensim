//! Objective transport of retained two-axis constitutive history.
//!
//! This is a change of orientation, not a constitutive time step. The existing
//! contact-frame convention supplies the minimum rotation between normals,
//! followed by the caller's common material spin about the current normal.
//! Changing only the tangent gauge re-expresses the same physical history.
//! Relative torsional slip is a separate motion and must not be supplied here.
//!
//! Finite normal rotations follow the supplied path: composing intermediate
//! frames need not equal a single endpoint rotation. The mechanical owner must
//! resolve that path and common spin to its required temporal accuracy.

use super::{DryTangentialState, TangentialLaw};
use crate::partial_slip::{PartialSlipError, TangentFrame};
use crate::{DryInterfaceSystemCard, FrictionLaw, TriboError, finite, nonnegative_finite};
use core::fmt;

/// A checked proper isometry from the previous tangent coordinates to current
/// coordinates. Its construction reuses the finite-patch contact convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentialHistoryRotation {
    cosine_sine: [f64; 2],
}

impl TangentialHistoryRotation {
    /// Minimum normal rotation, then physical common spin, expressed in the
    /// current gauge. Antipodal normals require resolved intermediate geometry.
    pub fn between(
        previous: TangentFrame,
        current: TangentFrame,
        common_spin_angle_rad: f64,
    ) -> Result<Self, TangentialTransportError> {
        finite(common_spin_angle_rad, "common_material_spin_angle_rad")?;
        Ok(Self {
            cosine_sine: current.history_rotation_from(previous, common_spin_angle_rad)?,
        })
    }

    /// Cosine and sine of the coordinate isometry; neither is a slip increment.
    #[must_use]
    pub const fn cosine_sine(self) -> [f64; 2] {
        self.cosine_sine
    }

    /// Re-express one finite tangent vector without projection or renormalizing
    /// the vector. Work-conjugate vector components use the same isometry.
    pub fn apply(self, vector: [f64; 2]) -> Result<[f64; 2], TriboError> {
        for value in vector {
            finite(value, "transported_tangent_vector")?;
        }
        if self.cosine_sine == [1.0, 0.0] {
            return Ok(vector);
        }
        let [c, s] = self.cosine_sine;
        let result = [c * vector[0] - s * vector[1], s * vector[0] + c * vector[1]];
        for value in result {
            finite(value, "transported_tangent_vector")?;
        }
        Ok(result)
    }
}

/// Optional capability for laws whose retained history supports objective
/// two-axis rotation. Scalar memory laws deliberately do not implement it.
pub trait ObjectiveTangentialLaw: TangentialLaw {
    /// Return the same original history with every tangent vector transported.
    /// Preserve source ownership, nonvector state, parameters and all clocks;
    /// identity rotation must return an exactly equal state. No relative slip,
    /// constitutive relaxation, parameter continuation, heat or work occurs.
    /// The framed API checks finite nonnegative storage and its preservation.
    fn transport_state(
        &self,
        state: &Self::State,
        rotation: TangentialHistoryRotation,
    ) -> Result<Self::State, TriboError>;
}

/// Explicit numerical tolerance for a zero-work isometric transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentialTransportConfig {
    /// Nonnegative absolute roundoff allowance [J]; zero is admitted.
    pub absolute_energy_tolerance_j: f64,
    /// Relative roundoff allowance in `[0, 1e-8]`; zero is admitted.
    pub relative_energy_tolerance: f64,
}

impl TangentialTransportConfig {
    /// Both tolerances may be zero to require exact storage preservation.
    pub fn validate(self) -> Result<(), TangentialTransportError> {
        if !self.absolute_energy_tolerance_j.is_finite() || self.absolute_energy_tolerance_j < 0.0 {
            return Err(TangentialTransportError::Invalid(
                "absolute_energy_tolerance_j",
            ));
        }
        if !self.relative_energy_tolerance.is_finite()
            || !(0.0..=1e-8).contains(&self.relative_energy_tolerance)
        {
            return Err(TangentialTransportError::Invalid(
                "relative_energy_tolerance",
            ));
        }
        Ok(())
    }
}

/// Native accepted history tied to the frame in which its vector components
/// were accepted. The native constitutive owner retains source responsibility.
#[derive(Clone, Debug, PartialEq)]
pub struct FramedTangentialState<S> {
    state: S,
    frame: TangentFrame,
}

impl<S> FramedTangentialState<S> {
    /// Attach a frame to a caller-declared native checkpoint. This validates the
    /// original law/state; it cannot infer a missing historical frame.
    pub fn new<L: ObjectiveTangentialLaw<State = S>>(
        law: &L,
        state: S,
        frame: TangentFrame,
    ) -> Result<Self, TangentialTransportError>
    where
        S: PartialEq,
    {
        let identity = TangentialHistoryRotation {
            cosine_sine: [1.0, 0.0],
        };
        if law.transport_state(&state, identity)? != state {
            return Err(TangentialTransportError::Invalid(
                "identity_transport_changed_state",
            ));
        }
        nonnegative_finite(law.stored_energy_j(&state), "transport_stored_energy_j")?;
        Ok(Self { state, frame })
    }

    /// Original native checkpoint in the associated accepted frame.
    #[must_use]
    pub const fn state(&self) -> &S {
        &self.state
    }

    /// Frame associated with this checkpoint, not a caller's latest predictor.
    #[must_use]
    pub const fn frame(&self) -> TangentFrame {
        self.frame
    }

    /// Consume the wrapper when the coupled owner publishes a native step.
    #[must_use]
    pub fn into_state(self) -> S {
        self.state
    }

    /// Propose orientation transport without mutating the accepted checkpoint.
    /// The reported storage residual is numerical error, never frictional heat.
    /// Coupled owners must include it in their own total energy acceptance.
    pub fn transport<L: ObjectiveTangentialLaw<State = S>>(
        &self,
        law: &L,
        current: TangentFrame,
        common_spin_angle_rad: f64,
        config: TangentialTransportConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TangentialTransportStep<S>, TangentialTransportError> {
        if cancelled() {
            return Err(TangentialTransportError::Cancelled);
        }
        config.validate()?;
        let rotation =
            TangentialHistoryRotation::between(self.frame, current, common_spin_angle_rad)?;
        let before_j = law.stored_energy_j(&self.state);
        nonnegative_finite(before_j, "transport_old_energy_j")?;
        if cancelled() {
            return Err(TangentialTransportError::Cancelled);
        }
        let state = law.transport_state(&self.state, rotation)?;
        let after_j = law.stored_energy_j(&state);
        nonnegative_finite(after_j, "transport_new_energy_j")?;
        let stored_energy_residual_j = after_j - before_j;
        let tolerance_j = config.absolute_energy_tolerance_j
            + config.relative_energy_tolerance * before_j.max(after_j);
        finite(tolerance_j, "transport_energy_tolerance_j")?;
        if stored_energy_residual_j.abs() > tolerance_j {
            return Err(TangentialTransportError::EnergyDrift {
                before_j,
                after_j,
                tolerance_j,
            });
        }
        if cancelled() {
            return Err(TangentialTransportError::Cancelled);
        }
        Ok(TangentialTransportStep {
            checkpoint: Self {
                state,
                frame: current,
            },
            stored_energy_residual_j,
        })
    }
}

/// Candidate orientation transfer. It has no interval, work or heat port.
#[derive(Clone, Debug, PartialEq)]
pub struct TangentialTransportStep<S> {
    /// Same physical native history expressed in the new accepted frame.
    pub checkpoint: FramedTangentialState<S>,
    /// `H(transformed state) - H(original state)` [J], retained as roundoff.
    pub stored_energy_residual_j: f64,
}

/// Refusal leaves both the original native state and its frame unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum TangentialTransportError {
    /// Malformed explicit tolerance or unsupported identity behavior.
    Invalid(&'static str),
    /// The existing contact-frame convention refused the orientation change.
    Frame(PartialSlipError),
    /// Original constitutive owner, state or finite arithmetic refused.
    Law(TriboError),
    /// Caller cancellation before candidate publication.
    Cancelled,
    /// The caller's numerical storage-preservation budget was exceeded.
    EnergyDrift {
        /// Original reversible storage [J].
        before_j: f64,
        /// Candidate reversible storage [J].
        after_j: f64,
        /// Caller absolute plus relative allowance [J].
        tolerance_j: f64,
    },
}

impl From<TriboError> for TangentialTransportError {
    fn from(value: TriboError) -> Self {
        Self::Law(value)
    }
}
impl From<PartialSlipError> for TangentialTransportError {
    fn from(value: PartialSlipError) -> Self {
        Self::Frame(value)
    }
}
impl fmt::Display for TangentialTransportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(field) => write!(f, "invalid tangential transport {field}"),
            Self::Frame(error) => write!(f, "tangential frame transport: {error}"),
            Self::Law(error) => write!(f, "tangential history transport: {error}"),
            Self::Cancelled => f.write_str("tangential history transport cancelled"),
            Self::EnergyDrift {
                before_j,
                after_j,
                tolerance_j,
            } => write!(
                f,
                "tangential transport changed storage from {before_j} to {after_j} J (tolerance {tolerance_j} J)"
            ),
        }
    }
}
impl std::error::Error for TangentialTransportError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Frame(error) => Some(error),
            Self::Law(error) => Some(error),
            _ => None,
        }
    }
}

impl ObjectiveTangentialLaw for FrictionLaw {
    fn transport_state(
        &self,
        state: &DryTangentialState,
        rotation: TangentialHistoryRotation,
    ) -> Result<DryTangentialState, TriboError> {
        self.validate()?;
        if let Some(temperature) = state.contact_temperature_k {
            crate::positive_finite(temperature, "transport_contact_temperature_k")?;
        }
        Ok(DryTangentialState {
            slip_m: rotation.apply(state.slip_m)?,
            contact_temperature_k: state.contact_temperature_k,
        })
    }
}

impl ObjectiveTangentialLaw for DryInterfaceSystemCard {
    fn transport_state(
        &self,
        state: &DryTangentialState,
        rotation: TangentialHistoryRotation,
    ) -> Result<DryTangentialState, TriboError> {
        // Orientation changes neither the immutable source card nor its last
        // explicit thermal boundary; a later physical step rechecks its domain.
        self.friction_law().transport_state(state, rotation)
    }
}
