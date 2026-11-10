//! Source-owned nonlinear storage on the composition kernel.
//!
//! The electrical curves, their finite support, energies, discrete efforts
//! and derivatives belong to `fs-phs`. These adapters only advance the
//! conserved charge or flux from the port input and expose that same law as
//! a [`Component`]. Each finite-step effort is the curve's integral mean;
//! endpoint voltage or current is obtained from [`NonlinearCapacitor::curve`]
//! or [`NonlinearInductor::curve`]. They need not equal the integral mean.
//!
//! Curves retain their authored piecewise-linear knots and refuse
//! extrapolation. At a coincident knot, tangents use the owner's one-sided
//! Newton representative; they do not claim a unique derivative there.
//! No hysteresis, leakage, temperature evolution or material provenance is
//! inferred from an authored table. The `constitutive` feature enables this
//! module without adding `fs-phs` to the default kernel's dependency closure.

use core::fmt;
use core::num::NonZeroUsize;

use fs_exec::{Cancelled, Cx};
use fs_phs::circuit::CircuitError;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};
use fs_ports::{PortKind, StableId};

pub use fs_phs::circuit::capacitance::{CapacitorChargeCurve, MAX_CAPACITOR_KNOTS};
pub use fs_phs::circuit::inductance::{InductorFluxCurve, MAX_INDUCTOR_KNOTS};

/// A nonlinear storage adapter refused an input or an original curve query.
#[derive(Debug)]
pub enum ConstitutiveError {
    /// The original source-law refusal, unchanged.
    Curve(CircuitError),
    /// An invalid port shape, timestep or unrepresentable arithmetic.
    Invalid {
        /// The rejected quantity or operation.
        what: &'static str,
    },
    /// The caller requested cancellation.
    Cancelled,
}

impl fmt::Display for ConstitutiveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Curve(error) => write!(f, "constitutive curve refused: {error:?}"),
            Self::Invalid { what } => write!(f, "constitutive storage refused: {what}"),
            Self::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for ConstitutiveError {}

impl From<CircuitError> for ConstitutiveError {
    fn from(error: CircuitError) -> Self {
        Self::Curve(error)
    }
}

impl From<Cancelled> for ConstitutiveError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

fn invalid(what: &'static str) -> ConstitutiveError {
    ConstitutiveError::Invalid { what }
}

fn scalar(values: &[f64]) -> Result<f64, ConstitutiveError> {
    match values {
        [value] if value.is_finite() => Ok(*value),
        _ => Err(invalid("one finite scalar port value required")),
    }
}

fn endpoint(old: f64, dt_s: f64, inputs: &[f64]) -> Result<f64, ConstitutiveError> {
    if !(dt_s.is_finite() && dt_s > 0.0) {
        return Err(invalid("positive finite timestep required"));
    }
    let next = old + dt_s * scalar(inputs)?;
    if next.is_finite() {
        Ok(next)
    } else {
        Err(invalid("conserved coordinate is not finite"))
    }
}

fn port(role: PortRole) -> [PortSpec; 1] {
    [PortSpec::power(
        StableId::new("p").expect("static port name"),
        PortKind::ElectricalVoltageCurrent,
        role,
        NonZeroUsize::MIN,
    )
    .expect("electrical power port")]
}

/// A reversible capacitor using the original bounded voltage/charge law.
///
/// State is charge [C]. Port `p` receives current [A] into the capacitor and
/// returns the curve's work-conjugate mean voltage [V]. Its energy is the
/// exact integral of voltage with respect to charge, including crossed knots.
#[derive(Clone, Debug)]
pub struct NonlinearCapacitor {
    curve: CapacitorChargeCurve,
    ports: [PortSpec; 1],
}

/// A reversible inductor using the original bounded current/flux law.
///
/// State is flux linkage [Wb]. Port `p` receives voltage [V] and returns the
/// curve's work-conjugate mean current [A] into the inductor. Its energy is
/// the exact integral of current with respect to flux, including crossed knots.
#[derive(Clone, Debug)]
pub struct NonlinearInductor {
    curve: InductorFluxCurve,
    ports: [PortSpec; 1],
}

// Charge and flux have the same conservative update. The source-owner
// methods supplied here keep their different units and conjugate outputs.
macro_rules! storage_adapter {
    ($name:ident, $curve:ty, $role:ident, $mean:ident, $derivative:ident) => {
        impl $name {
            /// Wrap an already admitted original curve without resampling it.
            #[must_use]
            pub fn new(curve: $curve) -> Self {
                Self { curve, ports: port(PortRole::$role) }
            }

            /// The original curve, including endpoint values and source support.
            #[must_use]
            pub fn curve(&self) -> &$curve {
                &self.curve
            }

            fn action(
                &self,
                old: &f64,
                dt_s: f64,
                inputs: &[f64],
                direction: &[f64],
                output: &mut [f64],
                cx: &Cx<'_>,
            ) -> Result<(), ConstitutiveError> {
                cx.checkpoint()?;
                if output.len() != 1 {
                    return Err(invalid("one scalar tangent output required"));
                }
                let next = endpoint(*old, dt_s, inputs)?;
                let result = dt_s * self.curve.$derivative(*old, next)? * scalar(direction)?;
                if !result.is_finite() {
                    return Err(invalid("port tangent is not finite"));
                }
                cx.checkpoint()?;
                output[0] = result;
                Ok(())
            }
        }

        impl Component for $name {
            type State = f64;
            type Error = ConstitutiveError;

            fn ports(&self) -> &[PortSpec] {
                &self.ports
            }

            fn class(&self) -> Class {
                Class::Monotone1d
            }

            fn stored_energy_j(&self, state: &f64) -> Result<f64, ConstitutiveError> {
                Ok(self.curve.stored_energy_j(*state)?)
            }

            fn propose(
                &self,
                old: &f64,
                dt_s: f64,
                inputs: &[f64],
                new: &mut f64,
                outputs: &mut [f64],
                cx: &Cx<'_>,
            ) -> Result<Proposal, ConstitutiveError> {
                cx.checkpoint()?;
                if outputs.len() != 1 {
                    return Err(invalid("one scalar port output required"));
                }
                let next = endpoint(*old, dt_s, inputs)?;
                let conjugate = self.curve.$mean(*old, next)?;
                let change = self.curve.energy_change_j(*old, next)?;
                let initial_energy = self.curve.stored_energy_j(*old)?;
                let final_energy = self.curve.stored_energy_j(next)?;
                let rounding = 64.0 * f64::EPSILON;
                let allowance = rounding * initial_energy.abs()
                    + rounding * final_energy.abs()
                    + rounding * change.abs()
                    + f64::MIN_POSITIVE;
                cx.checkpoint()?;
                *new = next;
                outputs[0] = conjugate;
                Ok(Proposal::Step(Receipt::lossless(change, allowance)))
            }
        }

        impl PortTangent for $name {
            fn jvp(
                &self,
                old: &f64,
                dt_s: f64,
                inputs: &[f64],
                d_inputs: &[f64],
                d_outputs: &mut [f64],
                cx: &Cx<'_>,
            ) -> Result<(), ConstitutiveError> {
                self.action(old, dt_s, inputs, d_inputs, d_outputs, cx)
            }

            fn vjp(
                &self,
                old: &f64,
                dt_s: f64,
                inputs: &[f64],
                bar_outputs: &[f64],
                bar_inputs: &mut [f64],
                cx: &Cx<'_>,
            ) -> Result<(), ConstitutiveError> {
                // The input/output Jacobian is one scalar, equal to its transpose.
                self.action(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
            }
        }
    };
}

storage_adapter!(
    NonlinearCapacitor,
    CapacitorChargeCurve,
    FlowIn,
    discrete_voltage_v,
    discrete_voltage_end_derivative
);
storage_adapter!(
    NonlinearInductor,
    InductorFluxCurve,
    EffortIn,
    discrete_current_a,
    discrete_current_end_derivative
);
