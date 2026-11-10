//! Ideal lumped elements for one-dimensional mechanical, rotational,
//! electrical and thermal networks.
//!
//! Each element is a [`Component`] with a closed energy account and the
//! tangent of its port map, so a device made of masses, springs, dampers,
//! heat capacities and conductances can be written down as data and stepped
//! by [`crate::System`] with no element-specific code in the driver.
//! Parameters are constants in SI units. Laws that depend on the state of a
//! material or come from material cards belong to the physics owners.
//!
//! # Domains
//!
//! [`Mass`], [`Spring`], [`Damper`] and [`Force`] are written for a point on
//! a line, and their `in_domain` constructors give the same four laws in any
//! other power domain of [`PortKind`], where effort times flow is power:
//!
//! | element | translation | rotation | electrical (series loop) |
//! |---|---|---|---|
//! | [`Mass`] | mass [kg] | inertia [kg m^2] | inductance [H] |
//! | [`Spring`] | stiffness [N/m] | stiffness [N m/rad] | 1 / capacitance [1/F] |
//! | [`Damper`] | damping [N s/m] | damping [N m s/rad] | resistance [ohm] |
//! | [`Force`] | force [N] | torque [N m] | source voltage [V] |
//!
//! The flow is velocity, angular velocity or current; a spring's state is
//! the displacement, twist or charge. A [`Gyrator`] joins two domains without
//! loss: it is the ideal DC machine between a current loop and a rotor.
//!
//! # Signs
//!
//! Every port value is positive into the element that owns it (the kernel's
//! rule). A [`Mass`] is joined to what it carries by a common-flow junction:
//! it returns its velocity and receives the sum of the members' efforts as
//! its net force, and each member receives minus the mass's velocity. A
//! two-port element ([`Spring`], [`Damper`]) whose port `a` is on mass A and
//! port `b` on mass B therefore sees its extension `x_B - x_A` grow at
//! `f_a - f_b`, and returns the force it exerts on A at `a` and on B at `b`.
//! A port held still from outside receives zero flow.
//!
//! Heat ports carry temperature [K] and heat rate [W]. A [`HeatCapacity`]
//! returns its temperature and receives the heat rate into it; an element
//! that sits at a temperature (damper, conductance, heater) receives the
//! temperature and returns the heat rate into itself, which is negative when
//! heat leaves it.
//!
//! # Integration
//!
//! Every element is the trapezoidal (implicit midpoint) rule on its own state
//! under the step-mean port values the kernel holds: a mass returns the mean
//! of its old and new velocity, a heat capacity the mean of its old and new
//! temperature, a spring the discrete gradient of its potential. Joined in a
//! loop, the kernel solves the interface equations, so the composed device
//! is integrated by the same rule as one monolithic system: second order,
//! A-stable, and energy-exact for the lossless parts.

use core::f64::consts::PI;
use core::fmt;
use core::num::NonZeroUsize;

use fs_exec::{Cancelled, Cx};
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};
use fs_ports::{PortKind, StableId};

use crate::ledger::Sum;

/// Why an element was not built or refused a step.
#[derive(Debug, Clone, PartialEq)]
pub enum ElementError {
    /// A parameter is not finite or is outside its admitted range.
    Parameter {
        /// Parameter name, with its unit.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },
    /// A waveform table is empty or its times are not finite and increasing.
    Table {
        /// What is wrong with it.
        what: &'static str,
    },
    /// A temperature received or reached is not positive and finite.
    Temperature {
        /// The rejected temperature [K].
        kelvin: f64,
    },
    /// A flow or effort received is not finite.
    Input {
        /// The rejected value.
        value: f64,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl fmt::Display for ElementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parameter { name, value } => write!(f, "parameter {name} = {value} is not admitted"),
            Self::Table { what } => write!(f, "waveform table: {what}"),
            Self::Temperature { kelvin } => {
                write!(f, "temperature {kelvin} K is not positive and finite")
            }
            Self::Input { value } => write!(f, "input {value} is not finite"),
            Self::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for ElementError {}

impl From<Cancelled> for ElementError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

fn port(name: &str, kind: Option<PortKind>, role: PortRole) -> PortSpec {
    let id = StableId::new(name).expect("element port names are valid stable ids");
    match kind {
        Some(kind) => PortSpec::power(id, kind, role, NonZeroUsize::MIN),
        None => PortSpec::heat(id, role),
    }
    .expect("element ports are valid power or heat ports")
}

fn heat(name: &str, role: PortRole) -> PortSpec {
    port(name, None, role)
}

fn positive(name: &'static str, value: f64) -> Result<f64, ElementError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(ElementError::Parameter { name, value })
    }
}

fn nonnegative(name: &'static str, value: f64) -> Result<f64, ElementError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(ElementError::Parameter { name, value })
    }
}

fn finite(name: &'static str, value: f64) -> Result<f64, ElementError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ElementError::Parameter { name, value })
    }
}

fn inputs_finite(inputs: &[f64]) -> Result<(), ElementError> {
    match inputs.iter().find(|value| !value.is_finite()) {
        Some(&value) => Err(ElementError::Input { value }),
        None => Ok(()),
    }
}

fn temperature(kelvin: f64) -> Result<f64, ElementError> {
    if kelvin.is_finite() && kelvin > 0.0 {
        Ok(kelvin)
    } else {
        Err(ElementError::Temperature { kelvin })
    }
}

/// Roundoff allowance for an energy balance whose terms have size `scale`.
fn roundoff(scale: f64) -> f64 {
    64.0 * f64::EPSILON * scale + f64::MIN_POSITIVE
}

/// A value prescribed as a function of time. Built only by its constructors,
/// which refuse values that are not finite.
#[derive(Debug, Clone, PartialEq)]
pub struct Waveform(Shape);

#[derive(Debug, Clone, PartialEq)]
enum Shape {
    /// The same value at every time.
    Constant(f64),
    /// `offset + amplitude sin(2 pi frequency t + phase)`.
    Sine {
        amplitude: f64,
        frequency_hz: f64,
        phase_rad: f64,
        offset: f64,
    },
    /// Linear between `(time [s], value)` knots with increasing times, and
    /// the end values held before the first knot and after the last.
    Table(Vec<(f64, f64)>),
}

impl Waveform {
    /// A constant value.
    ///
    /// # Errors
    /// A value that is not finite.
    pub fn constant(value: f64) -> Result<Self, ElementError> {
        Ok(Self(Shape::Constant(finite("value", value)?)))
    }

    /// `offset + amplitude sin(2 pi frequency_hz t + phase_rad)`.
    ///
    /// # Errors
    /// A parameter that is not finite, or a negative frequency.
    pub fn sine(
        amplitude: f64,
        frequency_hz: f64,
        phase_rad: f64,
        offset: f64,
    ) -> Result<Self, ElementError> {
        Ok(Self(Shape::Sine {
            amplitude: finite("amplitude", amplitude)?,
            frequency_hz: nonnegative("frequency_hz", frequency_hz)?,
            phase_rad: finite("phase_rad", phase_rad)?,
            offset: finite("offset", offset)?,
        }))
    }

    /// Linear between `(time [s], value)` knots, with the end values held
    /// before the first knot and after the last.
    ///
    /// # Errors
    /// An empty table, a value or time that is not finite, or times that do
    /// not strictly increase.
    pub fn table(knots: Vec<(f64, f64)>) -> Result<Self, ElementError> {
        if knots.is_empty() {
            return Err(ElementError::Table { what: "it has no knots" });
        }
        if knots.iter().any(|(t, v)| !t.is_finite() || !v.is_finite()) {
            return Err(ElementError::Table {
                what: "a time or value is not finite",
            });
        }
        if knots.windows(2).any(|pair| !(pair[1].0 > pair[0].0)) {
            return Err(ElementError::Table {
                what: "times do not strictly increase",
            });
        }
        Ok(Self(Shape::Table(knots)))
    }

    /// The lowest value the waveform takes at any time.
    #[must_use]
    pub fn lowest(&self) -> f64 {
        match &self.0 {
            Shape::Constant(value) => *value,
            Shape::Sine {
                amplitude,
                frequency_hz,
                phase_rad,
                offset,
            } => {
                if *frequency_hz == 0.0 {
                    offset + amplitude * phase_rad.sin()
                } else {
                    offset - amplitude.abs()
                }
            }
            Shape::Table(knots) => knots.iter().map(|knot| knot.1).fold(f64::INFINITY, f64::min),
        }
    }

    /// The value at time `t_s`.
    #[must_use]
    pub fn value_at(&self, t_s: f64) -> f64 {
        match &self.0 {
            Shape::Constant(value) => *value,
            Shape::Sine {
                amplitude,
                frequency_hz,
                phase_rad,
                offset,
            } => offset + amplitude * (2.0 * PI * frequency_hz * t_s + phase_rad).sin(),
            Shape::Table(knots) => table_value(knots, t_s),
        }
    }

    /// The exact mean over `[t0_s, t0_s + dt_s]`; the value at `t0_s` when
    /// `dt_s` is not positive.
    #[must_use]
    pub fn mean(&self, t0_s: f64, dt_s: f64) -> f64 {
        if !(dt_s > 0.0) {
            return self.value_at(t0_s);
        }
        match &self.0 {
            Shape::Constant(value) => *value,
            Shape::Sine {
                amplitude,
                frequency_hz,
                phase_rad,
                offset,
            } => {
                // mean of sin(w t + p) over [t0, t0 + dt]
                //   = sin(w (t0 + dt / 2) + p) sin(w dt / 2) / (w dt / 2),
                // written without the cancelling difference of cosines.
                let omega = 2.0 * PI * frequency_hz;
                let half = 0.5 * omega * dt_s;
                let sinc = if half == 0.0 { 1.0 } else { half.sin() / half };
                offset + amplitude * (omega * (t0_s + 0.5 * dt_s) + phase_rad).sin() * sinc
            }
            Shape::Table(knots) => table_integral(knots, t0_s, t0_s + dt_s) / dt_s,
        }
    }
}

fn table_value(knots: &[(f64, f64)], t: f64) -> f64 {
    let (first, last) = (knots[0], knots[knots.len() - 1]);
    if t <= first.0 {
        return first.1;
    }
    if t >= last.0 {
        return last.1;
    }
    // First knot strictly after t; it exists and is not the first.
    let upper = knots.partition_point(|knot| knot.0 <= t);
    let ((t0, v0), (t1, v1)) = (knots[upper - 1], knots[upper]);
    v0 + (v1 - v0) * (t - t0) / (t1 - t0)
}

/// Integral of the table over `[start, end]`, piece by piece. Each piece is
/// linear, so the trapezoid of its end values is exact.
fn table_integral(knots: &[(f64, f64)], start: f64, end: f64) -> f64 {
    let piece = |a: f64, b: f64| -> f64 {
        let (lo, hi) = (a.max(start), b.min(end));
        if hi > lo {
            (hi - lo) * 0.5 * (table_value(knots, lo) + table_value(knots, hi))
        } else {
            0.0
        }
    };
    let (first, last) = (knots[0].0, knots[knots.len() - 1].0);
    let mut total = piece(f64::NEG_INFINITY, first) + piece(last, f64::INFINITY);
    // Segments that can overlap [start, end].
    let from = knots.partition_point(|knot| knot.0 <= start).saturating_sub(1);
    for pair in knots[from..].windows(2) {
        if pair[0].0 >= end {
            break;
        }
        total += piece(pair[0].0, pair[1].0);
    }
    total
}

/// Position and velocity of a [`Mass`]. In another domain the same two
/// fields hold its displacement and flow: angle [rad] and angular velocity
/// [rad/s], or charge passed [C] and current [A].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Motion {
    /// Position [m].
    pub position_m: f64,
    /// Velocity [m/s].
    pub velocity_m_s: f64,
}

/// A point mass on a line. Port `p` receives the net force [N] and returns
/// the mean velocity of the step [m/s].
#[derive(Debug, Clone)]
pub struct Mass {
    mass_kg: f64,
    ports: [PortSpec; 1],
}

impl Mass {
    /// A mass of `mass_kg`.
    ///
    /// # Errors
    /// A mass that is not positive and finite.
    pub fn new(mass_kg: f64) -> Result<Self, ElementError> {
        Self::in_domain(PortKind::MechanicalForceVelocity, mass_kg)
    }

    /// The same law in the power domain `kind`: `inertance` is momentum per
    /// unit flow there, a moment of inertia [kg m^2] for rotation or an
    /// inductance [H] for a current loop.
    ///
    /// # Errors
    /// An inertance that is not positive and finite.
    pub fn in_domain(kind: PortKind, inertance: f64) -> Result<Self, ElementError> {
        Ok(Self {
            mass_kg: positive("mass_kg", inertance)?,
            ports: [port("p", Some(kind), PortRole::EffortIn)],
        })
    }

    /// Mass [kg].
    #[must_use]
    pub fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
}

impl Component for Mass {
    type State = Motion;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Motion) -> Result<f64, ElementError> {
        Ok(0.5 * self.mass_kg * state.velocity_m_s * state.velocity_m_s)
    }

    fn propose(
        &self,
        old: &Motion,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Motion,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        let force_n = inputs[0];
        let velocity = old.velocity_m_s + dt_s * force_n / self.mass_kg;
        let mean = 0.5 * (old.velocity_m_s + velocity);
        *new = Motion {
            position_m: old.position_m + dt_s * mean,
            velocity_m_s: velocity,
        };
        outputs[0] = mean;
        let stored_change_j = self.mass_kg * mean * (velocity - old.velocity_m_s);
        let scale =
            self.stored_energy_j(old)? + self.stored_energy_j(new)? + (dt_s * force_n * mean).abs();
        Ok(Proposal::Step(Receipt::lossless(
            stored_change_j,
            roundoff(scale),
        )))
    }
}

impl PortTangent for Mass {
    fn jvp(
        &self,
        _old: &Motion,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        d_outputs[0] = 0.5 * dt_s / self.mass_kg * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        _old: &Motion,
        dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        bar_inputs[0] = 0.5 * dt_s / self.mass_kg * bar_outputs[0];
        Ok(())
    }
}

/// A spring between ports `a` and `b` with force
/// `stiffness e + cubic e^3` at extension `e` [m]. The state is the extension
/// from the unloaded length. Each port receives a flow [m/s] and returns the
/// force [N] the spring exerts at that end.
///
/// The force of a step is the discrete gradient of the potential between the
/// old and the new extension, so the stored change equals the port work for
/// any step, also with the cubic term.
#[derive(Debug, Clone)]
pub struct Spring {
    stiffness_n_m: f64,
    cubic_n_m3: f64,
    ports: [PortSpec; 2],
}

impl Spring {
    /// A spring with linear stiffness `stiffness_n_m` and cubic coefficient
    /// `cubic_n_m3` (positive hardens, negative softens; zero is linear).
    ///
    /// # Errors
    /// A stiffness that is not positive and finite, or a cubic coefficient
    /// that is not finite.
    pub fn new(stiffness_n_m: f64, cubic_n_m3: f64) -> Result<Self, ElementError> {
        Self::in_domain(PortKind::MechanicalForceVelocity, stiffness_n_m, cubic_n_m3)
    }

    /// The same law in the power domain `kind`: effort per unit displacement
    /// there, a torsional stiffness [N m/rad] for rotation or the reciprocal
    /// of a capacitance [1/F] in a current loop (the state is then the
    /// charge on the capacitor).
    ///
    /// # Errors
    /// As for [`Self::new`].
    pub fn in_domain(kind: PortKind, stiffness: f64, cubic: f64) -> Result<Self, ElementError> {
        Ok(Self {
            stiffness_n_m: positive("stiffness_n_m", stiffness)?,
            cubic_n_m3: finite("cubic_n_m3", cubic)?,
            ports: [
                port("a", Some(kind), PortRole::FlowIn),
                port("b", Some(kind), PortRole::FlowIn),
            ],
        })
    }

    fn potential(&self, e: f64) -> f64 {
        0.5 * self.stiffness_n_m * e * e + 0.25 * self.cubic_n_m3 * e * e * e * e
    }

    fn potential_size(&self, e: f64) -> f64 {
        0.5 * self.stiffness_n_m * e * e + 0.25 * self.cubic_n_m3.abs() * e * e * e * e
    }

    /// `(V(b) - V(a)) / (b - a)`, written without the division.
    fn discrete_force(&self, a: f64, b: f64) -> f64 {
        0.5 * self.stiffness_n_m * (a + b) + 0.25 * self.cubic_n_m3 * (a + b) * (a * a + b * b)
    }

    /// Derivative of [`Self::discrete_force`] with respect to `b`.
    fn discrete_slope(&self, a: f64, b: f64) -> f64 {
        0.5 * self.stiffness_n_m + 0.25 * self.cubic_n_m3 * (a * a + 2.0 * a * b + 3.0 * b * b)
    }
}

impl Component for Spring {
    type State = f64;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        if self.cubic_n_m3 == 0.0 {
            Class::ExactLinear
        } else {
            Class::Smooth
        }
    }

    fn stored_energy_j(&self, extension_m: &f64) -> Result<f64, ElementError> {
        Ok(self.potential(*extension_m))
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        let (a, b) = (*old, old + dt_s * (inputs[0] - inputs[1]));
        *new = b;
        let force_n = self.discrete_force(a, b);
        outputs[0] = force_n;
        outputs[1] = -force_n;
        let stored_change_j = force_n * (b - a);
        let scale = self.potential_size(a) + self.potential_size(b) + stored_change_j.abs();
        Ok(Proposal::Step(Receipt::lossless(
            stored_change_j,
            roundoff(scale),
        )))
    }
}

impl PortTangent for Spring {
    fn jvp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let (a, b) = (*old, old + dt_s * (inputs[0] - inputs[1]));
        let d_force = self.discrete_slope(a, b) * dt_s * (d_inputs[0] - d_inputs[1]);
        d_outputs[0] = d_force;
        d_outputs[1] = -d_force;
        Ok(())
    }

    fn vjp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        // The Jacobian is s [[1, -1], [-1, 1]]: symmetric.
        self.jvp(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}

/// A linear viscous damper between ports `a` and `b`, force
/// `damping (f_a - f_b)`. The power it absorbs leaves as heat through port
/// `heat` at the temperature that port receives.
#[derive(Debug, Clone)]
pub struct Damper {
    damping_n_s_m: f64,
    ports: [PortSpec; 3],
}

impl Damper {
    /// A damper of `damping_n_s_m`.
    ///
    /// # Errors
    /// A damping coefficient that is not positive and finite.
    pub fn new(damping_n_s_m: f64) -> Result<Self, ElementError> {
        Self::in_domain(PortKind::MechanicalForceVelocity, damping_n_s_m)
    }

    /// The same law in the power domain `kind`: effort per unit flow there,
    /// a rotary damping [N m s/rad] or a resistance [ohm]. Its heat port is
    /// where the Joule or friction heat leaves.
    ///
    /// # Errors
    /// A coefficient that is not positive and finite.
    pub fn in_domain(kind: PortKind, resistance: f64) -> Result<Self, ElementError> {
        Ok(Self {
            damping_n_s_m: positive("damping_n_s_m", resistance)?,
            ports: [
                port("a", Some(kind), PortRole::FlowIn),
                port("b", Some(kind), PortRole::FlowIn),
                heat("heat", PortRole::EffortIn),
            ],
        })
    }
}

impl Component for Damper {
    type State = ();
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        let rate = inputs[0] - inputs[1];
        let kelvin = temperature(inputs[2])?;
        let force_n = self.damping_n_s_m * rate;
        let loss_w = force_n * rate;
        outputs[0] = force_n;
        outputs[1] = -force_n;
        outputs[2] = -loss_w;
        let dissipated_j = dt_s * loss_w;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dissipated_j / kelvin,
            energy_tolerance_j: roundoff(dt_s * (force_n * inputs[0]).abs()
                + dt_s * (force_n * inputs[1]).abs()
                + dissipated_j),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for Damper {
    fn jvp(
        &self,
        _old: &(),
        _dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let rate = inputs[0] - inputs[1];
        let d_rate = d_inputs[0] - d_inputs[1];
        d_outputs[0] = self.damping_n_s_m * d_rate;
        d_outputs[1] = -self.damping_n_s_m * d_rate;
        d_outputs[2] = -2.0 * self.damping_n_s_m * rate * d_rate;
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        _dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let rate = inputs[0] - inputs[1];
        let c = self.damping_n_s_m;
        let along = c * (bar_outputs[0] - bar_outputs[1]) - 2.0 * c * rate * bar_outputs[2];
        bar_inputs[0] = along;
        bar_inputs[1] = -along;
        bar_inputs[2] = 0.0;
        Ok(())
    }
}

/// Clock and energy account of a prescribed source ([`Force`], [`Heater`]).
///
/// The clock is summed with the same compensation as the kernel's own, so a
/// source started with the system reads the time a prescribed boundary
/// receives, bit for bit, however many steps have been taken.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SourceAccount {
    /// Time since the start [s].
    pub clock_s: Sum,
    /// Energy that has entered the source through its port [J]; negative when
    /// the source has delivered energy to the system.
    pub energy_in_j: f64,
}

impl SourceAccount {
    /// An account whose clock starts at `time_s`.
    #[must_use]
    pub fn starting_at(time_s: f64) -> Self {
        let mut clock_s = Sum::default();
        clock_s.add(time_s);
        Self {
            clock_s,
            energy_in_j: 0.0,
        }
    }

    /// Time on the source's clock [s].
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.clock_s.value()
    }
}

fn source_step(old: &SourceAccount, dt_s: f64, change_j: f64) -> (SourceAccount, Receipt) {
    let mut clock_s = old.clock_s;
    clock_s.add(dt_s);
    let new = SourceAccount {
        clock_s,
        energy_in_j: old.energy_in_j + change_j,
    };
    let scale = old.energy_in_j.abs() + new.energy_in_j.abs() + change_j.abs();
    (new, Receipt::lossless(change_j, roundoff(scale)))
}

/// A prescribed force [N] as a function of time, applied at port `p`. Its
/// stored energy is the account of the work it has exchanged, so a device
/// driven by it still closes its first law.
#[derive(Debug, Clone)]
pub struct Force {
    force_n: Waveform,
    ports: [PortSpec; 1],
}

impl Force {
    /// A force that follows `force_n`.
    #[must_use]
    pub fn new(force_n: Waveform) -> Self {
        Self::in_domain(PortKind::MechanicalForceVelocity, force_n)
    }

    /// The same source in the power domain `kind`: a prescribed torque
    /// [N m] or a source voltage [V].
    #[must_use]
    pub fn in_domain(kind: PortKind, effort: Waveform) -> Self {
        Self {
            force_n: effort,
            ports: [port("p", Some(kind), PortRole::FlowIn)],
        }
    }
}

impl Component for Force {
    type State = SourceAccount;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &SourceAccount) -> Result<f64, ElementError> {
        Ok(state.energy_in_j)
    }

    fn propose(
        &self,
        old: &SourceAccount,
        dt_s: f64,
        inputs: &[f64],
        new: &mut SourceAccount,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        let force_n = self.force_n.mean(old.time_s(), dt_s);
        outputs[0] = force_n;
        let (next, receipt) = source_step(old, dt_s, dt_s * force_n * inputs[0]);
        *new = next;
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for Force {
    fn jvp(
        &self,
        _old: &SourceAccount,
        _dt_s: f64,
        _inputs: &[f64],
        _d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        d_outputs[0] = 0.0;
        Ok(())
    }

    fn vjp(
        &self,
        _old: &SourceAccount,
        _dt_s: f64,
        _inputs: &[f64],
        _bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        bar_inputs[0] = 0.0;
        Ok(())
    }
}

/// A body of constant heat capacity [J/K]. Port `heat` receives the heat rate
/// into it [W] and returns its mean temperature over the step [K]. Its stored
/// energy is measured from absolute zero, `C T`.
#[derive(Debug, Clone)]
pub struct HeatCapacity {
    capacity_j_k: f64,
    ports: [PortSpec; 1],
}

impl HeatCapacity {
    /// A body of `capacity_j_k`.
    ///
    /// # Errors
    /// A capacity that is not positive and finite.
    pub fn new(capacity_j_k: f64) -> Result<Self, ElementError> {
        Ok(Self {
            capacity_j_k: positive("capacity_j_k", capacity_j_k)?,
            ports: [heat("heat", PortRole::FlowIn)],
        })
    }
}

impl Component for HeatCapacity {
    type State = f64;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, temperature_k: &f64) -> Result<f64, ElementError> {
        Ok(self.capacity_j_k * temperature_k)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        let heat_j = dt_s * inputs[0];
        let next = old + heat_j / self.capacity_j_k;
        if !(next > 0.0) {
            // A constant capacity holds only above absolute zero; a step that
            // would cross it is too long for the heat it carries.
            return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
        }
        *new = next;
        outputs[0] = 0.5 * (old + next);
        let scale = self.capacity_j_k * (old.abs() + next.abs());
        Ok(Proposal::Step(Receipt::lossless(heat_j, roundoff(scale))))
    }
}

impl PortTangent for HeatCapacity {
    fn jvp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        d_outputs[0] = 0.5 * dt_s / self.capacity_j_k * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        bar_inputs[0] = 0.5 * dt_s / self.capacity_j_k * bar_outputs[0];
        Ok(())
    }
}

/// A thermal conductance [W/K] between heat ports `a` and `b`. Heat
/// `G (T_a - T_b)` enters at `a` and leaves at `b`; it stores nothing and
/// generates `G (T_a - T_b)^2 / (T_a T_b)` of entropy per second.
#[derive(Debug, Clone)]
pub struct Conductance {
    conductance_w_k: f64,
    ports: [PortSpec; 2],
}

impl Conductance {
    /// A conductance of `conductance_w_k`.
    ///
    /// # Errors
    /// A conductance that is not positive and finite.
    pub fn new(conductance_w_k: f64) -> Result<Self, ElementError> {
        Ok(Self {
            conductance_w_k: positive("conductance_w_k", conductance_w_k)?,
            ports: [heat("a", PortRole::EffortIn), heat("b", PortRole::EffortIn)],
        })
    }
}

impl Component for Conductance {
    type State = ();
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        let (hot, cold) = (temperature(inputs[0])?, temperature(inputs[1])?);
        let rate_w = self.conductance_w_k * (hot - cold);
        outputs[0] = rate_w;
        outputs[1] = -rate_w;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s
                * self.conductance_w_k
                * (hot - cold)
                * (hot - cold)
                / (hot * cold),
            energy_tolerance_j: roundoff(dt_s * rate_w.abs()),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for Conductance {
    fn jvp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let d_rate = self.conductance_w_k * (d_inputs[0] - d_inputs[1]);
        d_outputs[0] = d_rate;
        d_outputs[1] = -d_rate;
        Ok(())
    }

    fn vjp(
        &self,
        old: &(),
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        // The Jacobian is G [[1, -1], [-1, 1]]: symmetric.
        self.jvp(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}

/// A Joule heater: prescribed electrical power [W] turned into heat at the
/// temperature port `heat` receives. Its stored energy is the account of the
/// energy it has delivered (negative), and the conversion generates
/// `P / T` of entropy per second.
#[derive(Debug, Clone)]
pub struct Heater {
    power_w: Waveform,
    ports: [PortSpec; 1],
}

impl Heater {
    /// A heater that follows `power_w`.
    ///
    /// # Errors
    /// A power that is negative at some time: Joule heating only turns work
    /// into heat.
    pub fn new(power_w: Waveform) -> Result<Self, ElementError> {
        nonnegative("power_w", power_w.lowest())?;
        Ok(Self {
            power_w,
            ports: [heat("heat", PortRole::EffortIn)],
        })
    }
}

impl Component for Heater {
    type State = SourceAccount;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &SourceAccount) -> Result<f64, ElementError> {
        Ok(state.energy_in_j)
    }

    fn propose(
        &self,
        old: &SourceAccount,
        dt_s: f64,
        inputs: &[f64],
        new: &mut SourceAccount,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        let kelvin = temperature(inputs[0])?;
        // Non-negative by construction; the mean of a non-negative waveform
        // is non-negative.
        let power_w = self.power_w.mean(old.time_s(), dt_s).max(0.0);
        outputs[0] = -power_w;
        let (next, mut receipt) = source_step(old, dt_s, -dt_s * power_w);
        *new = next;
        receipt.dissipated_j = dt_s * power_w;
        receipt.entropy_generated_j_per_k = dt_s * power_w / kelvin;
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for Heater {
    fn jvp(
        &self,
        _old: &SourceAccount,
        _dt_s: f64,
        _inputs: &[f64],
        _d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        d_outputs[0] = 0.0;
        Ok(())
    }

    fn vjp(
        &self,
        _old: &SourceAccount,
        _dt_s: f64,
        _inputs: &[f64],
        _bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        bar_inputs[0] = 0.0;
        Ok(())
    }
}

/// A lossless gyrator between ports `a` and `b`, which may be in different
/// power domains: the effort at each port is `ratio` times the flow at the
/// other, with opposite signs, so the power entering at one port leaves at
/// the other.
///
/// It is the ideal DC machine. With its constant as `ratio` [N m/A, the same
/// number in V s/rad], port `a` in a current loop and port `b` on a rotor,
/// both joined as members of common-flow junctions, a positive loop current
/// drives the rotor with torque `ratio i` and a positive rotor speed puts
/// the back-EMF `ratio w` against the loop. Each port receives a flow and
/// returns an effort: `e_a = ratio f_b` and `e_b = -ratio f_a`.
#[derive(Debug, Clone)]
pub struct Gyrator {
    ratio: f64,
    ports: [PortSpec; 2],
}

impl Gyrator {
    /// A gyrator of `ratio` between the power domains `kind_a` and `kind_b`.
    ///
    /// # Errors
    /// A ratio that is zero or not finite.
    pub fn new(kind_a: PortKind, kind_b: PortKind, ratio: f64) -> Result<Self, ElementError> {
        if !(ratio.is_finite() && ratio != 0.0) {
            return Err(ElementError::Parameter {
                name: "ratio",
                value: ratio,
            });
        }
        Ok(Self {
            ratio,
            ports: [
                port("a", Some(kind_a), PortRole::FlowIn),
                port("b", Some(kind_b), PortRole::FlowIn),
            ],
        })
    }
}

impl Component for Gyrator {
    type State = ();
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        outputs[0] = self.ratio * inputs[1];
        outputs[1] = -self.ratio * inputs[0];
        // The two port powers are equal and opposite products of the same
        // three numbers; what is left is their rounding.
        Ok(Proposal::Step(Receipt::lossless(
            0.0,
            roundoff(2.0 * dt_s * (self.ratio * inputs[0] * inputs[1]).abs()),
        )))
    }
}

impl PortTangent for Gyrator {
    fn jvp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        d_outputs[0] = self.ratio * d_inputs[1];
        d_outputs[1] = -self.ratio * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        bar_inputs[0] = -self.ratio * bar_outputs[1];
        bar_inputs[1] = self.ratio * bar_outputs[0];
        Ok(())
    }
}
