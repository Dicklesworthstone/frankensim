//! An ideal dry-friction slider between two ends, with constant static and
//! kinetic coefficients under a constant normal force.
//!
//! It is the element the non-smooth solve tier is exercised with, and
//! nothing more: friction that follows a material, a surface card, a
//! temperature or a slip history belongs to `fs-tribo`, whose law a
//! component of its own puts behind the same tier.
//!
//! # The law of a step
//!
//! The elements of this crate use the midpoint rule: a mass returns its
//! MEAN velocity over the step on its port, and its velocity at the end of
//! the step is twice that mean less the velocity it started with. The
//! slider keeps the slip it ended the last step with, `s0`, and is given a
//! load line for the mean slip of this one, `mean = r - c e`, with `r` the
//! released mean slip, `c > 0` the compliance of what it is joined to, and
//! `e` the friction force it returns (positive along the slip it resists).
//!
//! Two obvious laws fail on step means. "Stuck means the mean slip is zero"
//! makes the end slip `-s0`: a thrown mass never comes to rest, its velocity
//! flips sign every step and nothing is dissipated. "Stuck means the end
//! slip is zero" alone gives a negative loss `e * mean` in a step in which
//! the drive reverses the body. The law here takes the end-of-step arrest
//! where its loss is not negative, and the mean elsewhere:
//!
//! | state | condition | force `e` | the step ends |
//! |---|---|---|---|
//! | stuck, `s0 = 0` | `abs(r) <= c mu_s N` | `r / c` | stuck |
//! | stuck | otherwise | `mu_k N sign(r)` | sliding |
//! | sliding, `sg = sign(s0)`, `ea = (r - s0/2) / c` | `0 <= sg ea <= mu_k N` | `ea` | at zero slip exactly |
//! | sliding | `sg ea > mu_k N` | `mu_k N sg` | sliding the same way |
//! | sliding | `0 <= sg r < abs(s0)/2` | `0` | reversed, slowly |
//! | sliding | `-c mu_k N <= sg r < 0` | `r / c` | reversed |
//! | sliding | `sg r < -c mu_k N` | `-mu_k N sg` | reversed |
//!
//! - The loss `e * mean` is not negative in any row, and the first law is
//!   exact on step means: what the two ends give up leaves through `heat`.
//! - `abs(e) <= mu_s N` stuck and `<= mu_k N` sliding.
//! - A mass thrown along a surface ends at exactly zero velocity, with all
//!   of its kinetic energy given up as heat.
//! - In the sliding state `e` is continuous and does not decrease with `r`.
//!   It jumps at break-away when `mu_s > mu_k`, as the physics does.
//! - A step in which the slip reverses is first-order accurate. The fifth
//!   row gives the friction of such a step up altogether: the drive alone
//!   reverses the body there, any friction force of the old direction
//!   would have a negative loss on the step's mean, and one of the new
//!   direction is not yet called for by the mean slip.
//!
//! With no compliance (`c = 0`, both flows prescribed) the slip is the
//! released one and the force is `mu_k N sign(r)`.

use fs_exec::Cx;
use fs_ports::PortKind;
use fs_ports::component::{Class, Component, LoadLine, PortRole, PortSpec, Proposal, Receipt};

use super::{ElementError, heat, inputs_finite, port, positive, roundoff, temperature};

/// A dry-friction slider between ports `a` and `b`, with the loss leaving
/// through port `heat`.
///
/// `a` and `b` receive the velocities [m/s] of the two sides and return
/// the friction force [N] and its opposite; the slip is `a - b`. `heat`
/// receives a temperature [K] and returns the heat rate into the slider
/// [W], which is minus the loss. Its state is [`Slip`]. Added with
/// `SystemBuilder::add_on_line`.
#[derive(Debug, Clone)]
pub struct Slider {
    /// Largest force held at zero slip [N].
    hold_n: f64,
    /// Force while sliding [N].
    slide_n: f64,
    ports: [PortSpec; 3],
}

/// The slip a [`Slider`] ended its last step with [m/s]; zero when stuck.
///
/// It is the flow port `a` receives less the flow port `b` receives, with
/// the signs the ports receive them in. A member of a common-flow junction
/// receives MINUS the flow of the port that sets it, so a slider whose `a`
/// is joined to a mass moving at `v`, with `b` held at rest, starts with a
/// slip of `-v`. A slider started with a slip that is not the one its ports
/// carry applies its sliding rows the wrong way round; nothing can check
/// that for it.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Slip {
    /// Flow received at `a` less flow received at `b`, at the end of the
    /// step.
    pub slip_m_s: f64,
}

/// Which row of the law a step took (see the module's table, top to
/// bottom).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SliderRow {
    /// Stuck, and held.
    Held,
    /// Stuck, and broke away.
    BrokeAway,
    /// Sliding, and arrested at zero slip.
    Arrested,
    /// Sliding on the same way.
    Slid,
    /// Sliding, reversed by the drive alone; no friction in the step.
    Released,
    /// Sliding, reversed with zero mean slip.
    Turned,
    /// Sliding, reversed against friction of the new direction.
    Reversed,
}

impl Slider {
    /// A slider pressed with `normal_n` [N], holding up to
    /// `static_coefficient * normal_n` at rest and resisting with
    /// `kinetic_coefficient * normal_n` while it slides.
    ///
    /// # Errors
    /// A normal force or a coefficient that is not positive and finite, or
    /// a kinetic coefficient above the static one.
    pub fn new(
        normal_n: f64,
        static_coefficient: f64,
        kinetic_coefficient: f64,
    ) -> Result<Self, ElementError> {
        let normal_n = positive("normal_n", normal_n)?;
        let holds = positive("static_coefficient", static_coefficient)?;
        let slides = positive("kinetic_coefficient", kinetic_coefficient)?;
        if slides > holds {
            return Err(ElementError::Parameter {
                name: "kinetic_coefficient",
                value: kinetic_coefficient,
            });
        }
        let kind = Some(PortKind::MechanicalForceVelocity);
        Ok(Self {
            hold_n: positive("static friction force", holds * normal_n)?,
            slide_n: positive("kinetic friction force", slides * normal_n)?,
            ports: [
                port("a", kind, PortRole::FlowIn),
                port("b", kind, PortRole::FlowIn),
                heat("heat", PortRole::EffortIn),
            ],
        })
    }

    /// The force, the mean slip and the row of one step: `start` the slip
    /// the step begins with, `released` the released mean slip, `c` the
    /// compliance of the slip to the force.
    ///
    /// A `c` that is not positive is no compliance: both sides are given,
    /// or what the slider is joined to is rigid to within rounding, which
    /// can leave a hair below zero. The slip is then what it is, and the
    /// force is the sliding one. The mean slip returned is still the point
    /// of the line as given, `released - c * force`, so that the loss is
    /// accounted at the flows the ports carry.
    #[must_use]
    pub fn on_line(&self, start: f64, released: f64, c: f64) -> (f64, f64, SliderRow) {
        let (hold, slide) = (self.hold_n, self.slide_n);
        if !(c > 0.0) {
            let force = if released == 0.0 { 0.0 } else { slide * released.signum() };
            let row = if start == 0.0 { SliderRow::BrokeAway } else { SliderRow::Slid };
            return (force, released - c * force, row);
        }
        if start == 0.0 {
            return if released.abs() <= c * hold {
                (released / c, 0.0, SliderRow::Held)
            } else {
                let force = slide * released.signum();
                (force, released - c * force, SliderRow::BrokeAway)
            };
        }
        let way = start.signum();
        let arrest = (released - 0.5 * start) / c;
        let along = way * arrest;
        if along >= 0.0 && along <= slide {
            (arrest, 0.5 * start, SliderRow::Arrested)
        } else if along > slide {
            let force = slide * way;
            (force, released - c * force, SliderRow::Slid)
        } else if way * released >= 0.0 {
            (0.0, released, SliderRow::Released)
        } else if way * released >= -c * slide {
            (released / c, 0.0, SliderRow::Turned)
        } else {
            let force = -slide * way;
            (force, released - c * force, SliderRow::Reversed)
        }
    }
}

impl Component for Slider {
    type State = Slip;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, _slip: &Slip) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.propose_on_line(old, dt_s, inputs, &[0.0; 4], new, outputs, cx)
    }
}

impl LoadLine for Slider {
    fn line_ports(&self) -> &[usize] {
        &[0, 1]
    }

    fn propose_on_line(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        inputs_finite(compliance)?;
        let kelvin = temperature(inputs[2])?;
        // The force acts as (e, -e) on the two ends, so a unit of it takes
        // (1, -1) C (1, -1)^T away from the slip.
        let c = compliance[0] - compliance[1] - compliance[2] + compliance[3];
        let released = inputs[0] - inputs[1];
        let (force, mean, _) = self.on_line(old.slip_m_s, released, c);
        new.slip_m_s = 2.0 * mean - old.slip_m_s;
        let loss_w = force * mean;
        outputs[0] = force;
        outputs[1] = -force;
        outputs[2] = -loss_w;
        let dissipated_j = dt_s * loss_w;
        // The two ends' powers are each as large as their own velocities
        // make them; their sum is rounded at that size, not at the loss's.
        let at_ends = inputs[0].abs()
            + inputs[1].abs()
            + force.abs() * compliance.iter().map(|entry| entry.abs()).sum::<f64>();
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dissipated_j / kelvin,
            energy_tolerance_j: roundoff(dt_s * force.abs() * at_ends + dissipated_j),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}
